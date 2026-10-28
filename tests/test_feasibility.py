"""Published scaling fixtures, coupled power closure, screening and sensitivity."""

import csv
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

from confinement import ReferenceTokamak, greenwald_density, ipb98y2_time, martin_threshold_mw
from feasibility import (StudyConfig, StudyGrid, break_even_roots, evaluate_point,
                         one_at_a_time, run_study, sensitivity_corners)

ROOT = Path(__file__).resolve().parents[1]


class PublishedScalingTests(unittest.TestCase):
    def test_ukaea_ipb98_unit_fixture(self):
        # Official PROCESS test: all fit-unit inputs equal 1 -> 0.0562 seconds.
        self.assertAlmostEqual(ipb98y2_time(1, 1, 1e19, 1, 1, 1, 1, 1), .0562, places=12)

    def test_ukaea_martin_nominal_fixture(self):
        # PROCESS test_calculate_martin08_nominal: n20=1, B=5 T, S=100 m², M=2.
        self.assertAlmostEqual(martin_threshold_mw(1e20, 5, 100, 2), 13.542309512892546, places=11)
        self.assertAlmostEqual(martin_threshold_mw(1e20, 5, 100, 4), 13.542309512892546 / 2, places=11)

    def test_power_density_and_h_scaling(self):
        args = [15, 5.3, 1e20, 100, 6.2, 1.7, 2 / 6.2, 2.5]
        base = ipb98y2_time(*args)
        doubled = args.copy()
        doubled[3] *= 2
        self.assertAlmostEqual(ipb98y2_time(*doubled) / base, 2**(-.69))
        doubled = args.copy()
        doubled[2] *= 10
        self.assertAlmostEqual(ipb98y2_time(*doubled) / base, 10**.41)
        self.assertAlmostEqual(ipb98y2_time(*args, h_factor=1.2) / base, 1.2)

    def test_greenwald_units_and_circular_geometry(self):
        self.assertAlmostEqual(greenwald_density(15, 2) / (15e20 / (4 * math.pi)), 1)
        device = ReferenceTokamak(elongation=1)
        self.assertAlmostEqual(device.volume_m3 / (2 * math.pi**2 * 6.2 * 4), 1)
        self.assertAlmostEqual(device.surface_m2 / (4 * math.pi**2 * 6.2 * 2), 1)

    def test_invalid_reference_inputs(self):
        for kwargs in ({'field_t': 0}, {'current_ma': math.nan}, {'minor_radius_m': 0},
                       {'major_radius_m': 5, 'minor_radius_m': 2}, {'elongation': .9}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                ReferenceTokamak(**kwargs)
        for power in (0, -1, math.inf, math.nan):
            with self.subTest(power=power), self.assertRaises(ValueError):
                ipb98y2_time(15, 5.3, 1e20, power, 6.2, 1.7, 2 / 6.2, 2.5)


class CoupledBalanceTests(unittest.TestCase):
    def test_scaling_and_plasma_power_close(self):
        for convention in ('transport', 'total'):
            config = StudyConfig(radiation_convention=convention)
            p = evaluate_point(10, 1e20, 1.1, .5, config)
            self.assertTrue(p['passes_screens'], p['exclusion_reasons'])
            self.assertAlmostEqual(p['scaling_confinement_s'] * p['scaled_loss_power_w'] / p['stored_energy_j'], 1, places=12)
            self.assertAlmostEqual((p['alpha_heating_w'] + p['required_heating_w']) / p['total_loss_w'], 1, places=12)
            self.assertAlmostEqual((p['transport_loss_w'] + p['radiation_loss_w']) / p['total_loss_w'], 1, places=12)
            self.assertAlmostEqual(p['transport_confinement_s'] * p['transport_loss_w'] / p['stored_energy_j'], 1, places=12)

    def test_radiation_is_not_counted_twice(self):
        a = evaluate_point(10, 1e20, 1.1, .5, StudyConfig(radiation_convention='transport'))
        b = evaluate_point(10, 1e20, 1.1, .5, StudyConfig(radiation_convention='total'))
        self.assertEqual(a['scaled_loss_power_w'], b['scaled_loss_power_w'])
        self.assertAlmostEqual((a['total_loss_w'] - b['total_loss_w']) / a['radiation_loss_w'], 1)
        self.assertAlmostEqual((a['transport_loss_w'] - b['transport_loss_w']) / a['radiation_loss_w'], 1)
        self.assertGreater(b['transport_confinement_s'], b['scaling_confinement_s'])

    def test_operating_power_is_distinct_from_pulse_energy(self):
        config = StudyConfig()
        point = evaluate_point(10, 1e20, 1, .5, config)
        self.assertTrue(point['passes_screens'])
        self.assertGreater(point['net_electric_w'], 0)
        self.assertLess(point['net_electric_energy_j'], 0)
        self.assertGreater(point['startup_payback_s'], config.pulse_duration_s)

    def test_startup_duration_and_plant_energy_accounting(self):
        config = StudyConfig()
        point = evaluate_point(10, 1e20, 1.1, .5, config)
        gross = config.thermal_efficiency * config.heat_recovery * (point['fusion_power_w'] + point['required_heating_w'])
        demand = point['required_heating_w'] / config.heater_efficiency + config.plant_load_mw * 1e6
        self.assertAlmostEqual((gross - demand) / point['net_electric_w'], 1)
        self.assertAlmostEqual(point['gross_electric_w'] / gross, 1)
        self.assertAlmostEqual((point['gross_electric_w'] - point['heater_electric_w'] - point['plant_load_w'])
                               / point['net_electric_w'], 1)
        self.assertAlmostEqual((point['net_electric_w'] * config.pulse_duration_s - point['net_electric_energy_j'])
                               / point['startup_electric_energy_j'], 1)
        twice = evaluate_point(10, 1e20, 1.1, .5, replace(config, startup_multiplier=2))
        self.assertEqual(twice['net_electric_w'], point['net_electric_w'])
        self.assertAlmostEqual((point['net_electric_energy_j'] - twice['net_electric_energy_j']) / point['startup_electric_energy_j'], 1)
        longer = evaluate_point(10, 1e20, 1.1, .5, replace(config, pulse_duration_s=20))
        self.assertAlmostEqual((longer['net_electric_energy_j'] - point['net_electric_energy_j'])
                               / (point['net_electric_w'] * 10), 1)

    def test_self_heating_surplus_is_excluded(self):
        point = evaluate_point(10, 1e20, 2, .5)
        self.assertLess(point['required_heating_w'], 0)
        self.assertIn('self_heating_surplus_no_fixed_point', point['exclusion_reasons'])
        self.assertIsNone(point['net_electric_w'])
        self.assertIsNone(point['net_electric_energy_j'])

    def test_density_pressure_hmode_capacity_and_radiation_screens(self):
        cases = (
            (evaluate_point(10, 1.5e20, 1, .5), 'greenwald_density'),
            (evaluate_point(20, 1e20, 1, .5), 'thermal_beta'),
            (evaluate_point(5, 1e20, 1, .5), 'hmode_threshold'),
            (evaluate_point(10, 1e20, 1.1, .5, StudyConfig(heating_capacity_mw=0)), 'heating_capacity'),
            (evaluate_point(10, 1e20, 1.1, .5, StudyConfig(radiation_fraction_limit=0)), 'radiation_fraction'),
            (evaluate_point(10, 1e20, 1.1, .5, StudyConfig(radiation_convention='total', extra_radiation_mw=1e5)),
             'radiation_exceeds_total_loss'),
        )
        for point, reason in cases:
            with self.subTest(reason=reason):
                self.assertIn(reason, point['exclusion_reasons'])
                self.assertFalse(point['passes_screens'])
                self.assertIsNone(point['net_electric_energy_j'])

    def test_no_dt_fuel_cannot_produce_positive_electricity(self):
        point = evaluate_point(10, 1e20, 1, 0, StudyConfig(heating_capacity_mw=1e6))
        self.assertTrue(point['passes_screens'])
        self.assertEqual(point['fusion_power_w'], 0)
        self.assertLess(point['net_electric_w'], 0)

    def test_stronger_h_does_not_leave_heating_fixed(self):
        a = evaluate_point(10, 1e20, 1, .5)
        b = evaluate_point(10, 1e20, 1.1, .5)
        self.assertLess(b['transport_loss_w'], a['transport_loss_w'])
        self.assertLess(b['required_heating_w'], a['required_heating_w'])
        self.assertAlmostEqual(b['transport_loss_w'] / a['transport_loss_w'], (1.1)**(-1 / .31), places=12)


class StudyAndSensitivityTests(unittest.TestCase):
    def test_cartesian_grid_and_filtering(self):
        grid = StudyGrid((5, 10), (1e20,), (1, 1.1), (.5,))
        study = run_study(grid=grid)
        self.assertEqual(len(study['rows']), 4)
        eligible = [r for r in study['rows'] if r['passes_screens']]
        self.assertEqual(study['summary']['eligible_points'], len(eligible))
        self.assertTrue(study['best_eligible_point']['passes_screens'])
        self.assertEqual(study['best_eligible_point']['net_electric_energy_j'], max(r['net_electric_energy_j'] for r in eligible))
        self.assertTrue(all(r['sensitivity'] is None for r in study['rows'] if not r['passes_screens']))

    def test_joint_corners_report_invalid_states(self):
        config = StudyConfig()
        baseline = evaluate_point(12.5, 1.15e20, 1.2, .5, config)
        self.assertTrue(baseline['passes_screens'])
        corners = sensitivity_corners(baseline, config)
        self.assertEqual(corners['scenario_count'], 128)
        self.assertFalse(corners['all_corners_positive'])
        self.assertLess(corners['eligible_scenarios'], corners['scenario_count'])
        self.assertIn('self_heating_surplus_no_fixed_point', corners['corner_exclusion_counts'])
        one_at_time = one_at_a_time(baseline, config)
        self.assertEqual(len(one_at_time), 13)
        self.assertTrue(any(r['exclusion_reasons'] for r in one_at_time))

    def test_all_corner_positive_can_pass_with_zero_excursion(self):
        config = StudyConfig(sensitivity_fraction=0)
        baseline = evaluate_point(10, 1e20, 1.1, .5, config)
        corners = sensitivity_corners(baseline, config)
        self.assertTrue(corners['all_corners_positive'])
        self.assertEqual(corners['eligible_scenarios'], 128)
        self.assertEqual(corners['positive_energy_scenarios'], 128)

    def test_break_even_refinement_and_excluded_gaps(self):
        config = StudyConfig()
        left = evaluate_point(10, 1.15e20, 1.2, .25, config)
        right = evaluate_point(12.5, 1.15e20, 1.2, .25, config)
        roots = break_even_roots([left, right], config)
        self.assertEqual(len(roots), 1)
        self.assertGreater(roots[0]['temperature_kev'], 10)
        self.assertLess(roots[0]['temperature_kev'], 12.5)
        self.assertLess(abs(roots[0]['net_electric_energy_j']), 1)
        gap = evaluate_point(11, 1.15e20, 1.2, .25, config)
        gap = {**gap, 'passes_screens': False, 'net_electric_energy_j': None}
        self.assertEqual(break_even_roots([left, gap, right], config), [])

    def test_empty_eligible_window_is_a_valid_result(self):
        study = run_study(grid=StudyGrid((5,), (1e20,), (1,), (.5,)))
        self.assertEqual(study['summary']['eligible_points'], 0)
        self.assertIsNone(study['best_eligible_point'])
        self.assertEqual(study['best_point_sensitivity'], [])
        self.assertEqual(study['break_even_roots'], [])

    def test_config_and_grid_validation(self):
        for kwargs in ({'heater_efficiency': 0}, {'thermal_efficiency': 2}, {'pulse_duration_s': 0},
                       {'startup_multiplier': .5}, {'sensitivity_fraction': .6}, {'beta_n_limit': 0},
                       {'hmode_margin': .9}, {'radiation_convention': 'unknown'}, {'plant_load_mw': math.nan}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                StudyConfig(**kwargs)
        for kwargs in ({'temperatures_kev': ()}, {'temperatures_kev': (31,)}, {'densities_m3': (0,)},
                       {'h_factors': (0,)}, {'tritium_fractions': (2,)}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                StudyGrid(**kwargs)
        config = StudyConfig(**asdict(StudyConfig()))
        self.assertIsInstance(config.device, ReferenceTokamak)
        self.assertEqual(StudyGrid((10, 10), (1e20,), (1,), (.5,)).temperatures_kev, (10,))


class StudyCliTests(unittest.TestCase):
    def cli(self, directory, *args):
        return subprocess.run([sys.executable, str(ROOT / 'ignis.py'), 'study', '--no-browser',
                               '--temperatures-kev', '5', '10', '--densities', '1e20',
                               '--h-factors', '1.1', '--tritium-fractions', '.5', *args],
                              cwd=directory, capture_output=True, text=True, timeout=30)

    def test_reports_filters_and_reproducible_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.cli(directory, '--output', 'out/report.html', '--csv', 'out/study.csv', '--json', 'out/study.json')
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(Path(directory, 'out/study.json').read_text())
            self.assertEqual(len(payload['rows']), 2)
            reproduced = run_study(StudyConfig(**payload['config']), StudyGrid(**payload['grid']))
            self.assertEqual(reproduced['summary'], payload['summary'])
            self.assertEqual(reproduced['rows'], payload['rows'])
            report = Path(directory, 'out/report.html').read_text()
            self.assertIn('PROTOTYPE 05', report)
            self.assertIn('id="h-select"', report)
            self.assertIn('id="f-select"', report)
            self.assertIn('hSelect.addEventListener', report)
            self.assertIn('not probabilities', report)
            self.assertIn('Highest-ranked eligible samples', report)
            svgs = re.findall(r'<svg\b.*?</svg>', report, re.S)
            self.assertEqual(len(svgs), 1)
            ET.fromstring(svgs[0])
            with Path(directory, 'out/study.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2)
            self.assertEqual(float(rows[0]['input_current_ma']), 15)
            self.assertEqual(rows[0]['net_electric_energy_j'], '')
            self.assertIn('hmode_threshold', rows[0]['exclusion_reasons'])
            self.assertNotRegex(Path(directory, 'out/study.json').read_text(), r'\b(?:NaN|Infinity)\b')

    def test_invalid_cli_writes_nothing(self):
        for args in (('--temperatures-kev', '31'), ('--densities', '0'), ('--h-factors', '0'),
                     ('--heater-efficiency', '0'), ('--major-radius-m', '5'),
                     ('--output', 'same', '--json', './same'), ('--pulse-duration-s', '0')):
            with self.subTest(args=args), tempfile.TemporaryDirectory() as directory:
                result = self.cli(directory, *args)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                self.assertEqual(list(Path(directory).iterdir()), [])


if __name__ == '__main__':
    unittest.main()
