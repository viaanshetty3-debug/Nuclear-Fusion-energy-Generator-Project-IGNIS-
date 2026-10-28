"""Independent reference checks, conservation laws, and end-to-end model exports."""

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
from unittest.mock import patch
import xml.etree.ElementTree as ET

from fusion import (ALPHA_ENERGY_J, NEUTRON_ENERGY_J, PlantConfig, dt_reactivity,
                    evaluate_case)
from fusion_report import reference_comparison, render_fusion_sections

ROOT = Path(__file__).resolve().parents[1]


class ReactivityTests(unittest.TestCase):
    def test_independent_nrl_benchmarks(self):
        comparisons = reference_comparison()
        self.assertEqual([r['temperature_kev'] for r in comparisons], [1, 2, 5, 10, 20, 50, 100])
        for row in comparisons:
            # NRL is an independent, older, rounded table, not Bosch-Hale output.
            # Low-temperature differences are explicitly published in our report.
            t = row['temperature_kev']
            tolerance = 0.30 if t == 1 else 0.20 if t == 2 else 0.06
            with self.subTest(temperature=t):
                self.assertLess(abs(row['relative_difference']), tolerance)
        self.assertGreater(abs(comparisons[0]['relative_difference']), 0.20)

    def test_fit_range_and_finite_values(self):
        for t in (0.2, 1, 10, 100):
            self.assertTrue(math.isfinite(dt_reactivity(t)))
            self.assertGreater(dt_reactivity(t), 0)
        for t in (0, -1, 0.1999, 100.0001, math.inf, math.nan):
            with self.subTest(t=t), self.assertRaises(ValueError):
                dt_reactivity(t)


class PowerBalanceTests(unittest.TestCase):
    def test_known_ten_kev_case(self):
        case = evaluate_case(10, 1e20)
        # NRL 10 keV <sigmav> ~1.1e-22 m³/s -> R~2.75e17/s for a 1 m³ 50/50 plasma.
        self.assertAlmostEqual(case.reactions_s / 2.75e17, 1, delta=0.04)
        self.assertAlmostEqual(case.fusion_power_w / 775454.0, 1, delta=0.04)
        # Independent SI checks: U=3*n*kT*V, and NRL Eq30 in its original cgs units.
        self.assertAlmostEqual(case.stored_energy_j, 480652.9902, places=6)
        cgs_brem_w_cm3 = 1.69e-32 * (1e14)**2 * math.sqrt(10000)
        self.assertAlmostEqual(case.bremsstrahlung_loss_w, cgs_brem_w_cm3 * 1e6, places=8)
        self.assertLess(case.net_electric_w, 0)
        self.assertGreater(case.plasma_gain_q, 1)  # Q>1 alone does not imply net electricity.

    def test_distinct_reactants_and_reaction_energy(self):
        with patch('fusion.dt_reactivity', return_value=4e-22):
            case = evaluate_case(10, 1e20, PlantConfig(volume_m3=2))
        self.assertAlmostEqual(case.reaction_rate_m3_s / 1e18, 1)
        self.assertAlmostEqual(case.reactions_s / 2e18, 1)
        self.assertAlmostEqual(case.alpha_power_w / (case.reactions_s * ALPHA_ENERGY_J), 1)
        self.assertAlmostEqual(case.neutron_power_w / (case.reactions_s * NEUTRON_ENERGY_J), 1)
        self.assertAlmostEqual(case.fusion_power_w / (case.alpha_power_w + case.neutron_power_w), 1)

    def test_density_volume_and_transport_scaling(self):
        config = PlantConfig()
        base = evaluate_case(10, 1e20, config)
        dense = evaluate_case(10, 2e20, config)
        larger = evaluate_case(10, 1e20, replace(config, volume_m3=3))
        longer = evaluate_case(10, 1e20, replace(config, confinement_s=2))
        for field in ('reactions_s', 'fusion_power_w', 'bremsstrahlung_loss_w'):
            self.assertAlmostEqual(getattr(dense, field) / getattr(base, field), 4)
        self.assertAlmostEqual(dense.transport_loss_w / base.transport_loss_w, 2)
        for field in ('reactions_s', 'fusion_power_w', 'stored_energy_j',
                      'transport_loss_w', 'bremsstrahlung_loss_w', 'required_heating_w'):
            self.assertAlmostEqual(getattr(larger, field) / getattr(base, field), 3)
        self.assertAlmostEqual(longer.transport_loss_w / base.transport_loss_w, 0.5)
        self.assertEqual(longer.fusion_power_w, base.fusion_power_w)
        self.assertLess(longer.required_heating_w, base.required_heating_w)

    def test_mixture_symmetry_maximum_and_pure_endpoints(self):
        balanced = evaluate_case(10, 1e20)
        a = evaluate_case(10, 1e20, PlantConfig(tritium_fraction=0.2))
        b = evaluate_case(10, 1e20, PlantConfig(tritium_fraction=0.8))
        self.assertAlmostEqual(a.fusion_power_w / b.fusion_power_w, 1)
        self.assertGreater(balanced.fusion_power_w, a.fusion_power_w)
        for fraction in (0, 1):
            case = evaluate_case(10, 1e20, PlantConfig(tritium_fraction=fraction))
            self.assertEqual(case.reactions_s, 0)
            self.assertEqual(case.plasma_gain_q, 0)
            self.assertEqual(case.required_heating_w, case.total_loss_w)

    def test_plasma_and_plant_energy_closure(self):
        for density in (1e19, 1e20, 1e21):
            for deposition in (0, 0.5, 1):
                with self.subTest(density=density, deposition=deposition):
                    config = PlantConfig(alpha_deposition=deposition, other_loss_mw=0.01)
                    c = evaluate_case(20, density, config)
                    plasma_input = c.alpha_heating_w + c.required_heating_w
                    plasma_output = c.total_loss_w + c.extra_exhaust_w
                    self.assertAlmostEqual(plasma_input / plasma_output, 1, places=12)
                    cycle_source = (c.neutron_power_w + c.alpha_power_w - c.alpha_heating_w
                                    + c.total_loss_w + c.extra_exhaust_w)
                    self.assertAlmostEqual(c.recovered_heat_w / (config.heat_recovery * cycle_source), 1)
                    self.assertAlmostEqual(c.recovered_heat_w / (
                        config.heat_recovery * (c.fusion_power_w + c.required_heating_w)), 1)
                    self.assertAlmostEqual(c.net_electric_w,
                                           c.gross_electric_w - c.heater_electric_w - c.plant_load_w)
                    self.assertAlmostEqual(c.heater_electric_w,
                                           c.required_heating_w / config.heater_efficiency)

    def test_self_heating_surplus_requires_exhaust(self):
        c = evaluate_case(20, 1e21)
        self.assertEqual(c.required_heating_w, 0)
        self.assertIsNone(c.plasma_gain_q)
        self.assertGreater(c.extra_exhaust_w, 0)
        self.assertEqual(c.extra_exhaust_w, c.self_heating_margin_w)
        self.assertGreater(c.net_electric_w, 0)

    def test_zero_recovery_and_increased_electrical_load(self):
        zero = evaluate_case(10, 1e20, PlantConfig(heat_recovery=0))
        self.assertEqual(zero.gross_electric_w, 0)
        self.assertEqual(zero.net_electric_w, -zero.heater_electric_w - zero.plant_load_w)
        base = evaluate_case(10, 1e20)
        load = evaluate_case(10, 1e20, PlantConfig(plant_load_mw=2))
        self.assertAlmostEqual(load.net_electric_w - base.net_electric_w, -1e6)
        bad_heater = evaluate_case(10, 1e20, PlantConfig(heater_efficiency=0.25))
        self.assertAlmostEqual(bad_heater.heater_electric_w / base.heater_electric_w, 2)
        self.assertLess(bad_heater.net_electric_w, base.net_electric_w)

    def test_additional_losses_increase_heating(self):
        base = evaluate_case(10, 1e20)
        lossy = evaluate_case(10, 1e20, PlantConfig(other_loss_mw=1))
        self.assertAlmostEqual(lossy.required_heating_w - base.required_heating_w, 1e6)
        self.assertLess(lossy.net_electric_w, base.net_electric_w)

    def test_invalid_configurations(self):
        for name in asdict(PlantConfig()):
            for value in (-1, math.nan, math.inf):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    PlantConfig(**{name: value})
        for name in ('tritium_fraction', 'alpha_deposition', 'heater_efficiency',
                     'thermal_efficiency', 'heat_recovery'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                PlantConfig(**{name: 1.01})
        for name in ('volume_m3', 'confinement_s', 'heater_efficiency'):
            with self.subTest(name=name), self.assertRaises(ValueError):
                PlantConfig(**{name: 0})
        for temperature, density in ((0.1, 1e20), (31, 1e20), (10, 0), (10, math.nan),
                                     (10, 1e308), (10, 1e-300)):
            with self.subTest(t=temperature, n=density), self.assertRaises(ValueError):
                evaluate_case(temperature, density)

    def test_plot_boundaries_and_zero_yield(self):
        config = PlantConfig(tritium_fraction=0)
        cases = [evaluate_case(t, 1e20, config) for t in (0.2, 30)]
        report = render_fusion_sections(cases, config)
        self.assertIn('Gross fusion power', report)
        for svg in re.findall(r'<svg\b.*?</svg>', report, flags=re.S):
            ET.fromstring(svg)
            self.assertNotRegex(svg, r'(?:nan|inf),')


class FusionCliTests(unittest.TestCase):
    def run_cli(self, directory, *args):
        return subprocess.run([sys.executable, str(ROOT / 'ignis.py'), '--no-browser', *args],
                              cwd=directory, capture_output=True, text=True, timeout=20)

    def test_cartesian_sweep_and_reproducible_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_cli(directory, '--temperatures', '1e8', '2e8', '--densities',
                                  '1e20', '2e20', '--volume-m3', '10', '--confinement-s', '2',
                                  '--plant-load-mw', '0.1', '--fusion-csv', 'out/cases.csv',
                                  '--json', 'out/scenario.json', '--output', 'out/report.html')
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(Path(directory, 'out/scenario.json').read_text())
            self.assertEqual(len(payload['cases']), 4)
            self.assertEqual(payload['config']['volume_m3'], 10)
            self.assertEqual(payload['config']['confinement_s'], 2)
            for row in payload['cases']:
                direct = evaluate_case(row['temperature_kev'], row['density_m3'],
                                       PlantConfig(**payload['config'])).to_dict()
                self.assertEqual(row, direct)
            with Path(directory, 'out/cases.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 4)
            self.assertEqual(float(rows[0]['volume_m3']), 10)
            self.assertAlmostEqual(float(rows[0]['net_electric_w']), payload['cases'][0]['net_electric_w'])
            report = Path(directory, 'out/report.html').read_text()
            self.assertIn('PROTOTYPE 03', report)
            self.assertIn('Net electricity (MW)', report)
            self.assertIn('Published benchmark comparison', report)
            self.assertIn('24.67%', report)
            self.assertIn('reactor feasibility', report)
            svgs = re.findall(r'<svg\b.*?</svg>', report, re.S)
            self.assertEqual(len(svgs), 3)
            for svg in svgs:
                self.assertIn('role="img"', svg)
                ET.fromstring(svg)

    def test_strict_json_handles_zero_auxiliary_heating(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_cli(directory, '--temperatures', '2e8', '--densities', '1e21',
                                  '--fusion-csv', 'cases.csv', '--json', 'scenario.json')
            self.assertEqual(result.returncode, 0, result.stderr)
            raw = Path(directory, 'scenario.json').read_text()
            self.assertNotIn('Infinity', raw)
            payload = json.loads(raw)
            self.assertIsNone(payload['cases'][0]['plasma_gain_q'])
            self.assertEqual(payload['cases'][0]['required_heating_w'], 0)
            with Path(directory, 'cases.csv').open(newline='') as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row['plasma_gain_q'], '')

    def test_invalid_full_model_inputs_write_no_files(self):
        invalid = (
            ('--temperatures', '1e9'), ('--temperatures', '1000'),
            ('--densities', '0'), ('--densities', 'nan'), ('--densities', '1e308'),
            ('--confinement-s', '0'), ('--heater-efficiency', '0'),
            ('--thermal-efficiency', '1.1'), ('--tritium-fraction', '2'),
            ('--volume-m3', '0'), ('--plant-load-mw', '-1'), ('--alpha-deposition', 'nan'),
            ('--fusion-csv', 'same', '--json', './same'),
            ('--thermal-only', '--json', 'scenario.json'),
            ('--thermal-only', '--densities', '1e20'),
        )
        for args in invalid:
            with self.subTest(args=args), tempfile.TemporaryDirectory() as directory:
                result = self.run_cli(directory, *args)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                self.assertEqual(list(Path(directory).iterdir()), [])

    def test_thermal_only_retains_broader_temperature_domain(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_cli(directory, '--thermal-only', '--temperatures', '1000')
            self.assertEqual(result.returncode, 0, result.stderr)
            report = Path(directory, 'ignis_results.html').read_text()
            self.assertNotIn('id="fusion-model"', report)
            self.assertNotIn('D-T target states', result.stdout)


if __name__ == '__main__':
    unittest.main()
