"""Reaction-network, energy-accounting, analytic-solution and convergence checks."""

from dataclasses import replace
import csv
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

from burn import BurnConfig, convergence_check, simulate, snapshot
from burn_report import dd_reference_comparison
from fusion import KEV_J, dd_reactivities

ROOT = Path(__file__).resolve().parents[1]


class BranchReactivityTests(unittest.TestCase):
    def test_independent_nrl_dd_reference(self):
        rows = dd_reference_comparison()
        self.assertEqual([r['temperature_kev'] for r in rows], [1, 2, 5, 10, 20, 50, 100])
        for row in rows:
            t = row['temperature_kev']
            tolerance = .40 if t == 1 else .20 if t == 2 else .10
            with self.subTest(temperature=t):
                self.assertLess(abs(row['relative_difference']), tolerance)
        self.assertGreater(rows[0]['relative_difference'], .30)

    def test_branch_ranges_and_positive_values(self):
        for temperature in (.2, 1, 10, 100):
            for value in dd_reactivities(temperature):
                self.assertTrue(math.isfinite(value))
                self.assertGreater(value, 0)
        for bad in (0, -.1, .199, 100.1, math.nan, math.inf):
            with self.subTest(temperature=bad), self.assertRaises(ValueError):
                dd_reactivities(bad)

    def test_identical_reactant_factor_and_product_energies(self):
        config = BurnConfig(volume_m3=3)
        state = (0, 0, 0, 30, 0, 0)
        with patch('burn.dd_reactivities', return_value=(2e-24, 6e-24)):
            row = snapshot(0, state, config)
        self.assertAlmostEqual(row['dd_tritium_reactions_s'] / 3e16, 1)
        self.assertAlmostEqual(row['dd_helium3_reactions_s'] / 9e16, 1)
        self.assertEqual(row['dt_reactions_s'], 0)
        charged = (3e16 * 4030 + 9e16 * 820) * KEV_J
        neutron = 9e16 * 2450 * KEV_J
        self.assertAlmostEqual(row['charged_power_w'] / charged, 1)
        self.assertAlmostEqual(row['neutron_power_w'] / neutron, 1)
        self.assertAlmostEqual(row['fusion_power_w'] / (charged + neutron), 1)


class BurnPhysicsTests(unittest.TestCase):
    def test_pure_deuterium_builds_tritium_then_burns_dt(self):
        result = simulate(BurnConfig(samples=31))
        initial, final = result.rows[0], result.rows[-1]
        self.assertEqual(result.status, 'completed')
        self.assertEqual(initial['tritium_density_m3'], 0)
        self.assertEqual(initial['dt_reactions_s'], 0)
        self.assertGreater(final['tritium_density_m3'], 0)
        self.assertGreater(final['dt_reactions'], 0)
        self.assertGreater(final['helium3_density_m3'], 0)
        self.assertGreater(final['helium4_density_m3'], 0)
        self.assertGreater(final['proton_density_m3'], 0)
        self.assertLess(final['temperature_kev'], initial['temperature_kev'])
        self.assertLess(final['deuterium_density_m3'], initial['deuterium_density_m3'])
        self.assertEqual(initial['net_electric_energy_j'], -initial['startup_electric_energy_j'])

    def test_charge_nucleons_and_stoichiometry(self):
        for fraction in (0, .5):
            config = BurnConfig(tritium_fraction=fraction, volume_m3=2, samples=21)
            result = simulate(config)
            n, v = config.density_m3, config.volume_m3
            for row in result.rows:
                with self.subTest(fraction=fraction, time=row['time_s']):
                    d = row['deuterium_density_m3']
                    t = row['tritium_density_m3']
                    p = row['proton_density_m3']
                    he3, he4 = row['helium3_density_m3'], row['helium4_density_m3']
                    neutron = row['escaped_neutron_density_m3']
                    charge = d + t + p + 2 * he3 + 2 * he4
                    nucleons = 2 * d + 3 * t + p + 3 * he3 + 4 * he4 + neutron
                    self.assertAlmostEqual(charge / n, 1, places=12)
                    self.assertAlmostEqual(nucleons / (n * (2 + fraction)), 1, places=12)
                    self.assertAlmostEqual(row['electron_density_m3'] / charge, 1, places=12)
                    self.assertAlmostEqual((p * v - row['dd_tritium_reactions']) / n, 0, places=12)
                    self.assertAlmostEqual((he3 * v - row['dd_helium3_reactions']) / n, 0, places=12)
                    self.assertAlmostEqual((he4 * v - row['dt_reactions']) / n, 0, places=12)
                    expected_d = n * (1 - fraction) - (row['dt_reactions'] + 2 * row['dd_tritium_reactions']
                                                     + 2 * row['dd_helium3_reactions']) / v
                    self.assertAlmostEqual((d - expected_d) / n, 0, places=12)
                    self.assertTrue(all(row[key] >= 0 for key in row if key.endswith('_density_m3')))

    def test_transient_energy_conservation_and_startup_cost(self):
        config = BurnConfig(tritium_fraction=.2, charged_deposition=.7, heating_mw=.05,
                            other_loss_mw=.01, volume_m3=3, samples=21)
        result = simulate(config)
        initial = result.rows[0]
        for row in result.rows:
            scale = max(initial['stored_energy_j'], row['fusion_energy_j'], row['aux_energy_j'])
            # Whole-system thermal accounting is independent of steady-state assumptions.
            inputs = initial['stored_energy_j'] + row['fusion_energy_j'] + row['aux_energy_j']
            outputs = row['stored_energy_j'] + row['outgoing_heat_energy_j']
            self.assertAlmostEqual((inputs - outputs) / scale, 0, places=11)
            gross_energy = config.thermal_efficiency * config.heat_recovery * row['outgoing_heat_energy_j']
            electrical_use = row['aux_energy_j'] / config.heater_efficiency + config.plant_load_mw * 1e6 * row['time_s']
            self.assertAlmostEqual((row['operating_net_electric_energy_j'] - gross_energy + electrical_use) / scale,
                                   0, places=11)
            self.assertAlmostEqual((row['net_electric_energy_j'] - row['operating_net_electric_energy_j']
                                    + initial['stored_energy_j'] / config.heater_efficiency) / scale, 0, places=11)
        self.assertLess(result.max_energy_relative_error, 1e-11)

    def test_radiation_accounts_for_helium_charge_squared(self):
        config = BurnConfig(tritium_fraction=.5)
        # Physically consistent hypothetical state after DT and both DD branches.
        state = (.1, .03, .02, 20, 0, 0)
        row = snapshot(0, state, config)
        sum_z2_n = (row['deuterium_density_m3'] + row['tritium_density_m3']
                    + row['proton_density_m3'] + 4 * row['helium3_density_m3'] + 4 * row['helium4_density_m3'])
        expected = (1.69e-38 * row['electron_density_m3'] * sum_z2_n
                    * math.sqrt(row['temperature_kev'] * 1000) * config.volume_m3)
        self.assertAlmostEqual(row['bremsstrahlung_loss_w'] / expected, 1, places=12)
        heat_capacity_energy = 1.5 * (row['ion_density_m3'] + row['electron_density_m3']) * row['temperature_kev'] * KEV_J
        self.assertAlmostEqual(row['stored_energy_j'] / (heat_capacity_energy * config.volume_m3), 1)

    def test_volume_scaling_preserves_trajectory_with_extensive_inputs(self):
        base = BurnConfig(heating_mw=.03, other_loss_mw=.01, plant_load_mw=.1, samples=21)
        a = simulate(base).rows[-1]
        b = simulate(replace(base, volume_m3=3, heating_mw=.09, other_loss_mw=.03,
                             plant_load_mw=.3)).rows[-1]
        self.assertAlmostEqual(a['temperature_kev'], b['temperature_kev'], places=10)
        self.assertAlmostEqual(a['tritium_density_m3'] / b['tritium_density_m3'], 1, places=10)
        for key in ('fusion_energy_j', 'stored_energy_j', 'net_electric_energy_j'):
            self.assertAlmostEqual(b[key] / a[key], 3, places=10)

    def test_no_fusion_cooling_cannot_produce_positive_cumulative_net_energy(self):
        config = BurnConfig(tritium_fraction=1, plant_load_mw=0, duration_s=1, samples=11,
                            heater_efficiency=1, heat_recovery=1, thermal_efficiency=1)
        result = simulate(config)
        for row in result.rows:
            self.assertEqual(row['fusion_energy_j'], 0)
            self.assertLess(row['net_electric_energy_j'], 0)
            self.assertAlmostEqual(row['net_electric_energy_j'] / row['stored_energy_j'], -1, places=11)


class NumericalTests(unittest.TestCase):
    def test_analytic_radiative_transport_cooling(self):
        config = BurnConfig(tritium_fraction=1, samples=2, max_step_s=.02, rtol=1e-10)
        result = simulate(config)
        beta = 1.69e-38 * config.density_m3 * math.sqrt(1000 / 3) / KEV_J
        root_e = ((math.sqrt(3 * config.temperature_kev) + beta * config.confinement_s)
                  * math.exp(-config.duration_s / (2 * config.confinement_s)) - beta * config.confinement_s)
        exact_temperature = root_e**2 / 3
        self.assertAlmostEqual(result.rows[-1]['temperature_kev'] / exact_temperature, 1, places=9)

    def test_constant_rate_dt_solution_and_fourth_order_refinement(self):
        config = BurnConfig(tritium_fraction=.5, charged_deposition=0, confinement_s=1e30,
                            duration_s=2, samples=2, rtol=.01, max_step_s=.4)
        # With equal reactants and constant <sigmav>, d(t)=d0/(1+n0*sv*d0*t).
        exact = .5 / (1 + config.density_m3 * 1e-20 * .5 * config.duration_s) * config.density_m3
        with patch('burn.dt_reactivity', return_value=1e-20), patch('burn.dd_reactivities', return_value=(0, 0)):
            coarse = simulate(config)
            fine = simulate(replace(config, max_step_s=.2))
        a = abs(coarse.rows[-1]['deuterium_density_m3'] / exact - 1)
        b = abs(fine.rows[-1]['deuterium_density_m3'] / exact - 1)
        self.assertLess(b, 1e-7)
        self.assertGreater(a / b, 12)
        self.assertLess(a / b, 20)

    def test_constant_rate_dd_solution(self):
        config = BurnConfig(charged_deposition=0, duration_s=2, samples=2, max_step_s=.1, rtol=1e-9)
        sv = 4e-22
        exact = config.density_m3 / (1 + config.density_m3 * sv * config.duration_s)
        with patch('burn.dt_reactivity', return_value=0), patch('burn.dd_reactivities', return_value=(sv, 0)):
            result = simulate(config)
        final = result.rows[-1]
        self.assertAlmostEqual(final['deuterium_density_m3'] / exact, 1, places=9)
        self.assertAlmostEqual(final['tritium_density_m3'] / ((config.density_m3 - exact) / 2), 1, places=8)
        self.assertEqual(final['helium3_density_m3'], 0)

    def test_convergence_check_really_uses_finer_steps(self):
        result, check = convergence_check(BurnConfig(samples=21, max_step_s=.1))
        self.assertTrue(check['passed'])
        self.assertTrue(check['same_termination'])
        self.assertGreater(check['refined_accepted_steps'], result.accepted_steps)
        for comparison in check['comparisons'].values():
            self.assertTrue(comparison['passed'])

    def test_temperature_limits_terminate_without_extrapolation(self):
        for config, status, temperature in (
            (BurnConfig(tritium_fraction=1, duration_s=20, samples=21), 'temperature_lower_bound', .2),
            (BurnConfig(heating_mw=10, duration_s=5, samples=21), 'temperature_upper_bound', 30),
        ):
            with self.subTest(status=status):
                result, check = convergence_check(config)
                self.assertEqual(result.status, status)
                self.assertLess(result.rows[-1]['time_s'], config.duration_s)
                self.assertAlmostEqual(result.rows[-1]['temperature_kev'], temperature, delta=1.1e-8)
                self.assertTrue(check['passed'])
                for row in result.rows:
                    self.assertGreaterEqual(row['temperature_kev'], .2)
                    self.assertLessEqual(row['temperature_kev'], 30)

    def test_initial_bound_terminates_at_zero_time(self):
        result = simulate(BurnConfig(temperature_kev=.2))
        self.assertEqual(result.status, 'temperature_lower_bound')
        self.assertEqual(result.rows[-1]['time_s'], 0)
        self.assertEqual(result.accepted_steps, 0)

    def test_invalid_inputs(self):
        invalid = ({'density_m3': 0}, {'temperature_kev': .1}, {'temperature_kev': 31},
                   {'duration_s': 0}, {'confinement_s': -1}, {'charged_deposition': 1.1},
                   {'heater_efficiency': 0}, {'thermal_efficiency': 2}, {'tritium_fraction': -1},
                   {'max_step_s': 0}, {'rtol': 1e-15}, {'samples': 1}, {'samples': 10002},
                   {'samples': 5.5}, {'heating_mw': math.nan}, {'density_m3': math.inf},
                   {'density_m3': 1e308, 'volume_m3': 1e308})
        for kwargs in invalid:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                BurnConfig(**kwargs)


class BurnCliTests(unittest.TestCase):
    def cli(self, directory, *args):
        return subprocess.run([sys.executable, str(ROOT / 'ignis.py'), 'burn', '--no-browser', *args],
                              cwd=directory, capture_output=True, text=True, timeout=30)

    def test_report_and_exports(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.cli(directory, '--samples', '21', '--output', 'out/report.html',
                              '--csv', 'out/burn.csv', '--json', 'out/burn.json')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('Convergence: PASS', result.stdout)
            payload = json.loads(Path(directory, 'out/burn.json').read_text())
            self.assertEqual(payload['model'], 'closed_volume_dd_dt_burn')
            self.assertEqual(payload['status'], 'completed')
            self.assertEqual(len(payload['rows']), 21)
            self.assertEqual(payload['config']['tritium_fraction'], 0)
            self.assertTrue(payload['convergence']['passed'])
            direct = simulate(BurnConfig(**payload['config']))
            self.assertEqual(payload['rows'][-1], direct.rows[-1])
            with Path(directory, 'out/burn.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 21)
            self.assertEqual(float(rows[0]['initial_duration_s']), 1)
            self.assertEqual(float(rows[-1]['net_electric_energy_j']), payload['rows'][-1]['net_electric_energy_j'])
            report = Path(directory, 'out/report.html').read_text()
            self.assertIn('PROTOTYPE 04', report)
            self.assertIn('including startup', report)
            self.assertIn('Independent D–D reference comparison', report)
            svgs = re.findall(r'<svg\b.*?</svg>', report, re.S)
            self.assertEqual(len(svgs), 5)
            for svg in svgs:
                self.assertIn('role="img"', svg)
                ET.fromstring(svg)
            self.assertNotRegex(Path(directory, 'out/burn.json').read_text(), r'\b(?:NaN|Infinity)\b')

    def test_explicit_temperature_heating_and_boundary_export(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.cli(directory, '--temperature-k', '1e8', '--heating-mw', '10',
                              '--duration-s', '5', '--samples', '11', '--json', 'burn.json')
            self.assertEqual(result.returncode, 0, result.stderr)
            payload = json.loads(Path(directory, 'burn.json').read_text())
            self.assertEqual(payload['status'], 'temperature_upper_bound')
            self.assertLess(payload['rows'][-1]['time_s'], 5)
            self.assertAlmostEqual(payload['config']['temperature_kev'], 8.61733326214518, places=12)

    def test_invalid_cli_writes_nothing(self):
        for args in (('--temperature-k', '0'), ('--duration-s', '0'), ('--heating-mw', '-1'),
                     ('--tritium-fraction', '2'), ('--density-m3', 'nan'), ('--samples', '1'),
                     ('--output', 'same', '--csv', './same'), ('--temperature-k', '1e9')):
            with self.subTest(args=args), tempfile.TemporaryDirectory() as directory:
                result = self.cli(directory, *args)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn('Traceback', result.stderr)
                self.assertEqual(list(Path(directory).iterdir()), [])

    def test_output_error_is_readable(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, 'blocker').write_text('existing')
            result = self.cli(directory, '--samples', '2', '--duration-s', '.01', '--output', 'blocker/report.html')
            self.assertEqual(result.returncode, 1)
            self.assertIn('Cannot write results', result.stderr)
            self.assertNotIn('Traceback', result.stderr)


if __name__ == '__main__':
    unittest.main()
