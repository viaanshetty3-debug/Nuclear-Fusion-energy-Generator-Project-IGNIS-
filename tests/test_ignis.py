"""Physics invariants and public CLI behavior; no external test dependencies."""

import contextlib
import csv
import io
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import ignis

SCRIPT = Path(ignis.__file__).resolve()


def integrate(function, start, stop, intervals=4000):
    """Composite Simpson quadrature, independent of the plot sampling grid."""
    step = (stop - start) / intervals
    return step / 3 * math.fsum(
        (1 if i in (0, intervals) else 4 if i % 2 else 2)
        * function(start + step * i) for i in range(intervals + 1)
    )


class PhysicsTests(unittest.TestCase):
    def test_si_conversion(self):
        self.assertAlmostEqual(ignis.thermal_energy(100_000_000), 8.61733326214518, places=12)

    def test_normalization_and_mean(self):
        # E = kT u² eliminates the square-root endpoint singularity.
        for temperature in (1e4, 5e7, 1e8, 2e8):
            with self.subTest(temperature=temperature):
                scale = ignis.thermal_energy(temperature)
                def density(u):
                    return ignis.probability(scale * u * u, temperature) * 2 * scale * u
                area = integrate(density, 0, math.sqrt(50))
                mean = integrate(lambda u: scale * u * u * density(u), 0, math.sqrt(50))
                self.assertAlmostEqual(area, 1.0, places=10)
                self.assertAlmostEqual(mean / scale, 1.5, places=10)

    def test_energy_density_mode(self):
        for temperature in ignis.TEMPERATURES:
            scale = ignis.thermal_energy(temperature)
            peak = ignis.probability(scale / 2, temperature)
            self.assertGreater(peak, ignis.probability(scale * 0.49, temperature))
            self.assertGreater(peak, ignis.probability(scale * 0.51, temperature))
            self.assertEqual(ignis.probability(0, temperature), 0)

    def test_temperature_scaling(self):
        for energy in (0.1, 10, 50):
            self.assertAlmostEqual(
                2 * ignis.probability(2 * energy, 2e8),
                ignis.probability(energy, 1e8), places=12,
            )

    def test_tail_matches_independent_quadrature(self):
        temperature = 1e8
        scale = ignis.thermal_energy(temperature)
        for x in (0, 0.01, 0.5, 1, 5, 20):
            with self.subTest(x=x):
                numeric = integrate(
                    lambda u: ignis.probability(scale * u * u, temperature) * 2 * scale * u,
                    math.sqrt(x), math.sqrt(50),
                )
                self.assertAlmostEqual(ignis.fraction_above(x * scale, temperature), numeric, places=10)
        self.assertEqual(ignis.fraction_above(0, temperature), 1)
        self.assertEqual(ignis.fraction_above(1e308, temperature), 0)
        self.assertEqual(ignis.probability(1e308, temperature), 0)

    def test_plot_grid_resolves_separated_temperatures(self):
        for temperature, energies, values in ignis.sample_curves((1e3, 1e9)):
            area = math.fsum((b - a) * (p + q) / 2
                             for a, b, p, q in zip(energies, energies[1:], values, values[1:]))
            self.assertAlmostEqual(area, 1, places=5)
            self.assertAlmostEqual(energies[-1] / ignis.thermal_energy(temperature), 20)
            self.assertLess(ignis.fraction_above(energies[-1], temperature), 2e-8)

    def test_invalid_physics_inputs(self):
        for bad in (-1, math.nan, math.inf, -math.inf):
            with self.subTest(value=bad):
                with self.assertRaises(ValueError):
                    ignis.probability(bad, 1e8)
                with self.assertRaises(ValueError):
                    ignis.fraction_above(bad, 1e8)
        for bad in (0, -1, math.nan, math.inf):
            with self.subTest(temperature=bad):
                with self.assertRaises(ValueError):
                    ignis.thermal_energy(bad)
        for samples in (0, 50, 100002, 51.5):
            with self.assertRaises(ValueError):
                ignis.sample_curves((1e8,), samples)
        with self.assertRaises(ValueError):
            ignis.sample_curves(())


class CommandLineTests(unittest.TestCase):
    def run_cli(self, directory, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *args], cwd=directory,
                              capture_output=True, text=True, timeout=20)

    def test_default_report(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_cli(directory, "--no-browser")
            self.assertEqual(result.returncode, 0, result.stderr)
            report = (Path(directory) / "ignis_results.html").read_text(encoding="utf-8")
            self.assertIn("Mean energy (keV)", report)
            self.assertIn("6.463", report)
            self.assertIn("12.926", report)
            svg = report[report.index("<svg"):report.index("</svg>") + len("</svg>")]
            root = ET.fromstring(svg)
            self.assertEqual(len(root.findall("polyline")), 3)
            self.assertIn("role=\"img\"", svg)

    def test_custom_report_and_csv(self):
        with tempfile.TemporaryDirectory() as directory:
            result = self.run_cli(directory, "--temperatures", "1e8", "2e8", "--samples", "101",
                                  "--threshold-kev", "0", "--no-browser",
                                  "--output", "nested/report.html", "--csv", "nested/curves.csv")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("100.000000%", result.stdout)
            path = Path(directory) / "nested/curves.csv"
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 202)
            self.assertEqual(set(rows[0]), {"temperature_K", "energy_keV", "probability_density_per_keV"})
            self.assertEqual(float(rows[0]["energy_keV"]), 0)
            self.assertEqual(float(rows[101]["temperature_K"]), 2e8)
            self.assertTrue((Path(directory) / "nested/report.html").is_file())

    def test_invalid_cli_writes_no_files(self):
        invalid = (
            ("--temperatures", "0"), ("--temperatures", "nan"),
            ("--temperatures", "inf"), ("--temperatures", "1e-320"),
            ("--temperatures", "1e309"),
            ("--threshold-kev", "-1"), ("--threshold-kev", "nan"),
            ("--samples", "50"), ("--samples", "100002"),
            ("--output", "same", "--csv", "./same"),
        )
        for args in invalid:
            with self.subTest(args=args), tempfile.TemporaryDirectory() as directory:
                result = self.run_cli(directory, "--no-browser", *args)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertNotIn("Traceback", result.stderr)
                self.assertEqual(list(Path(directory).iterdir()), [])

    def test_output_error_is_readable(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "blocker").write_text("existing", encoding="utf-8")
            result = self.run_cli(directory, "--no-browser", "--output", "blocker/report.html")
            self.assertEqual(result.returncode, 1)
            self.assertIn("Cannot write results", result.stderr)
            self.assertNotIn("Traceback", result.stderr)

    def test_browser_failure_does_not_discard_report(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.html"
            with patch("ignis.webbrowser.open", return_value=False) as browser, \
                    contextlib.redirect_stdout(io.StringIO()) as output:
                self.assertEqual(ignis.main(["--output", str(path)]), 0)
            browser.assert_called_once_with(path.as_uri())
            self.assertTrue(path.is_file())
            self.assertIn("Open the HTML", output.getvalue())


if __name__ == "__main__":
    unittest.main()
