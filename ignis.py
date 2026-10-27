"""IGNIS thermal ion-energy and D-T power-balance explorer. Uses only the Python standard library."""

import argparse
from dataclasses import asdict
import csv
import math
from pathlib import Path
import webbrowser

from fusion import PlantConfig, evaluate_case
from fusion_report import render_fusion_sections, write_summary_csv, write_json

# Exact SI defining constants: joules per kelvin / joules per keV.
K = 1.380649e-23 / (1000 * 1.602176634e-19)
TEMPERATURES = (50_000_000, 100_000_000, 200_000_000)
COLORS = ("#2563eb", "#b45309", "#15803d", "#9333ea", "#be123c", "#0e7490")


def _finite(value, name, *, positive=False):
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        bound = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {bound}")
    return value


def thermal_energy(temperature):
    """Return kT in keV for a finite, positive temperature in kelvin."""
    _finite(temperature, "temperature", positive=True)
    scale = K * temperature
    _finite(scale, "thermal energy", positive=True)
    return scale


def probability(energy, temperature):
    """Maxwellian single-ion kinetic-energy density in 1/keV."""
    _finite(energy, "energy")
    scale = thermal_energy(temperature)
    x = energy / scale
    if energy == 0 or not math.isfinite(x):
        return 0.0
    # Log form avoids inf * 0 for very large energy/temperature ratios.
    return math.exp(math.log(2 / math.sqrt(math.pi)) + 0.5 * (math.log(energy) - math.log(scale))
                    - x - math.log(scale))


def fraction_above(energy, temperature):
    """Analytic probability P(E >= energy); this is not a fusion probability."""
    _finite(energy, "energy")
    x = energy / thermal_energy(temperature)
    if not math.isfinite(x):
        return 0.0
    root = math.sqrt(x)
    return min(1.0, math.erfc(root) + 2 / math.sqrt(math.pi) * root * math.exp(-x))


def sample_curves(temperatures, samples=1001):
    """Sample each curve on its own quadratic grid through 20 kT.

    Independent grids resolve cold curves even when temperatures differ widely.
    Extra resolution near E=0 captures the square-root behavior of the density.
    """
    if not temperatures:
        raise ValueError("at least one temperature is required")
    if not isinstance(samples, int) or not 51 <= samples <= 100_001:
        raise ValueError("samples must be an integer between 51 and 100001")
    curves = []
    for temperature in temperatures:
        scale = thermal_energy(temperature)
        xmax = 20 * scale
        _finite(xmax, "plot energy range", positive=True)
        energies = [xmax * (i / (samples - 1)) ** 2 for i in range(samples)]
        curves.append((temperature, energies,
                       [probability(e, temperature) for e in energies]))
    return curves


def render_html(curves, threshold):
    """Build a standalone report with inline SVG and a summary table."""
    _finite(threshold, "threshold")
    xmax = max(energies[-1] for _, energies, _ in curves)
    ymax = 1.1 * max(max(values) for _, _, values in curves)
    svg = [
        '<svg viewBox="0 0 960 520" role="img" aria-labelledby="plot-title plot-desc">',
        '<title id="plot-title">Maxwellian ion kinetic-energy distributions</title>',
        '<desc id="plot-desc">Energy in keV versus probability density in inverse keV. '
        'Temperatures and exact distribution summaries appear in the table below.</desc>',
        '<path d="M100 35 V435 H920" fill="none" stroke="currentColor"/>',
    ]
    for i in range(6):
        x, y = 100 + 164 * i, 435 - 80 * i
        svg.extend([
            f'<path d="M100 {y} H920" stroke="#e2e8f0"/>',
            f'<text x="{x}" y="460" text-anchor="middle">{xmax * i / 5:.3g}</text>',
            f'<text x="90" y="{y + 5}" text-anchor="end">{ymax * i / 5:.3g}</text>',
        ])
    if threshold <= xmax:
        x = 100 + 820 * (threshold / xmax)
        svg.append(f'<path d="M{x:.2f} 35 V435" stroke="#475569" '
                   'stroke-dasharray="5 5"><title>Selected energy threshold</title></path>')
    rows = []
    for i, (temperature, energies, values) in enumerate(curves):
        color = COLORS[i % len(COLORS)]
        points = " ".join(f"{100 + 820 * (e / xmax):.2f},{435 - 400 * (p / ymax):.2f}"
                          for e, p in zip(energies, values))
        svg.append(f'<polyline points="{points}" fill="none" stroke="{color}" '
                   f'stroke-width="2"><title>{temperature:,.9g} K</title></polyline>')
        scale = thermal_energy(temperature)
        rows.append(
            f'<tr><th scope="row"><span style="color:{color}">●</span> '
            f'{temperature:,.9g}</th><td>{scale:.6g}</td><td>{1.5 * scale:.6g}</td>'
            f'<td>{fraction_above(threshold, temperature):.6%}</td></tr>'
        )
    svg.extend([
        '<text x="510" y="500" text-anchor="middle">Ion kinetic energy (keV)</text>',
        '<text transform="translate(22 235) rotate(-90)" text-anchor="middle">'
        'Probability density (1/keV)</text></svg>',
    ])
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>IGNIS · Thermal ion energy</title>
<style>
body{{font:17px/1.6 system-ui,sans-serif;color:#172033;background:#f8fafc;
max-width:1100px;margin:32px auto;padding:0 20px}}
h1{{line-height:1.15}} .eyebrow{{color:#b45309;letter-spacing:.15em;font-weight:700}}
.panel{{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:20px;margin:24px 0}}
svg{{width:100%;height:auto}}svg text{{font-size:15px;fill:currentColor}}
.table-scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%}}
th,td{{padding:12px;text-align:left;border-bottom:1px solid #e2e8f0}}
code{{font-size:.95em}}a{{color:#1d4ed8}}
</style></head><body>
<p class="eyebrow">IGNIS / PROTOTYPE 02</p>
<h1>Temperature and ion energy</h1>
<p>Explore the energy spread of classical, nonrelativistic Maxwellian ions.</p>
<section class="panel" aria-label="Energy distribution plot">{"".join(svg)}
<p>Each curve extends to 20 kT of its own temperature, covering more than
99.999998% of the distribution. The full curve has area 1.
The dashed line marks {threshold:.6g} keV when it lies within the plot.</p>
</section>
<section class="panel"><h2>Distribution summary</h2><div class="table-scroll">
<table><caption>Fraction of individual ions with kinetic energy at or above
{threshold:.6g} keV</caption><thead><tr><th scope="col">Temperature (K)</th>
<th scope="col">kT (keV)</th><th scope="col">Mean energy (keV)</th>
<th scope="col">Fraction ≥ {threshold:.6g} keV</th></tr></thead>
<tbody>{"".join(rows)}</tbody></table></div></section>
<section class="panel"><h2>What this model describes</h2>
<p>The mean ion energy is <code>3/2 kT</code>; doubling temperature doubles
 the mean energy. The density peaks at <code>E = kT/2</code>.</p>
<p>The threshold fraction is a thermal population statistic. Fusion rates require
relative collision energies and reaction cross sections. The thermal panel describes individual ions; the separate D–T module uses
a relative-velocity averaged reactivity for reaction and power calculations.</p>
<p>For very different temperatures, the shared linear axes can visually compress
one curve. The summary and CSV retain each distribution's values.</p>
<p>References: <a href="https://www.nist.gov/si-redefinition/meet-constants">NIST SI constants</a>;
<a href="https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf">2023 NRL Plasma Formulary</a>.</p>
</section></body></html>
'''


def write_csv(path, curves):
    """Export long-form samples with explicit units, one row per energy point."""
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("temperature_K", "energy_keV", "probability_density_per_keV"))
        for temperature, energies, values in curves:
            writer.writerows((temperature, e, p) for e, p in zip(energies, values))


def _positive(text):
    try:
        return _finite(float(text), "value", positive=True)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _nonnegative(text):
    try:
        return _finite(float(text), "value")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--temperatures", type=_positive, nargs="+", default=TEMPERATURES,
                        metavar="K", help="ion temperatures in kelvin (default: 5e7 1e8 2e8)")
    parser.add_argument("--threshold-kev", type=_nonnegative, default=10.0,
                        help="energy for the thermal tail statistic (default: 10)")
    parser.add_argument("--samples", type=int, default=1001,
                        help="points per thermal curve, 51–100001 (default: 1001)")
    parser.add_argument("--output", type=Path, default=Path("ignis_results.html"),
                        help="HTML output path, relative to the current directory")
    parser.add_argument("--csv", type=Path, metavar="PATH", help="export thermal curve samples")
    parser.add_argument("--no-browser", action="store_true", help="generate without opening a browser")
    parser.add_argument("--thermal-only", action="store_true", help="run the original ion-energy model only")
    group = parser.add_argument_group("D-T steady-state scenario (default mode; kT 0.2–30 keV)")
    group.add_argument("--densities", type=_positive, nargs="+", metavar="M^-3",
                       help="total ion densities (default: 5e19 1e20 2e20); every T × density is evaluated")
    descriptions = {
        "volume_m3": "uniform plasma volume (m^3)",
        "tritium_fraction": "tritium number fraction of D+T ions, 0–1",
        "confinement_s": "transport energy confinement time (s), excluding radiation",
        "alpha_deposition": "fraction of alpha power heating the plasma, 0–1",
        "heater_efficiency": "wall electricity to deposited heating efficiency, >0 to 1",
        "thermal_efficiency": "recovered heat to electricity efficiency, 0–1",
        "heat_recovery": "fraction of outgoing heat recovered for generation, 0–1",
        "plant_load_mw": "other plant electrical demand (MW)",
        "other_loss_mw": "additional total plasma loss allowance (MW)",
    }
    defaults = asdict(PlantConfig())
    for name, default in defaults.items():
        group.add_argument("--" + name.replace("_", "-"), type=_nonnegative,
                           help=f"{descriptions[name]} (default: {default:g})")
    group.add_argument("--fusion-csv", type=Path, metavar="PATH", help="export all case inputs and power budgets")
    group.add_argument("--json", type=Path, metavar="PATH", help="export scenario, cases and reference comparison")
    args = parser.parse_args(argv)
    cases, config = None, None
    try:
        paths = [p.resolve() for p in (args.output, args.csv, args.fusion_csv, args.json) if p is not None]
        if len(paths) != len(set(paths)):
            raise ValueError("all output paths must be different")
        config_inputs = {name: getattr(args, name) for name in defaults if getattr(args, name) is not None}
        if args.thermal_only and (config_inputs or args.densities or args.fusion_csv or args.json):
            raise ValueError("D-T scenario options cannot be used with --thermal-only")
        curves = sample_curves(args.temperatures, args.samples)
        report = render_html(curves, args.threshold_kev)
        if not args.thermal_only:
            config = PlantConfig(**config_inputs)
            densities = args.densities if args.densities is not None else (5e19, 1e20, 2e20)
            cases = [evaluate_case(thermal_energy(t), density, config)
                     for t in args.temperatures for density in densities]
            report = report.replace("IGNIS / PROTOTYPE 02", "IGNIS / PROTOTYPE 03")
            report = report.replace("<title>IGNIS · Thermal ion energy</title>",
                                    "<title>IGNIS · D–T fusion and energy balance</title>")
            report = report.replace("<h1>Temperature and ion energy</h1>",
                                    '<h1>D–T fusion and energy balance</h1>'
                                    '<p><a href="#fusion-model">Go to fusion power and energy budgets</a></p>')
            report = report.replace("</body>", render_fusion_sections(cases, config) + "</body>")
    except (ValueError, OverflowError) as exc:
        parser.error(str(exc))
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
        if args.csv:
            args.csv.parent.mkdir(parents=True, exist_ok=True)
            write_csv(args.csv, curves)
        if args.fusion_csv:
            args.fusion_csv.parent.mkdir(parents=True, exist_ok=True)
            write_summary_csv(args.fusion_csv, cases, config)
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            write_json(args.json, cases, config)
    except OSError as exc:
        parser.exit(1, f"Cannot write results: {exc}\n")
    for temperature in args.temperatures:
        print(f"{temperature:,.9g} K: mean energy = {1.5 * thermal_energy(temperature):.6g} keV; "
              f"fraction >= {args.threshold_kev:g} keV = "
              f"{fraction_above(args.threshold_kev, temperature):.6%}")
    if cases:
        print("D-T target states; conditional steady-state estimates (power in MW):")
        for case in cases:
            print(f"  kT={case.temperature_kev:.5g} keV, n={case.density_m3:.4g} m^-3: "
                  f"rate={case.reactions_s:.5g}/s, fusion={case.fusion_power_w / 1e6:.5g}, "
                  f"heating={case.required_heating_w / 1e6:.5g}, "
                  f"extra exhaust={case.extra_exhaust_w / 1e6:.5g}, "
                  f"net electric={case.net_electric_w / 1e6:.5g}")
    for label, path in (("HTML", args.output), ("Thermal CSV", args.csv),
                        ("Fusion CSV", args.fusion_csv), ("JSON", args.json)):
        if path:
            print(f"{label}: {path.resolve()}")
    if not args.no_browser:
        try:
            opened = webbrowser.open(args.output.resolve().as_uri())
        except webbrowser.Error:
            opened = False
        if not opened:
            print("Open the HTML path above in a browser to view the report.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
