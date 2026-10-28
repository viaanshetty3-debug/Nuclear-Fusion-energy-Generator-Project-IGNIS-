"""Standalone SVG report and finite CSV/JSON exports for burn trajectories."""

from dataclasses import asdict
import csv
import html
import json
from pathlib import Path

from fusion import NRL_URL, dd_reactivities

COLORS = ('#2563eb', '#b45309', '#15803d', '#9333ea', '#be123c')


def dd_reference_comparison():
    path = Path(__file__).resolve().parent / 'data' / 'nrl_dd_reference.csv'
    with path.open(encoding='utf-8', newline='') as stream:
        reference = list(csv.DictReader(stream))
    result = []
    for row in reference:
        temperature = float(row['temperature_kev'])
        nrl = float(row['reactivity_cm3_s']) * 1e-6
        fitted = sum(dd_reactivities(temperature))
        result.append({'temperature_kev': temperature, 'nrl_dd_reactivity_m3_s': nrl,
                       'bosch_hale_dd_reactivity_m3_s': fitted, 'relative_difference': fitted / nrl - 1})
    return result


def _table(headers, rows):
    header = ''.join(f'<th scope="col">{html.escape(str(h))}</th>' for h in headers)
    body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(v))}</td>' for v in row) + '</tr>' for row in rows)
    return f'<div class="table-scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'


def _chart(rows, series, title, unit, identifier):
    xmax = max(rows[-1]['time_s'], 1e-12)
    curves = [(label, [function(row) for row in rows]) for label, function in series]
    values = [v for _, points in curves for v in points]
    low, high = min(0.0, min(values)), max(0.0, max(values))
    pad = (high - low) * 0.08 if high != low else 1.0
    low, high = low - pad, high + pad
    def x(time):
        return 110 + 780 * time / xmax
    def y(value):
        return 385 - 335 * (value - low) / (high - low)
    svg = [f'<svg viewBox="0 0 950 475" role="img" aria-labelledby="{identifier}-title">',
           f'<title id="{identifier}-title">{html.escape(title)}</title>']
    for i in range(6):
        value = low + (high - low) * i / 5
        time = xmax * i / 5
        svg.append(f'<path d="M110 {y(value):.2f} H890" stroke="#e2e8f0"/>'
                   f'<text x="100" y="{y(value) + 5:.2f}" text-anchor="end">{value:.3g}</text>'
                   f'<text x="{x(time):.2f}" y="415" text-anchor="middle">{time:.3g}</text>')
    svg.append(f'<path d="M110 50 V385 H890 M110 {y(0):.2f} H890" stroke="#64748b" fill="none"/>')
    legend = []
    for i, (label, points) in enumerate(curves):
        color = COLORS[i % len(COLORS)]
        coords = ' '.join(f"{x(row['time_s']):.2f},{y(value):.2f}" for row, value in zip(rows, points))
        svg.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="2">'
                   f'<title>{html.escape(label)}</title></polyline>')
        legend.append(f'<span style="color:{color}">● {html.escape(label)}</span>')
    svg.append('<text x="500" y="458" text-anchor="middle">Time (s)</text>'
               f'<text transform="translate(25 220) rotate(-90)" text-anchor="middle">{html.escape(unit)}</text></svg>')
    return f'<section class="panel"><h2>{html.escape(title)}</h2>' + ''.join(svg) + '<p>' + ' &nbsp; '.join(legend) + '</p></section>'


def render_burn_report(result, convergence):
    config, rows = result.config, result.rows
    first, last = rows[0], rows[-1]
    labels = {
        'temperature_kev': 'Initial ion/electron kT (keV)', 'density_m3': 'Initial total fuel ion density (m⁻³)',
        'tritium_fraction': 'Initial tritium number fraction', 'volume_m3': 'Fixed volume (m³)',
        'confinement_s': 'Transport confinement time (s)', 'charged_deposition': 'Charged birth-energy deposition fraction',
        'heating_mw': 'Constant deposited auxiliary heating (MW)', 'heater_efficiency': 'Heater efficiency',
        'thermal_efficiency': 'Heat-to-electricity efficiency', 'heat_recovery': 'Heat recovery fraction',
        'plant_load_mw': 'Other plant electrical load (MW)', 'other_loss_mw': 'Other plasma loss (MW)',
        'duration_s': 'Requested duration (s)', 'max_step_s': 'Maximum integration step (s)',
        'rtol': 'Relative integration tolerance', 'samples': 'Requested output samples',
    }
    inputs = _table(('Input', 'Value'), [(labels[k], f'{v:.6g}') for k, v in asdict(config).items()])
    summary = _table(('Result', 'Value'), [
        ('Termination', result.status.replace('_', ' ')), ('Simulated duration (s)', f"{last['time_s']:.8g}"),
        ('Final temperature (keV)', f"{last['temperature_kev']:.6g}"),
        ('D–D tritium births', f"{last['dd_tritium_reactions']:.6g}"),
        ('Subsequent D–T reactions', f"{last['dt_reactions']:.6g}"),
        ('Tritium remaining (m⁻³)', f"{last['tritium_density_m3']:.6g}"),
        ('Fusion energy (MJ)', f"{last['fusion_energy_j'] / 1e6:.6g}"),
        ('Initial thermal startup electricity (MJ)', f"{last['startup_electric_energy_j'] / 1e6:.6g}"),
        ('Operating net electrical energy (MJ)', f"{last['operating_net_electric_energy_j'] / 1e6:.6g}"),
        ('Net electrical energy including startup (MJ)', f"{last['net_electric_energy_j'] / 1e6:.6g}"),
    ])
    checks = _table(('Check', 'Result'), [
        ('Maximum relative plasma energy residual', f'{result.max_energy_relative_error:.3g}'),
        ('Maximum relative charge residual', f'{result.max_charge_relative_error:.3g}'),
        ('Maximum relative nucleon residual (including escaped neutrons)', f'{result.max_nucleon_relative_error:.3g}'),
        ('Accepted / rejected steps', f'{result.accepted_steps} / {result.rejected_steps}'),
        ('Tighter-run convergence', 'PASS' if convergence['passed'] else 'FAIL'),
        ('Matching termination reason', str(convergence['same_termination'])),
        ('Refined maximum step (s)', f"{convergence['refined_max_step_s']:.5g}"),
        ('Refined relative tolerance', f"{convergence['refined_rtol']:.5g}"),
    ])
    comparisons = _table(('Final quantity', 'Absolute difference', 'Allowed difference', 'Result'), [
        (key, f"{values['absolute_difference']:.5g}", f"{values['tolerance']:.5g}",
         'PASS' if values['passed'] else 'FAIL') for key, values in convergence['comparisons'].items()])
    references = _table(('kT (keV)', 'NRL DD total (m³/s)', 'Bosch–Hale DD total (m³/s)', 'Difference'), [
        (f"{r['temperature_kev']:g}", f"{r['nrl_dd_reactivity_m3_s']:.5g}",
         f"{r['bosch_hale_dd_reactivity_m3_s']:.5g}", f"{r['relative_difference']:+.2%}")
        for r in dd_reference_comparison()])
    n = config.density_m3
    plots = _chart(rows, [('Temperature', lambda r: r['temperature_kev'])],
                   'Plasma temperature', 'kT (keV)', 'burn-temperature')
    plots += _chart(rows, [('Deuterium', lambda r: r['deuterium_density_m3'] / n),
                          ('Tritium', lambda r: r['tritium_density_m3'] / n)],
                    'Fuel inventory', 'Fraction of initial ion count', 'burn-fuel')
    plots += _chart(rows, [(name, lambda r, key=key: r[key] / n) for name, key in (
        ('Tritium remaining', 'tritium_density_m3'), ('Protons', 'proton_density_m3'),
        ('Helium-3', 'helium3_density_m3'), ('Helium-4 ash', 'helium4_density_m3'))],
        'Retained products and tritium', 'Count per initial ion', 'burn-products')
    plots += _chart(rows, [
        ('Total fusion', lambda r: r['fusion_power_w'] / 1e6),
        ('D–T', lambda r: r['dt_reactions_s'] * 17600 * 1.602176634e-16 / 1e6),
        ('D–D → T+p', lambda r: r['dd_tritium_reactions_s'] * 4030 * 1.602176634e-16 / 1e6),
        ('D–D → He3+n', lambda r: r['dd_helium3_reactions_s'] * 3270 * 1.602176634e-16 / 1e6)],
        'Fusion power over time', 'Power (MW)', 'burn-power')
    plots += _chart(rows, [('Fusion energy', lambda r: r['fusion_energy_j'] / 1e6),
                          ('Net electricity including startup', lambda r: r['net_electric_energy_j'] / 1e6)],
                    'Cumulative energy', 'Energy (MJ)', 'burn-energy')
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>IGNIS · Time-dependent burn</title>
<style>body{{font:17px/1.6 system-ui,sans-serif;color:#172033;background:#f8fafc;max-width:1100px;
margin:32px auto;padding:0 20px}}h1{{line-height:1.15}}a{{color:#1d4ed8}}
.eyebrow{{color:#b45309;font-weight:700;letter-spacing:.15em}}.panel{{background:white;border:1px solid #e2e8f0;
border-radius:14px;padding:20px;margin:24px 0}}svg{{width:100%;height:auto}}svg text{{font-size:15px;fill:currentColor}}
.table-scroll{{overflow-x:auto}}table{{width:100%;border-collapse:collapse}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #e2e8f0}}</style>
</head><body><p class="eyebrow">IGNIS / PROTOTYPE 04</p><h1>A fuel burn over time</h1>
<p>Closed, fixed-volume plasma with evolving D–D and D–T reactions, product buildup,
and one shared ion/electron temperature. Constant auxiliary heating is prescribed;
temperature evolves from the energy balance.</p>
<section class="panel"><h2>Burn outcome</h2>{summary}
<p>A temperature-bound termination ends the trajectory at the last valid state,
within 10⁻⁸ keV of the 0.2–30 keV model boundary. It does not extrapolate the reaction or loss model.
A completed run means the requested duration was simulated, not that a device is feasible.</p></section>
{plots}<section class="panel"><h2>Inputs and assumptions</h2>{inputs}
<p>Both D–D branches use independent Bosch–Hale reactivities with the identical-reactant
factor 1/2. D–T uses distinct reactants. Charged products are retained and assumed to
join the thermal population immediately; their deposited birth-energy fraction is
an input. Neutrons escape and are included in cumulative nucleon accounting.</p>
<p>All ions are fully ionized; electron density follows charge neutrality, and thermal
capacity and bremsstrahlung change as helium accumulates. Transport removes energy
without prescribed particle loss. There is no refuelling or ash extraction.</p>
<p>Electrical output follows outgoing neutron/charged birth energy plus plasma heat
losses, with the selected recovery and conversion efficiencies. Stored plasma heat
is not automatically credited at the end. The startup cost is the initial thermal
inventory divided by heater efficiency; startup losses, ionization and equipment
preparation would add costs. Net energy starts negative by that startup allowance.</p>
<p>The reaction network omits D–He3, T–T and other secondary channels. Fast-ion slowing,
nonthermal reactions, separate electron/ion temperatures, compression, spatial
transport, particle escape, confinement feasibility, stability and engineering
limits are outside the model. Fixed MeV birth-energy partitions neglect small
thermal-energy corrections. Results are conditional on these approximations.</p></section>
<section class="panel"><h2>Conservation and numerical convergence</h2>{checks}{comparisons}
<p>The comparison reruns the burn with a smaller maximum step and a tighter error
tolerance (down to the 10⁻¹¹ floor). Displayed limits combine relative and quantity-specific absolute tolerances;
numerical agreement does not establish physical accuracy.</p>
<h3>Independent D–D reference comparison</h3>{references}
<p>The summed branch reactivities differ from the older rounded NRL table by about
34% at 1 keV, 16% at 2 keV, and under 10% at the listed 5–100 keV points. These
cross-checks do not assign a physical uncertainty to the fit or validate the burn model.</p>
<p>Sources: <a href="{NRL_URL}#page=45">NRL reaction products</a>,
<a href="{NRL_URL}#page=59">NRL radiation</a>,
<a href="https://doi.org/10.1088/0029-5515/32/4/i07">Bosch &amp; Hale (1992)</a>,
<a href="https://github.com/ukaea/PROCESS/blob/main/process/models/physics/fusion_reactions.py">UKAEA branch coefficients</a>.</p>
</section></body></html>'''


def write_burn_csv(path, result):
    config = {'initial_' + k: v for k, v in asdict(result.config).items()}
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[*config, *result.rows[0]])
        writer.writeheader()
        writer.writerows({**config, **row} for row in result.rows)


def write_burn_json(path, result, convergence):
    payload = {'schema_version': 1, 'model': 'closed_volume_dd_dt_burn', **result.to_dict(),
               'convergence': convergence, 'dd_reference_comparison': dd_reference_comparison(),
               'assumptions': ['Closed fixed volume; no refuelling or particle transport.',
                               'Instant thermalization of retained charged products; single Te=Ti.',
                               'Escaped neutrons included in cumulative nucleon inventory.',
                               'Fixed auxiliary heating; no temperature thermostat.',
                               'Net electricity charges minimum initial thermal startup energy.',
                               'Network only DT, DD->T+p, DD->He3+n; other reactions omitted.'],
               'sources': [NRL_URL, 'https://doi.org/10.1088/0029-5515/32/4/i07']}
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')
