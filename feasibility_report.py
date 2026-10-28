"""Offline, filterable screening maps plus complete reproducible exports."""

import csv
import html
import json

from confinement import BETA_SOURCE, CONFINEMENT_SOURCE, DENSITY_SOURCE, HMODE_SOURCE

REASONS = {
    'radiation_exceeds_total_loss': 'Radiation exceeds total loss',
    'self_heating_surplus_no_fixed_point': 'Self-heating surplus',
    'heating_capacity': 'Heating capacity', 'greenwald_density': 'Density limit',
    'thermal_beta': 'Thermal pressure limit', 'hmode_threshold': 'H-mode threshold',
    'radiation_fraction': 'Radiation fraction',
}


def _table(headers, rows):
    header = ''.join(f'<th scope="col">{html.escape(str(h))}</th>' for h in headers)
    body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(v))}</td>' for v in row) + '</tr>' for row in rows)
    return f'<div class="table-scroll"><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>'


def _format(value, divisor=1):
    return '—' if value is None else f'{value / divisor:.5g}'


def _heatmap(points, h, fraction, identifier):
    temperatures = sorted({p['temperature_kev'] for p in points})
    densities = sorted({p['density_m3'] for p in points}, reverse=True)
    width = max(720, 130 + 95 * len(temperatures))
    cell_width = (width - 130) / len(temperatures)
    height = 110 + 70 * len(densities)
    svg = [f'<svg viewBox="0 0 {width} {height}" role="img" aria-labelledby="{identifier}-title">',
           f'<title id="{identifier}-title">Sampled net pulse energy: H={h:g}, T fraction={fraction:g}</title>']
    for i, density in enumerate(densities):
        svg.append(f'<text x="95" y="{60 + 70 * i}" text-anchor="end">{density / 1e20:.4g}</text>')
        for j, temperature in enumerate(temperatures):
            point = next(p for p in points if p['temperature_kev'] == temperature and p['density_m3'] == density)
            x, y = 110 + j * cell_width, 25 + i * 70
            energy = point['net_electric_energy_j']
            color = '#e2e8f0' if energy is None else '#bbf7d0' if energy > 0 else '#fecdd3'
            label = 'excluded' if energy is None else _format(energy, 1e6)
            details = ('; '.join(REASONS[r] for r in point['exclusion_reasons']) if energy is None
                       else f"Net pulse energy {energy / 1e6:.5g} MJ; net power {point['net_electric_w'] / 1e6:.5g} MW")
            title = html.escape(f'kT={temperature:g} keV, n={density:.5g} m^-3: {details}')
            svg.append(f'<rect x="{x:.2f}" y="{y}" width="{cell_width - 3:.2f}" height="66" '
                       f'fill="{color}" stroke="white"><title>{title}</title></rect>'
                       f'<text x="{x + cell_width / 2:.2f}" y="{y + 39}" text-anchor="middle">{label}</text>')
    for j, temperature in enumerate(temperatures):
        svg.append(f'<text x="{110 + (j + .5) * cell_width:.2f}" y="{45 + 70 * len(densities)}" '
                   f'text-anchor="middle">{temperature:g}</text>')
    svg.append(f'<text x="{width / 2:.2f}" y="{height - 12}" text-anchor="middle">Temperature kT (keV; sampled columns)</text>'
               f'<text transform="translate(22 {25 + 35 * len(densities)}) rotate(-90)" '
               'text-anchor="middle">Ion density (10²⁰ m⁻³)</text></svg>')
    return ''.join(svg)


def render_study_report(study):
    rows, summary, config = study['rows'], study['summary'], study['config']
    best = study['best_eligible_point']
    selected_h = best['h_factor'] if best else study['grid']['h_factors'][0]
    selected_f = best['tritium_fraction'] if best else study['grid']['tritium_fractions'][0]
    h_options = ''.join(f'<option value="{h:g}"{" selected" if h == selected_h else ""}>{h:g}</option>'
                        for h in study['grid']['h_factors'])
    f_options = ''.join(f'<option value="{f:g}"{" selected" if f == selected_f else ""}>{f:g}</option>'
                        for f in study['grid']['tritium_fractions'])
    panels = []
    for i, h in enumerate(study['grid']['h_factors']):
        for j, fraction in enumerate(study['grid']['tritium_fractions']):
            points = [p for p in rows if p['h_factor'] == h and p['tritium_fraction'] == fraction]
            hidden = '' if h == selected_h and fraction == selected_f else ' hidden'
            panels.append(f'<div class="heatmap-panel" data-h="{h:g}" data-f="{fraction:g}"{hidden}>'
                          + _heatmap(points, h, fraction, f'map-{i}-{j}') + '</div>')
    counts = _table(('Screening result', 'Count'), [
        ('Sampled target states', summary['sampled_points']),
        ('Pass all selected baseline screens', summary['eligible_points']),
        ('Positive operating electrical power', summary['positive_power_points']),
        ('Positive flat-top energy after startup', summary['positive_pulse_energy_points']),
        ('Positive and admissible at every tested corner', summary['all_corner_positive_points']),
    ])
    limits = _table(('Exclusion reason', 'Affected points'),
                    [(REASONS[k], v) for k, v in summary['exclusion_counts'].items()])
    sorted_rows = sorted((r for r in rows if r['passes_screens']),
                         key=lambda row: row['net_electric_energy_j'], reverse=True)
    leaders = _table(('kT (keV)', 'n (10²⁰ m⁻³)', 'H', 'T fraction', 'Aux heat (MW)',
                      'Transport time (s)', 'Net power (MW)', 'Net energy (MJ)', 'Startup payback (s)'),
                     [(_format(p['temperature_kev']), _format(p['density_m3'], 1e20), _format(p['h_factor']),
                       _format(p['tritium_fraction']), _format(p['required_heating_w'], 1e6),
                       _format(p['transport_confinement_s']), _format(p['net_electric_w'], 1e6),
                       _format(p['net_electric_energy_j'], 1e6), _format(p['startup_payback_s']))
                      for p in sorted_rows[:10]]) if sorted_rows else '<p>No sampled point passes every selected screen.</p>'
    sensitivity = _table(('Changed assumption', 'Value', 'Net energy (MJ)', 'Change (MJ)', 'Status'), [
        (r['parameter'].replace('_', ' '), r['value'], _format(r['net_energy_j'], 1e6),
         _format(r['delta_energy_j'], 1e6), '; '.join(REASONS[s] for s in r['exclusion_reasons']) or 'passes screens')
        for r in study['best_point_sensitivity']]) if best else '<p>No eligible baseline for sensitivity analysis.</p>'
    corner_text = ''
    if best:
        c = best['sensitivity']
        corner_text = (f'<p>For the highest-ranked baseline, {c["eligible_scenarios"]} of '
                       f'{c["scenario_count"]} tested corner evaluations pass the screens; '
                       f'{c["positive_energy_scenarios"]} have positive pulse energy. '
                       f'The valid-corner energy range is {_format(c["worst_valid_energy_j"], 1e6)} to '
                       f'{_format(c["best_valid_energy_j"], 1e6)} MJ. A dash means no corner remained admissible.</p>')
    roots = _table(('kT (keV)', 'n (10²⁰ m⁻³)', 'H', 'T fraction', 'Energy residual (J)'), [
        (_format(r['temperature_kev']), _format(r['density_m3'], 1e20), _format(r['h_factor']),
         _format(r['tritium_fraction']), _format(r['net_electric_energy_j'])) for r in study['break_even_roots']
    ]) if study['break_even_roots'] else '<p>No eligible adjacent samples bracket a pulse-energy break-even crossing.</p>'
    all_cases = _table(('kT (keV)', 'n (10²⁰ m⁻³)', 'H', 'T fraction', 'Beta N', 'Greenwald fraction',
                       'Net energy (MJ)', 'Exclusions'), [
        (_format(p['temperature_kev']), _format(p['density_m3'], 1e20), _format(p['h_factor']),
         _format(p['tritium_fraction']), _format(p['beta_n']), _format(p['greenwald_fraction']),
         _format(p['net_electric_energy_j'], 1e6), '; '.join(REASONS[r] for r in p['exclusion_reasons']) or 'passes screens')
        for p in rows])
    labels = {
        'major_radius_m': 'Major radius (m)', 'minor_radius_m': 'Minor radius (m)', 'elongation': 'Effective elongation',
        'field_t': 'Toroidal field (T)', 'current_ma': 'Plasma current (MA)',
        'radiation_convention': 'Power convention for confinement scaling', 'alpha_deposition': 'Alpha deposition fraction',
        'heater_efficiency': 'Electricity-to-heating efficiency', 'thermal_efficiency': 'Heat-to-electricity efficiency',
        'heat_recovery': 'Outgoing heat recovery fraction', 'plant_load_mw': 'Other electrical plant load (MW)',
        'extra_radiation_mw': 'Extra core radiation allowance (MW)', 'heating_capacity_mw': 'Deposited heating capacity (MW)',
        'greenwald_fraction_limit': 'Maximum Greenwald density fraction', 'beta_n_limit': 'Maximum thermal normalized beta',
        'hmode_margin': 'Minimum separatrix power / L-H threshold', 'radiation_fraction_limit': 'Maximum core radiation fraction',
        'pulse_duration_s': 'Assumed flat-top duration (s)', 'startup_multiplier': 'Thermal startup cost multiplier',
        'sensitivity_fraction': 'Selected relative sensitivity excursion',
    }
    input_rows = [(labels[k], v) for k, v in {**config['device'], **{k: v for k, v in config.items() if k != 'device'}}.items()]
    inputs = _table(('Reference input', 'Value'), input_rows)
    source_links = ', '.join(f'<a href="{url}">{title}</a>' for title, url in (
        ('UKAEA confinement and radiation conventions', CONFINEMENT_SOURCE), ('UKAEA L-H threshold', HMODE_SOURCE),
        ('UKAEA density limit', DENSITY_SOURCE), ('UKAEA pressure screening', BETA_SOURCE)))
    return f'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>IGNIS · Confinement screening study</title>
<style>body{{font:17px/1.6 system-ui,sans-serif;color:#172033;background:#f8fafc;max-width:1150px;margin:32px auto;padding:0 20px}}
h1{{line-height:1.15}}a{{color:#1d4ed8}}.eyebrow{{color:#b45309;font-weight:700;letter-spacing:.15em}}
.panel{{background:white;border:1px solid #e2e8f0;border-radius:14px;padding:20px;margin:24px 0}}
svg{{width:100%;height:auto}}svg text{{font-size:14px;fill:currentColor}}.table-scroll{{overflow-x:auto}}
table{{width:100%;border-collapse:collapse}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #e2e8f0}}
select{{font:inherit;padding:8px;margin:8px 16px 8px 4px}}[hidden]{{display:none}}</style></head><body>
<p class="eyebrow">IGNIS / PROTOTYPE 05</p><h1>Confinement and break-even screening</h1>
<p>This study compares uniform D–T target states against an illustrative large
<strong>tokamak reference</strong>. It does not select or validate the confinement
mechanism of the original D–D fuel-cycle concept.</p>
<section class="panel"><h2>What the scan found</h2>{counts}
<p>Positive numerical output is conditional on the geometry, operating regime,
maintained fuel mixture and conversion assumptions. Passing these screens is not
proof of a stable equilibrium, an achievable device or engineering feasibility.</p></section>
<section class="panel"><h2>Explore the sampled operating window</h2>
<label>Confinement multiplier H <select id="h-select">{h_options}</select></label>
<label>Tritium number fraction <select id="f-select">{f_options}</select></label>
{''.join(panels)}<p>Cells show net flat-top electrical energy after startup, in MJ.
Green: positive; pink: nonpositive; gray: excluded. Hover a cell for the budget or
exclusion reason. Columns are sampled temperatures, without interpolation.</p>
<noscript><p>All sampled cases are also available in the full table below.</p></noscript></section>
<section class="panel"><h2>Highest-ranked eligible samples</h2>{leaders}
<p>This ranks the finite grid by net pulse energy. It is not a global optimization.
Startup payback is meaningful only for positive operating electrical power.</p></section>
<section class="panel"><h2>Which assumptions change the result?</h2>{corner_text}{sensitivity}
<p>The table changes one selected input at a time at the highest-ranked baseline.
The corner check varies H, thermal efficiency, heat recovery, heater efficiency,
plant load and startup cost, across both radiation conventions. Relative excursions
are {config['sensitivity_fraction']:.0%}; fractions are bounded by 1 and startup cost
by the minimum thermal allowance. Corner ranges include eligible corners only;
excluded corners fail the all-corner check. These are deterministic stress tests,
not probabilities, statistical confidence limits or guarantees between corners.</p></section>
<section class="panel"><h2>Bracketed pulse-energy break-even</h2>{roots}
<p>Roots are refined only between adjacent eligible temperature samples at the same
density, H and mixture. Excluded gaps are not bridged; unsampled crossings can be missed.</p></section>
<section class="panel"><h2>Why cases were excluded</h2>{limits}
<p>Counts overlap when a point fails multiple screens. Negative required auxiliary
heating is a self-heating surplus without a fixed point at the requested target
state; it is excluded rather than silently adding an exhaust controller.</p>
<details><summary>Show all sampled cases</summary>{all_cases}</details></section>
<section class="panel"><h2>Reference geometry and assumptions</h2>{inputs}
<p>IPB98(y,2) gives tau = A × P⁻⁰·⁶⁹. The coupled loss power is solved from
U = tau × P, so confinement is not chosen independently of heating. The default
radiation-corrected convention uses transported power in the scaling and adds
radiation separately. The alternate convention uses total loss power and subtracts
radiation to find transported power. Radiation is counted once in either case.</p>
<p>Martin's mass-adjusted L-H threshold, a Greenwald density screen, a thermal
normalized-beta ceiling, heating capacity and a selected radiation-fraction ceiling
restrict the comparison. The thresholds are necessary screens and modelling choices;
they do not establish the validity of an empirical extrapolation. The uncorrected
L-H reference is restricted to aspect ratio greater than 2.7.</p>
<p>Net pulse energy assumes a constant flat-top state and maintained fuel inventory:
net operating power × duration minus the selected initial thermal startup allowance.
It is not a time-integrated burn prediction. Fuel supply, tritium production, current
drive, ash, profiles, fast particles, ramp losses, equilibrium, MHD stability,
blanket/material limits and detailed plant engineering remain outside this study.
The earlier burn mode continues to model the closed-volume D–D/D–T network.</p>
<p>Sources: {source_links}. See FEASIBILITY_MODEL.md for equations, conventions,
applicability and reproducibility.</p></section>
<script>
const hSelect=document.getElementById('h-select');
const fSelect=document.getElementById('f-select');
function updateMap(){{for(const panel of document.querySelectorAll('.heatmap-panel')){{
 panel.hidden=panel.dataset.h!==hSelect.value||panel.dataset.f!==fSelect.value;
}}}}
hSelect.addEventListener('change',updateMap);fSelect.addEventListener('change',updateMap);
</script></body></html>'''


def write_study_csv(path, study):
    input_fields = {**study['config']['device'], **{k: v for k, v in study['config'].items() if k != 'device'}}
    input_fields = {'input_' + k: v for k, v in input_fields.items()}
    flat_rows = []
    for row in study['rows']:
        flat = {k: v for k, v in row.items() if k not in ('sensitivity', 'exclusion_reasons')}
        flat['exclusion_reasons'] = ';'.join(row['exclusion_reasons'])
        sensitivity = row['sensitivity'] or {}
        for key in ('scenario_count', 'eligible_scenarios', 'positive_energy_scenarios',
                    'worst_valid_energy_j', 'best_valid_energy_j', 'all_corners_positive'):
            flat['sensitivity_' + key] = sensitivity.get(key)
        flat_rows.append({**input_fields, **flat})
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat_rows[0]))
        writer.writeheader()
        writer.writerows(flat_rows)


def write_study_json(path, study):
    payload = {**study, 'sources': [CONFINEMENT_SOURCE, HMODE_SOURCE, DENSITY_SOURCE, BETA_SOURCE],
               'scope': 'Conditional uniform D-T tokamak reference screening; not a device feasibility verdict.'}
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')
