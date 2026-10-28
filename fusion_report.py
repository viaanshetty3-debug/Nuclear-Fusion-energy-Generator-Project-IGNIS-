"""Offline report sections and reproducible exports for the D–T model."""

from dataclasses import asdict
import csv
import html
import json
from pathlib import Path

from fusion import NRL_URL, dt_reactivity, evaluate_case

COLORS = ("#2563eb", "#b45309", "#15803d", "#9333ea", "#be123c", "#0e7490")


def reference_comparison():
    """Independent rounded NRL benchmarks; differences are reported, not fitted away."""
    path = Path(__file__).resolve().parent / "data" / "nrl_dt_reference.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    results = []
    for row in rows:
        t = float(row["temperature_kev"])
        reference = float(row["reactivity_cm3_s"]) * 1e-6
        fitted = dt_reactivity(t)
        results.append({"temperature_kev": t, "nrl_reactivity_m3_s": reference,
                        "bosch_hale_reactivity_m3_s": fitted,
                        "relative_difference": fitted / reference - 1})
    return results


def write_summary_csv(path, cases, config):
    inputs = asdict(config)
    fields = list(inputs) + list(cases[0].to_dict())
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({**inputs, **case.to_dict()} for case in cases)


def write_json(path, cases, config):
    payload = {
        "schema_version": 1,
        "model": "uniform_dt_steady_state",
        "assumptions": [
            "Pure fully ionized D-T, ne=ni, Te=Ti, spatially uniform, prescribed density.",
            "Bosch-Hale thermal reactivity, no use of a single-ion energy threshold.",
            "Transport confinement time excludes separately calculated radiation.",
            "Minimum auxiliary heating; excess alpha heating requires added exhaust.",
            "Heat cycle recovers a fraction of fusion plus deposited auxiliary power.",
            "No burn evolution, engineering feasibility, stability or startup dynamics.",
        ],
        "config": asdict(config),
        "cases": [case.to_dict() for case in cases],
        "nrl_reference_comparison": reference_comparison(),
        "q_null_meaning": "No auxiliary heating required in the prescribed balance; Q undefined.",
        "sources": [NRL_URL, "https://doi.org/10.1088/0029-5515/32/4/i07",
                    "https://ukaea.github.io/PROCESS/physics-models/fusion_reactions/plasma_bosch_hale/"],
    }
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _table(caption, headers, rows):
    head = "".join(f'<th scope="col">{html.escape(label)}</th>' for label in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row)
                   + "</tr>" for row in rows)
    return (f'<div class="table-scroll"><table><caption>{html.escape(caption)}</caption>'
            f'<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>')


def _power_chart(cases, config, field, title, identifier):
    temps = sorted({case.temperature_kev for case in cases})
    densities = list(dict.fromkeys(case.density_m3 for case in cases))
    low, high = temps[0], temps[-1]
    if low == high:
        low, high = max(0.2, low * 0.8), min(30.0, high * 1.2)
    grid = [low + (high - low) * i / 80 for i in range(81)]
    curves = [(density, [getattr(evaluate_case(t, density, config), field) / 1e6 for t in grid])
              for density in densities]
    values = [value for _, curve in curves for value in curve]
    ymin, ymax = min(0.0, min(values)), max(0.0, max(values))
    padding = 0.08 * (ymax - ymin) if ymax != ymin else 1.0
    ymin, ymax = ymin - padding, ymax + padding
    def x(t):
        return 110 + 790 * (t - low) / (high - low)
    def y(power):
        return 395 - 350 * (power - ymin) / (ymax - ymin)
    svg = [f'<svg viewBox="0 0 960 485" role="img" aria-labelledby="{identifier}-title">',
           f'<title id="{identifier}-title">{html.escape(title)}</title>']
    for i in range(6):
        power = ymin + (ymax - ymin) * i / 5
        temperature = low + (high - low) * i / 5
        svg.append(f'<path d="M110 {y(power):.2f} H900" stroke="#e2e8f0"/>'
                   f'<text x="100" y="{y(power) + 5:.2f}" text-anchor="end">{power:.3g}</text>'
                   f'<text x="{x(temperature):.2f}" y="425" text-anchor="middle">{temperature:.3g}</text>')
    svg.append(f'<path d="M110 45 V395 H900 M110 {y(0):.2f} H900" '
               'stroke="#64748b" fill="none"/>')
    legend = []
    for i, (density, curve) in enumerate(curves):
        color = COLORS[i % len(COLORS)]
        points = " ".join(f"{x(t):.2f},{y(power):.2f}" for t, power in zip(grid, curve))
        svg.append(f'<polyline fill="none" stroke="{color}" stroke-width="2" points="{points}">'
                   f'<title>Ion density {density:.4g} m⁻³</title></polyline>')
        legend.append(f'<span style="color:{color}">● {density:.4g} m⁻³</span>')
    svg.append('<text x="505" y="465" text-anchor="middle">Ion temperature kT (keV)</text>'
               '<text transform="translate(25 220) rotate(-90)" text-anchor="middle">'
               'Power (MW)</text></svg>')
    return f'<h3>{title}</h3>' + "".join(svg) + '<p>' + " &nbsp; ".join(legend) + '</p>'


def render_fusion_sections(cases, config):
    labels = {
        "volume_m3": "Plasma volume (m³)",
        "tritium_fraction": "Tritium number fraction",
        "confinement_s": "Transport confinement time (s)",
        "alpha_deposition": "Alpha deposition fraction",
        "heater_efficiency": "Heater efficiency",
        "thermal_efficiency": "Heat-to-electricity efficiency",
        "heat_recovery": "Heat recovery fraction",
        "plant_load_mw": "Other plant electrical load (MW)",
        "other_loss_mw": "Other plasma loss allowance (MW)",
    }
    input_rows = [(labels[name], f"{value:.6g}") for name, value in asdict(config).items()]
    inputs = _table("Scenario inputs (SI units; powers labelled MW)", ("Input", "Value"), input_rows)
    summary = _table(
        "Selected temperature × total ion-density cases; all powers in MW",
        ("kT (keV)", "Ion density (m⁻³)", "⟨σv⟩ (m³/s)", "Reactions/s", "Fusion (MW)",
         "Required heating (MW)", "Plasma Q", "Net electricity (MW)"),
        [(f"{c.temperature_kev:.5g}", f"{c.density_m3:.5g}", f"{c.reactivity_m3_s:.5g}",
          f"{c.reactions_s:.5g}", f"{c.fusion_power_w / 1e6:.5g}",
          f"{c.required_heating_w / 1e6:.5g}",
          f"{c.plasma_gain_q:.5g}" if c.plasma_gain_q is not None else "— (no auxiliary heat)",
          f"{c.net_electric_w / 1e6:.5g}") for c in cases],
    )
    breakdown = _table(
        "Power accounting for each selected case (MW); stored thermal energy in MJ",
        ("kT (keV)", "Ion density (m⁻³)", "Stored (MJ)", "Alpha", "Neutron", "Alpha deposited",
         "Transport loss", "Bremsstrahlung", "Other loss", "Extra exhaust required",
         "Recovered heat", "Gross electric", "Heater electric", "Plant load"),
        [(f"{c.temperature_kev:.5g}", f"{c.density_m3:.5g}",
          *[f"{value / 1e6:.5g}" for value in (
              c.stored_energy_j, c.alpha_power_w, c.neutron_power_w, c.alpha_heating_w,
              c.transport_loss_w, c.bremsstrahlung_loss_w, c.other_loss_w, c.extra_exhaust_w,
              c.recovered_heat_w, c.gross_electric_w, c.heater_electric_w, c.plant_load_w)])
         for c in cases],
    )
    references = _table(
        "Independent NRL 2023 D–T comparison (p45); cm³/s converted to m³/s",
        ("kT (keV)", "NRL (m³/s)", "Bosch–Hale (m³/s)", "Difference"),
        [(f"{r['temperature_kev']:g}", f"{r['nrl_reactivity_m3_s']:.5g}",
          f"{r['bosch_hale_reactivity_m3_s']:.5g}", f"{r['relative_difference']:+.2%}")
         for r in reference_comparison()],
    )
    charts = _power_chart(cases, config, "fusion_power_w", "Gross fusion power", "fusion-power")
    charts += _power_chart(cases, config, "net_electric_w", "Estimated net electrical power", "net-power")
    return f'''<section class="panel" id="fusion-model"><h2>D–T reactions and energy balance</h2>
<p>Uniform, fully ionized D–T with equal ion and electron temperatures. Density is
<strong>total fuel ion density</strong>; the tritium fraction splits it into D and T.
Volume and confinement are prescribed inputs. The model evaluates target states.</p>
{inputs}{summary}
<p>Positive net electricity is conditional on these inputs and the assumed heat
recovery and exhaust being achievable. It does not establish reactor feasibility.
Plasma Q uses fusion power divided by deposited auxiliary heating and is distinct
from plant net electricity. A dash means zero auxiliary heating, not a finite Q.</p>
<details><summary>Show the complete energy budget</summary>{breakdown}</details>
<p>Transport loss is stored ion-plus-electron energy divided by the selected
<strong>transport confinement time</strong>. Radiation is added separately.
Required heating offsets losses after alpha self-heating. If alpha heating exceeds
losses, the reported extra exhaust must remove the surplus to hold temperature;
this calculation does not predict stable ignition.</p>
{charts}
<p>Curves sample temperatures between the selected bounds for each density, keeping
all other inputs fixed. A wider scan does not assess whether those states can be reached.</p>
<h3>Model boundaries and power accounting</h3>
<p>Steady-state recovered heat = recovery fraction × (fusion power + deposited
auxiliary heating). Deposited heater energy is recovered once; heater inefficiency
is not recovered. Gross electricity applies the thermal conversion efficiency.
Net electricity subtracts heater wall power and the assumed plant electrical load.</p>
<p>Bremsstrahlung uses a nonrelativistic, optically thin, pure hydrogen-isotope
approximation. The full balance is limited to 0.2–30 keV; the reactivity fit alone
supports 0.2–100 keV. Other loss is a user-supplied allowance, not a radiation solver.
Density and mixture are fixed. Fuel depletion, D–D reactions, tritium production,
alpha ash, impurities, spatial profiles, magnetic stability, startup, compression,
blanket breeding and materials limits are not modelled.</p>
<h3>Published benchmark comparison</h3>{references}
<p>The independent NRL table and Bosch–Hale fit differ by up to about 25% at 1 keV;
the differences are shown rather than treated as numerical precision. For the
listed points from 5 to 100 keV, differences are below 6%.</p>
<p>Sources: <a href="{NRL_URL}#page=46">NRL 2023 fusion data</a>,
<a href="{NRL_URL}#page=59">NRL bremsstrahlung</a>,
<a href="https://doi.org/10.1088/0029-5515/32/4/i07">Bosch &amp; Hale (1992)</a>,
<a href="https://ukaea.github.io/PROCESS/physics-models/fusion_reactions/plasma_bosch_hale/">UKAEA reactivity documentation</a>,
<a href="https://ukaea.github.io/PROCESS/physics-models/plasma_confinement/">UKAEA transport accounting</a>,
<a href="https://www.iter.org/faqs">ITER explanation of plasma gain and plant power</a>.</p>
</section>'''
