# IGNIS

IGNIS Prototype 04 models **time-dependent D–D/D–T fuel burns** and retains the
steady-state D–T energy-balance explorer. It tracks deuterium depletion, tritium
production and consumption, helium/proton buildup, evolving temperature, fusion
output, and plant electrical energy.

A positive calculated net output is conditional on the specified plasma,
confinement, exhaust and conversion assumptions being achievable. This project
is a computational model, not a demonstrated generator or a reactor-design tool.

## Time-dependent burn

Python 3.10 or newer is sufficient; no third-party packages are required.

```sh
python ignis.py burn --no-browser --output results/burn_report.html \
  --csv results/burn.csv --json results/burn.json
```

Open `results/burn_report.html`. The default is a 1-second, initially pure-deuterium
burn at kT=10 keV, total ion density 1e20 m⁻³, fixed volume 1 m³ and transport
confinement time 1 second, with no auxiliary heating. Both D–D branches produce
tritium/protons and helium-3/neutrons; produced tritium can react through D–T.

To explore a seeded mixture and fixed auxiliary heating:

```sh
python ignis.py burn --temperature-k 1e8 --density-m3 1e20 \
  --tritium-fraction 0.1 --heating-mw 0.1 --duration-s 2 \
  --output results/seeded_burn.html --json results/seeded_burn.json --no-browser
```

The report includes five trajectory plots, charge/nucleon/energy residuals,
and a tighter integration run for numerical convergence. Temperature evolves
without a thermostat. Runs terminate at the last valid state near either
0.2 or 30 keV instead of extrapolating. The CLI reports the actual simulated
interval and termination reason. Cumulative net electricity subtracts a minimum
initial thermal startup cost; it starts negative.

In the default illustrative run, fusion releases about **1.19 kJ**, temperature
falls to **3.52 keV**, and net electricity including startup is about **−1.85 MJ**.
These are model outputs under the stated assumptions, not measured performance.

See [BURN_MODEL.md](BURN_MODEL.md) for the reaction network, transient energy
accounting, solver tolerances, validation and limitations. Use
`python ignis.py burn --help` for all inputs. To reproduce an export, pass the
`config` object in its JSON file to `burn.simulate(BurnConfig(**config))`.

## Steady-state and thermal modes

Python 3.10 or newer is sufficient. There are no third-party dependencies.

```sh
python ignis.py --no-browser
```

Open `ignis_results.html` in a browser. It is self-contained and works offline.
Without `--no-browser`, the script also attempts to open the browser automatically.
The default report contains nine cases: three temperatures × three densities,
thermal ion-energy curves, fusion and net-power plots, energy budgets, and a
published reference comparison.

Export a complete, reproducible scenario:

```sh
python ignis.py --temperatures 5e7 1e8 2e8 --densities 5e19 1e20 2e20 \
  --volume-m3 1 --tritium-fraction 0.5 --confinement-s 1 \
  --alpha-deposition 1 --heater-efficiency 0.5 \
  --thermal-efficiency 0.4 --heat-recovery 0.9 --plant-load-mw 1 \
  --output results/fusion_report.html --csv results/curves.csv \
  --fusion-csv results/cases.csv --json results/scenario.json --no-browser
```

Inputs in `--temperatures` are **kelvin**; the reaction model uses `kT` in keV.
Density is **total fuel ion density**, `n_D + n_T`, in m⁻³. Each temperature is
combined with every density. Use `python ignis.py --help` for all flags.

The original thermal-only mode remains available, including temperatures outside
the full D–T balance's range:

```sh
python ignis.py --thermal-only --temperatures 1000 10000 --no-browser
```

Fusion-specific options cannot be combined with `--thermal-only`.

## Steady-state scenario assumptions

| Input | Default | Meaning |
| --- | --- | --- |
| Temperatures | 5e7, 1e8, 2e8 K | Ion and electron temperatures are equal |
| Densities | 5e19, 1e20, 2e20 m⁻³ | Total D+T ion density; electron density is equal |
| Volume | 1 m³ | Uniform plasma volume |
| Tritium fraction | 0.5 | Number fraction of tritium in the fuel ions |
| Confinement | 1 s | **Transport** energy timescale, excluding radiation |
| Alpha deposition | 1 | Fraction of alpha power deposited in the plasma |
| Heater efficiency | 0.5 | Electricity-to-deposited-heat efficiency |
| Thermal efficiency | 0.4 | Recovered-heat-to-electricity efficiency |
| Heat recovery | 0.9 | Fraction of outgoing thermal power recovered |
| Plant load | 1 MW | Electrical demand beyond the plasma heater |
| Other plasma loss | 0 MW | User-specified extra loss, `--other-loss-mw` |

These are illustrative inputs, not inferred device performance. Fractions and
conversion efficiencies must be in `[0, 1]`; heater efficiency must be positive.
Volume, density and confinement time must be finite and positive. Both output
signs are supported; losses and loads are never silently clipped away.

## Physics

### Reaction rates and fusion power

The pure D–T model uses Maxwellian **relative-velocity averaged reactivity** from
the Bosch–Hale fit. This already averages over thermal collisions; an arbitrary
single-ion energy threshold is not used to determine which ions fuse.

```text
n_T = f_T n_i                    n_D = (1 - f_T) n_i
rate density = n_D n_T <σv>      [reactions m⁻³ s⁻¹]
reactions/s = rate density × V
P_alpha = reactions/s × 3.5 MeV
P_neutron = reactions/s × 14.1 MeV
P_fusion = P_alpha + P_neutron
```

D and T are distinct reactants, so there is no identical-particle factor of 1/2.
Internally, reactivity is m³/s, energies are joules and powers are watts.
`1 eV = 1.602176634e-19 J` exactly.

The reactivity function accepts `0.2 ≤ kT ≤ 100 keV`, with no extrapolation.
The full energy balance is deliberately restricted to `0.2 ≤ kT ≤ 30 keV` to
limit the scope of the simple nonrelativistic radiation approximation. The 30 keV
cap is a modelling choice, not a published accuracy guarantee.

### Heating, transport, and radiation

For a fully ionized, uniform pure D–T plasma with `n_e = n_i` and `T_e = T_i`:

```text
U = (3/2) (n_i + n_e) kT V = 3 n_i kT V
P_transport = U / tau_transport
P_bremsstrahlung = 1.69e-38 n_e n_i sqrt(T_e[eV]) V   [W]
P_losses = P_transport + P_bremsstrahlung + P_other
P_alpha_deposited = alpha_deposition × P_alpha
P_aux_required = max(0, P_losses - P_alpha_deposited)
P_extra_exhaust = max(0, P_alpha_deposited - P_losses)
```

Bremsstrahlung is the NRL hydrogen-like, optically thin expression converted to
SI with ionic charge `Z=1`. The confinement input covers **transport only**;
adding radiation to a confinement time that already includes radiation would
count that loss twice. `P_other` is a user allowance, not a separate prediction
of cyclotron or impurity radiation.

Auxiliary heating is automatically sized to the minimum required deposited
power. If alpha heating exceeds the baseline losses, holding the requested
temperature requires the reported additional exhaust. The model accounts for
that controlled heat removal; it does not demonstrate stable ignition or solve
for an equilibrium temperature. `U` gives the thermal energy inventory; startup
heating, ionization and compression work are not calculated.

### Plant power and net electricity

```text
P_alpha_deposited + P_aux_required = P_losses + P_extra_exhaust
P_recovered_heat = heat_recovery × (P_fusion + P_aux_required)
P_gross_electric = thermal_efficiency × P_recovered_heat
P_heater_electric = P_aux_required / heater_efficiency
P_net_electric = P_gross_electric - P_heater_electric - P_plant_load
Q_plasma = P_fusion / P_aux_required    (only when P_aux_required > 0)
```

At steady state, neutron power, escaping alpha power, plasma losses and extra
exhaust sum to `P_fusion + P_aux_required`. The heat cycle recovers a fraction of
that sum. Deposited auxiliary heat is counted once; heater inefficiency is not
recovered. Do **not** subtract plasma losses again from gross fusion power when
calculating plant electricity: here they are outgoing heat, subject to the heat
recovery factor. The electrical plant load is subtracted separately.

Plasma gain `Q > 1` does not imply positive net electrical output. When required
auxiliary heating is zero, the zero-denominator Q is displayed as a dash,
exported as JSON `null`, and left blank in CSV. Extra exhaust is still reported.
All electrical outputs assume the prescribed state and required heating/exhaust
can be maintained. Blanket multiplication and breeding energy are omitted.

### Original thermal distribution

The thermal panel retains the single-ion kinetic-energy model:

```text
k = 1.380649e-23 / (1000 × 1.602176634e-19) keV/K
θ = kT
f(E; T) = (2/√π) √(E/θ) exp(-E/θ) / θ
mean(E) = 3θ/2                  mode(E) = θ/2
P(E ≥ E₀) = erfc(√x) + (2/√π) √x exp(-x),  x = E₀/θ
```

Its density is in inverse keV and integrates to 1 over `[0, ∞)`. Each curve has
an independent quadratic grid through 20 kT; summary values are analytic.
The thermal tail statistic is separate from the fusion reaction calculation.

## Validation and published comparisons

```sh
python -m unittest discover -s tests -v
```

The suite also checks the time-dependent network against analytic cooling and
constant-reactivity depletion solutions, fourth-order step refinement, charge
and nucleon conservation, transient energy closure, startup accounting,
nonnegative inventories and model-boundary termination.

The steady-state suite checks thermal normalization and moments, independent published
reactivity benchmarks, cgs-to-SI radiation conversion, reaction energies,
mixture/density/volume/confinement scaling, plasma and plant energy conservation,
heating deficits and surpluses, negative electrical balances, range rejection,
and reproducibility of CLI/CSV/JSON exports. CI runs on Python 3.10 and 3.14.

`data/nrl_dt_reference.csv` transcribes seven independent D–T reference values
from **NRL 2023, printed p45**, retaining the original cm³/s units. The report
shows each value alongside the Bosch–Hale result after SI conversion. The
largest difference is about **24.7% at 1 keV**, then 14.5% at 2 keV; differences
at the listed 5–100 keV points are below 6%. Tests use explicit comparison bounds
of 30%, 20%, and 6% respectively. These bounds are cross-checks against a
separate approximate table, not uncertainty estimates or a precision claim.
No validation of confinement, plant efficiency, stability or achievable net
power is implied by agreement of reaction rates.

## Outputs and code

- `ignis.py`: CLI, thermal distribution and standalone report assembly.
- `fusion.py`: Bosch–Hale reactivity, validated inputs and steady-state accounting.
- `fusion_report.py`: density/temperature plots, budgets, reference comparison and exports.
- `burn.py`: closed-volume burn integration and convergence checks.
- `burn_cli.py` / `burn_report.py`: burn-mode CLI, trajectory reports and exports.
- `data/nrl_dt_reference.csv` / `data/nrl_dd_reference.csv`: independent benchmarks.
- `--csv`: thermal samples, `temperature_K`, `energy_keV`, `probability_density_per_keV`.
- `--fusion-csv`: one row per case, including all plant inputs and numerical results.
- `--json`: versioned scenario inputs, assumptions, all results and reference comparison.

Generated files are excluded from Git; recreate them with the commands above.

## Limits and future work

Steady-state mode prescribes density, mixture and temperature. Burn mode evolves
the fuel inventory and a shared temperature in a closed, fixed volume; confinement
and electrical conversion remain assumed inputs. Burn mode retains charged
products with instantaneous thermalization and includes only D–T and the two
D–D branches. Helium-3 and helium-4 accumulate without further burning or removal.

Neither mode predicts spatial profiles, compression, achievable confinement,
magnetic/inertial stability, fast-ion slowing, breeding, material limits or
specific equipment performance. The burn startup cost covers the initial thermal
inventory and heater efficiency; ionization, ramp losses and other startup loads
would add costs. Positive numerical output alone does not establish feasibility.

The original D–D-to-D–T chain is now represented as a limited thermal reaction
network: the tritium-producing branch is `D + D → T + proton`; the other main
branch produces `He-3 + neutron`. This does not simulate collecting tritium and
transferring it between separate devices.

Assessing a particular device would require choosing and validating a
confinement model, particle transport and ash removal, separate ion/electron
energy equations, additional reaction channels, and engineering constraints.

## Sources

- [Bosch & Hale (1992), Nuclear Fusion 32, 611–631](https://doi.org/10.1088/0029-5515/32/4/i07):
  thermal reactivity parametrization and coefficients.
- [UKAEA PROCESS reactivity documentation](https://ukaea.github.io/PROCESS/physics-models/fusion_reactions/plasma_bosch_hale/)
  and [official implementation](https://github.com/ukaea/PROCESS/blob/main/process/models/physics/fusion_reactions.py):
  independent check of the D–T and D–D branch coefficients, fit range and cm³/s-to-m³/s conversion.
- [2023 NRL Plasma Formulary](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf):
  fusion products on printed p44, reference reactivities on p45, bremsstrahlung
  on p58, Eq30. PDF viewer page numbers are one greater than printed numbers.
- [UKAEA confinement accounting](https://ukaea.github.io/PROCESS/physics-models/plasma_confinement/)
  and [core power balance](https://ukaea.github.io/PROCESS/physics-models/plasma_power_balance/):
  energy-over-confinement transport losses and deposited-alpha power accounting.
- [ITER FAQs](https://www.iter.org/faqs): distinction between plasma gain and plant net electricity.
- [NIST SI constants](https://www.nist.gov/si-redefinition/meet-constants): exact SI conversion constants.

Licensed under Apache 2.0; see [LICENSE](LICENSE).
