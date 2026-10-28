# Confinement and break-even reference study

Prototype 05 asks how a finite grid of uniform D–T target states behaves when
confinement is coupled to a published empirical scaling. It uses a **tokamak
reference**, chosen as a documented comparison case. It does not establish a
physical confinement mechanism for the original IGNIS D–D fuel-cycle concept.

## Run

```sh
python ignis.py study --no-browser --output results/feasibility_report.html \
  --csv results/feasibility.csv --json results/feasibility.json
```

The offline report has selectors for H and fuel mixture, a temperature/density
map, eligible-case rankings, break-even crossings, exclusion reasons and
sensitivity comparisons. Every output path is optional except the HTML report.
Use `python ignis.py study --help` for the complete input list.

Default grid:

- kT = 5, 7.5, 10, 12.5, 15, 20, 25, 30 keV.
- Total D+T ion density = 0.5, 0.75, 1.0, 1.15 × 10²⁰ m⁻³.
- Confinement multiplier H = 0.8, 1.0, 1.1, 1.2.
- Tritium number fraction = 0.1, 0.25, 0.5.

These 384 cases assume equal ion/electron temperatures and uniform profiles,
so line- and volume-averaged densities are approximated as equal. Reactions
use the existing D–T Maxwellian reactivity module. This study does not add D–D
power or tritium production; use `ignis.py burn` for that limited evolving network.

## Reference geometry

The illustrative device has major radius R=6.2 m, minor radius a=2 m,
effective elongation κ=1.7, toroidal field B=5.3 T and current I=15 MA. It is a
large benchmark geometry, not a proposed IGNIS hardware configuration.

For an elliptical torus:

```text
V = 2π² R a² κ
S = 2π R × ellipse_perimeter(a, κa)
```

The surface uses Ramanujan's second ellipse-perimeter approximation.
κ is consistent with the IPB definition `V/(2π² R a²)`. The uncorrected
Martin reference is restricted to `R/a > 2.7`; low-aspect-ratio corrections
are outside this implementation. Geometry and fit inputs must be finite and
positive, with κ ≥ 1. This restriction is not a claim that all remaining
geometries lie inside the experiments used to train either empirical fit.

## Coupled confinement and power

The implemented IPB98(y,2) relation is:

```text
tau = 0.0562 H I^0.93 B^0.15 n19^0.41 P_MW^-0.69
      × R^1.97 κ^0.78 (a/R)^0.58 M^0.19
```

Units are seconds, MA, T, density in 10¹⁹ m⁻³, MW, metres and average fuel mass
number `M=2+f_T`. H multiplies the empirical relation. The implementation is
checked against UKAEA PROCESS and its published unit-input test fixture.

At a specified temperature and density, the existing ion-plus-electron model
gives thermal energy U. Writing `tau = A P_MW^-0.69`, the coupled power is solved
analytically instead of treating confinement and heating as independent inputs:

```text
U_MJ = tau P_MW = A P_MW^0.31
P_MW = (U_MJ/A)^(1/0.31)
```

Because empirical confinement scalings have uncertain radiation treatment,
two conventions are supported, following the choices documented in PROCESS:

| Convention | Power used by the scaling | Total loss | Transport power |
| --- | --- | --- | --- |
| `transport` (default) | Transport power | Scaled power + core radiation | Scaled power |
| `total` | Total loss | Scaled power | Scaled power − core radiation |

Bremsstrahlung plus the user-specified extra core-radiation allowance form
`P_rad`. Radiation is counted once in either convention. Under `total`, a
nonpositive transport power is excluded. `scaling_confinement_s` and
`transport_confinement_s` are exported separately to prevent conflating them.

The signed required auxiliary heating is:

```text
P_aux = P_transport + P_rad - P_alpha_deposited
```

Negative P_aux is excluded as a self-heating surplus with no fixed point at the
requested target state and prescribed scaling. The code does not clamp it to
zero or add an unmodelled exhaust controller. It does not imply that the device
could not reach some other temperature or confinement condition.

## Implemented screens

The following are necessary reference screens and analyst-selected limits,
not a complete stability or feasibility assessment:

1. **Heating capacity:** deposited P_aux ≤ 200 MW by default, and P_aux ≥ 0.
2. **Density:** `n_G = 1e20 I_MA/(π a²)` m⁻³; `n/n_G ≤ 1` by default.
3. **Thermal pressure:** `p = 2 n_i kT`, `beta_t = 2 mu0 p/B²`, and
   `beta_N = 100 beta_t a B/I_MA ≤ 3` by default. This omits fast-particle pressure
   and does not solve a magnetic equilibrium or MHD stability problem.
4. **H-mode applicability:** transported power must be at least 1.3 times the
   nominal L-H threshold, using UKAEA's mass-adjusted Martin 2008 relation:
   `P_LH[MW] = 0.0488 n20^0.717 B^0.803 S^0.941 × (2/M)`.
   Crossing this threshold alone does not guarantee an H-mode state.
5. **Radiation fraction:** core radiation / total loss ≤ 0.3 by default. The 0.3
   ceiling is a chosen screening restriction for the simple radiation model,
   not a published universal boundary or an empirical-fit accuracy guarantee.
6. **Model scope:** kT is restricted to 0.2–30 keV by the existing radiation
   model. No reactivity or loss extrapolation is performed outside that range.

The density/pressure ceilings, heating capacity, H-mode margin and radiation
ceiling are explicit configuration inputs. Exclusion reasons may overlap.
Out-of-screen cases retain diagnostic algebraic budgets but have null eligible
net outputs. They cannot win rankings, count as positive cases, or bracket a
reported break-even crossing. CSV uses blank cells for null values.

## Electrical break-even

At a consistent steady target state:

```text
P_recovered = eta_recovery (P_fusion + P_aux)
P_gross_electric = eta_thermal P_recovered
P_net_electric = P_gross_electric - P_aux/eta_heater - P_plant_load
E_startup = startup_multiplier × U/eta_heater
E_flat_top_net = P_net_electric × duration - E_startup
startup_payback_time = E_startup/P_net_electric, only if P_net_electric > 0
```

The default flat-top duration is 10 s; heater efficiency is 0.5, thermal
conversion 0.4, recovery 0.9 and other plant electrical load 100 MW. These are
illustrative assumptions, not inferred engineering performance. The startup
multiplier cannot be below 1; actual preparation, ramp and equipment costs
can exceed this minimum thermal allowance.

The flat-top calculation assumes the target temperature, density and mixture
are maintained. It is **not** a trajectory from the closed-volume burn solver,
and does not account for the mechanisms or energy costs of maintaining fuel
composition, current drive, ash removal or temperature control.

A case can have positive operating power but insufficient duration to repay
startup. Bracketed flat-top-energy zeros are refined by bisection between
adjacent eligible temperature samples at fixed density, H and mixture. If a
midpoint becomes excluded, that crossing is discarded. The finite grid can
miss roots and cannot prove global optimality or the absence of any operating
window outside the samples.

## Sensitivity study

For each baseline case that passes the screens, the study evaluates all binary
corners of six selected inputs: H, thermal efficiency, heat recovery, heater
efficiency, plant load and startup multiplier. The default excursions are ±20%.
Efficiencies are bounded above by 1 and startup cost below by the minimum thermal
allowance. Both radiation conventions are tested, for 128 corner evaluations.
Some evaluations can coincide when an excursion is bounded or zero.

The report retains the number of eligible corners, positive corners, rejection
reasons and min/max energy **among eligible corners only**. No admissible corner
means the range is undefined, not zero. An excluded corner fails the all-corner
criterion. The highest-ranked baseline also gets a separate one-at-a-time table,
so a change in energy can be distinguished from leaving the screened regime.

These are deterministic scenario checks, not probability distributions or
confidence intervals. Passing all corners does not guarantee interior behavior,
and failing them does not prove that a controller could not shift to a different
target state. In the default study, some baseline points are positive, but no
individual sampled state passes and stays positive at every selected corner.

## Reproduce and verify

```python
import json
from feasibility import StudyConfig, StudyGrid, run_study

with open('results/feasibility.json') as stream:
    saved = json.load(stream)
reproduced = run_study(StudyConfig(**saved['config']), StudyGrid(**saved['grid']))
```

The JSON includes every case, inputs, uncertainty summaries, roots and
highest-ranked baseline. The CSV includes reference inputs, budgets and corner
summaries. Python 3.10+ and only the standard library are required.

```sh
python -m unittest discover -s tests -v
```

Tests cover the official PROCESS IPB98 and Martin numerical fixtures, unit
conversion, coupled power/scaling closure, both radiation conventions,
startup accounting, all implemented screens, exclusion-aware rankings and roots,
sensitivity transitions, a case that passes the all-corner check, invalid input
handling, finite exports and JSON reproducibility. Earlier thermal/burn tests
remain in the suite.

## Sources and limits

- [UKAEA IPB98 and radiation conventions](https://ukaea.github.io/PROCESS/physics-models/plasma_confinement/).
- [Official confinement implementation](https://github.com/ukaea/PROCESS/blob/main/process/models/physics/confinement_time.py)
  and [numeric fixtures](https://github.com/ukaea/PROCESS/blob/main/tests/unit/models/physics/test_confinement_time.py).
- [UKAEA L-H threshold and isotope correction](https://ukaea.github.io/PROCESS/physics-models/plasma_h_mode/)
  and [Martin numeric fixture](https://github.com/ukaea/PROCESS/blob/main/tests/unit/models/physics/test_l_h_transition.py).
- [UKAEA Greenwald implementation](https://github.com/ukaea/PROCESS/blob/main/process/models/physics/density_limit.py).
- [UKAEA beta definitions and ceilings](https://ukaea.github.io/PROCESS/physics-models/plasma_beta/plasma_beta/).
- [ITER reference discussion of H-mode power margin](https://www-pub.iaea.org/MTCD/Publications/PDF/ITER-EDA-DS-20.pdf).
- [NIST CODATA constants](https://physics.nist.gov/cuu/Constants/Table/allascii.txt), including vacuum permeability.

Published relations are empirical and regime-dependent. This code does not
verify inclusion in the original training datasets or validate extrapolation to
an arbitrary device. Profiles, equilibrium, safety factor, MHD stability,
nonthermal pressure, tritium supply/production, extra reaction channels, detailed
radiation, divertor heat fluxes, neutron transport, materials, magnets and plant
engineering remain unmodelled. The selected plant load is a lumped allowance,
not a current-drive or facility power calculation. This is a constrained
comparison study, not a validated IGNIS reactor design.
