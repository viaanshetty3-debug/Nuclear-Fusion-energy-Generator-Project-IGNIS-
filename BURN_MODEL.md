# Time-dependent burn model

Prototype 04 is a closed, fixed-volume thermal reaction-network model. It adds
fuel evolution and a transient energy equation to IGNIS's snapshot calculations.
It does not model separate tritium collection chambers or transfer hardware.

## Run and reproduce

```sh
python ignis.py burn --no-browser --output results/burn_report.html \
  --csv results/burn.csv --json results/burn.json
```

`python ignis.py burn --help` lists the inputs. The initial temperature flag is
`--temperature-k`, in kelvin; the API and exported configuration use kT in keV.
Defaults are kT=10 keV, n₀=1e20 m⁻³, no initial tritium, volume 1 m³,
transport confinement 1 s, no auxiliary heating and duration 1 s. Electrical
assumptions are heater efficiency 0.5, thermal efficiency 0.4, heat recovery 0.9,
and other plant load 1 MW. Charged-product energy deposition defaults to 1.

The JSON export includes the complete `BurnConfig`, samples, termination reason,
conservation diagnostics, a refined-run comparison and source references. It
can be reproduced directly:

```python
import json
from burn import BurnConfig, simulate

with open('results/burn.json') as stream:
    saved = json.load(stream)
result = simulate(BurnConfig(**saved['config']))
```

CSV includes initial configuration columns plus time samples. The offline HTML
contains temperature, fuel inventory, product buildup, fusion-power and cumulative
energy plots. Cumulative escaped neutron density is the number of neutrons that
have escaped divided by the fixed volume, not an in-plasma neutron population.

## Reaction network

| Channel | Retained charged products | Escaping product | Birth energy |
| --- | --- | --- | --- |
| D + D → T + p | T, proton | — | 1.01 + 3.02 MeV |
| D + D → He-3 + n | He-3 | neutron | 0.82 + 2.45 MeV |
| D + T → He-4 + n | He-4 | neutron | 3.5 + 14.1 MeV |

These fixed product energies come from NRL 2023, printed p44. Each D–D branch
uses its own Bosch–Hale Maxwellian reactivity. Branch fractions therefore depend
on temperature; the model does not force a 50/50 split. The coefficients are
checked against UKAEA PROCESS: its DD1 means He-3+n and DD2 means T+p.

For physical reaction-rate densities in reactions m⁻³ s⁻¹:

```text
r_DDp = 1/2 n_D² <σv>_DDp
r_DDn = 1/2 n_D² <σv>_DDn
r_DT  = n_D n_T <σv>_DT

dn_D/dt   = -2 r_DDp - 2 r_DDn - r_DT
dn_T/dt   = r_DDp - r_DT
dn_p/dt   = r_DDp
dn_He3/dt = r_DDn
dn_He4/dt = r_DT
dN_neutron_escaped/dt = (r_DDn + r_DT) V
```

The factor 1/2 counts identical D–D pairs once. It is not applied to D–T.
Three cumulative reaction extents reconstruct the five charged inventories and
the escaped-neutron count. This preserves stoichiometry; no negative species is
clipped and no fuel is replenished. Charge and nucleon totals are checked from
the reconstructed inventories, including neutrons that have left the volume.

## Temperature and thermal energy

All charged species are fully ionized and share one Maxwellian temperature with
the electrons. Products join this thermal population immediately. This is a
strong approximation: fast-product slowing, nonthermal reactions and separate
ion/electron equilibration are not calculated.

```text
n_i = n_D + n_T + n_p + n_He3 + n_He4
n_e = n_D + n_T + n_p + 2 n_He3 + 2 n_He4
U = 3/2 (n_i + n_e) kT V

P_charged = V [4.03 MeV r_DDp + 0.82 MeV r_DDn + 3.5 MeV r_DT]
P_neutron = V [2.45 MeV r_DDn + 14.1 MeV r_DT]
P_fusion = P_charged + P_neutron

P_transport = U / tau_transport
P_brem = 1.69e-38 n_e sqrt(T[eV]) [n_D+n_T+n_p+4n_He3+4n_He4] V

dU/dt = f_deposition P_charged + P_aux - P_transport - P_brem - P_other
```

The CLI's `--heating-mw` sets **constant deposited auxiliary heat**, rather than
sizing a heater to maintain a target temperature. Thermal capacity changes as
fuel reacts and helium accumulates. Temperature follows U and the retained ion
and electron counts. Radiation includes helium's charge-squared contribution.
Transport confinement removes energy without an imposed loss of particles;
particle confinement is idealized. The radiation expression is nonrelativistic
and optically thin. Fixed MeV birth partitions neglect small thermal corrections
to fusion-product energies.

No temperature thermostat or artificial exhaust is added in burn mode.
Source and loss models are evaluated only within `0.2 ≤ kT ≤ 30 keV`. The
solver terminates at the last valid state within 1e-8 keV of a boundary when the
temperature is moving outward. It also stops if the deuterium fraction falls to
1e-12 of the initial total ion density. The report records the actual elapsed
time and reason; a partial trajectory is not labelled a completed-duration burn.

## Transient electrical energy accounting

A changing thermal inventory makes the steady-state `P_fusion + P_aux` recovery
shortcut inappropriate. Instead, the model tracks actual outgoing power:

```text
P_out = P_neutron + (1-f_deposition) P_charged + P_transport + P_brem + P_other
P_gross_electric = eta_thermal eta_recovery P_out
P_net_electric = P_gross_electric - P_aux/eta_heater - P_plant

U(t) - U(0) + integral(P_out dt) = integral((P_fusion + P_aux) dt)
E_startup_electric = U(0)/eta_heater
E_net_electric(t) = integral(P_net_electric dt) - E_startup_electric
```

`operating_net_electric_energy_j` omits the initial startup charge;
`net_electric_energy_j` includes it and starts negative. The startup allowance
is a minimum thermal-inventory estimate, not a startup simulation. Ionization,
ramp losses, pumps/magnets during preparation and other equipment loads would
add costs. Remaining plasma thermal energy is not automatically credited at the
end. Thus recovering the initial plasma's stored heat cannot create apparent
positive cumulative net energy in a zero-fusion cooling test.

Heat recovery and all conversion efficiencies are prescribed effective inputs.
The model does not establish whether such recovery is achievable or whether
any particular confinement device can maintain the assumed timescale.

## Integration and convergence

The standard-library solver uses RK4 with step doubling. Each attempted interval
is evaluated once as a full step and once as two half steps. The two-half-step
solution is accepted when its estimated local error `(fine - coarse)/15` fits
`1e-12 + rtol × max(abs(old), abs(fine))` for every normalized state component.
Reaction extents are normalized by initial ion count and energy by initial ion
count × keV. Invalid intermediate inventories or temperatures reject the step
and reduce it; there is no out-of-domain reactivity extrapolation.

The maximum step, error tolerance and requested output sampling are separate.
Output-time alignment accounts for floating-point roundoff, avoiding spurious
near-zero steps. The convergence run uses a maximum step no larger than half
both the original maximum step and the output interval. It tightens rtol by up
to 16×, down to a floor of 1e-11, using the same requested output times.

The report compares final elapsed time, temperature, D/T inventories, fusion
energy and net electricity. Allowed differences are `20 × rtol` times the
quantity scale plus an absolute floor: 1e-8 × max(1, requested seconds) for time,
1e-7 keV for temperature, 1e-10 n₀ for species, and energy floors of order 1e-9
of the initial thermal or cumulative electrical scale. Both runs must have the
same termination reason. A failed convergence check produces reviewable outputs
and a nonzero exit code. Local error control alone is not a global accuracy bound.

## Verification

```sh
python -m unittest discover -s tests -v
```

Tests check the full suite of snapshot calculations plus:

- D–D identical-reactant counting and channel-specific product energies.
- Tritium formation from initially pure D, then nonzero D–T consumption.
- Independent charge/nucleon reconstruction and nonnegative inventories.
- Whole-system transient energy closure, including escaping neutron energy.
- Startup charging and a zero-fusion cooling case that cannot yield net energy.
- Changing thermal capacity and helium contributions to bremsstrahlung.
- An analytic transport-plus-bremsstrahlung cooling solution.
- Analytic constant-reactivity D–T and D–D depletion solutions.
- Fourth-order error reduction when the maximum step is halved.
- A genuinely finer convergence run and temperature-bound termination.
- Reproducibility and finite CSV/JSON/HTML outputs.

`data/nrl_dd_reference.csv` contains the independent NRL 2023 summed D–D
reactivities. The fitted branch sum differs by about 34% at 1 keV, 16% at 2 keV,
and less than 10% at the listed 5–100 keV points. Tests use explicit 40%, 20% and
10% comparison limits respectively, and the report displays actual differences.
These broad source comparisons do not quantify fit uncertainty or validate the
thermalization, confinement or electricity assumptions.

## Limits

The network intentionally includes only the two D–D branches and D–T. D–He3,
T–T and other secondary reactions are omitted; accumulated helium does not burn
or leave. There is no refuelling, particle transport, ash removal, volume change,
compression, separate electron/ion energy evolution, spatial structure, magnetic
stability, breeding, material limit or reactor equipment model. The simulated
D–D-to-D–T chain is a shared thermal population, not a tritium separation process.

References are linked in [README.md](README.md#sources). Relevant NRL locations
are printed p44 (products), p45 (reactivities), and p58 Eq30 (radiation). UKAEA's
official implementation confirms the Bosch–Hale coefficient/channel mapping.
