# IGNIS

IGNIS is an early computational project exploring the physics needed to study
fusion energy. Prototype 02 develops the original temperature-and-ion-energy
script into a reproducible thermal distribution explorer.

The current milestone models **single-ion kinetic energies** in a classical,
nonrelativistic Maxwellian population. It does not yet calculate fusion power
or predict whether a generator can produce net energy.

## Run

Use Python 3.10 or newer. No third-party packages are required.

```sh
python ignis.py
```

This creates `ignis_results.html` in the current directory and attempts to open
it in your browser. The report is self-contained and works offline.

Choose temperatures in kelvin, select an energy threshold in keV, and export
the sampled curves:

```sh
python ignis.py --temperatures 5e7 1e8 2e8 --threshold-kev 10 \
  --output results/report.html --csv results/curves.csv --no-browser
```

`--no-browser` is useful on Termux, remote machines, and CI. Open the generated
HTML file in a browser manually. `python ignis.py --help` lists all options.
The CSV uses long format: `temperature_K`, `energy_keV`, and
`probability_density_per_keV`.

## Physics and units

For temperature `T` in kelvin and kinetic energy `E` in keV:

```text
k = 1.380649e-23 / (1000 × 1.602176634e-19) keV/K
θ = kT
f(E; T) = (2/√π) √(E/θ) exp(-E/θ) / θ
mean(E) = 3θ/2
mode(E) = θ/2
P(E ≥ E₀) = erfc(√x) + (2/√π) √x exp(-x),  x = E₀/θ
```

The density has units of inverse keV and integrates to 1 over `[0, ∞)`.
The tail expression follows by integrating this density. It describes the
fraction of ions above an arbitrary energy threshold; it is **not a fusion
probability**. Fusion calculations need relative collision energies and
reaction cross sections.

| Temperature (K) | kT (keV) | Mean ion energy (keV) |
| --- | --- | --- |
| 50,000,000 | 4.30867 | 6.46300 |
| 100,000,000 | 8.61733 | 12.9260 |
| 200,000,000 | 17.2347 | 25.8520 |

Each curve uses its own quadratic grid from 0 to 20 kT. This resolves the
low-energy region and covers more than 99.999998% of the full distribution.
HTML curves share linear axes, so widely separated temperatures can make
colder curves appear narrow. The summary uses analytic values; CSV samples
are not renormalized over the finite plotting interval.

The model assumes an isotropic thermal population with no bulk drift. It
omits compression, confinement, electron physics, losses, reactions, and fuel
burnup. Input validation checks finite positive temperatures, not the physical
applicability of the Maxwellian approximation to a particular plasma.

References:

- [NIST: SI defining constants](https://www.nist.gov/si-redefinition/meet-constants)
  for the exact Boltzmann constant and elementary charge used in conversion.
- [2023 NRL Plasma Formulary](https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf),
  especially Collisions and Transport, for the Maxwellian velocity distribution.
  The energy density above follows by changing variables using `E = mv²/2`.

## Validate

```sh
python -m unittest discover -s tests -v
```

Tests check normalization, the mean and mode, temperature scaling, the
analytic tail against independent quadrature, invalid inputs, and CLI exports.
GitHub Actions runs the suite on Python 3.10 and 3.14.

## Next milestones

1. Add a separately validated relative collision-energy and fusion-reactivity
   module with cited cross-section data and explicit validity ranges.
2. Introduce density, confinement time, and documented energy-loss models.
3. Compare heating, losses, and fusion output in a transparent energy balance.

These are planned extensions, not capabilities of the current prototype.

Licensed under Apache 2.0; see [LICENSE](LICENSE).
