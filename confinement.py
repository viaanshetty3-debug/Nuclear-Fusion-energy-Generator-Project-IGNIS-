"""Published tokamak reference scalings; SI inputs except explicitly named MW/MA.

IPB98(y,2) and mass-adjusted Martin 2008 follow UKAEA PROCESS. These empirical
relations and a few screening limits do not establish device feasibility.
"""

from dataclasses import dataclass
import math

from fusion import _number

MU0 = 1.25663706127e-6  # NIST CODATA 2022, N A^-2 (measured, not exact).
CONFINEMENT_SOURCE = 'https://ukaea.github.io/PROCESS/physics-models/plasma_confinement/'
HMODE_SOURCE = 'https://ukaea.github.io/PROCESS/physics-models/plasma_h_mode/'
DENSITY_SOURCE = 'https://github.com/ukaea/PROCESS/blob/main/process/models/physics/density_limit.py'
BETA_SOURCE = 'https://ukaea.github.io/PROCESS/physics-models/plasma_beta/plasma_beta/'


def ipb98y2_time(current_ma, field_t, density_m3, loss_power_mw, major_radius_m,
                 elongation, inverse_aspect, fuel_mass_amu, h_factor=1.0):
    """IPB98(y,2) energy timescale in seconds; power convention must be specified.

    Density is line averaged in SI; the fit uses n/(1e19 m^-3). The effective
    elongation is V/(2π² R a²). Geometric applicability is screened separately.
    """
    inputs = locals().copy()
    for name, value in inputs.items():
        _number(value, name, positive=True)
    if inverse_aspect > 1:
        raise ValueError('inverse_aspect must be at most 1')
    value = (0.0562 * h_factor * current_ma**.93 * field_t**.15 * (density_m3 / 1e19)**.41
             * loss_power_mw**(-.69) * major_radius_m**1.97 * elongation**.78
             * inverse_aspect**.58 * fuel_mass_amu**.19)
    return _number(value, 'confinement time', positive=True)


def martin_threshold_mw(density_m3, field_t, surface_m2, fuel_mass_amu):
    """Martin 2008 nominal L-H threshold with UKAEA's 2/M isotope correction."""
    for name, value in locals().copy().items():
        _number(value, name, positive=True)
    return _number(.0488 * (density_m3 / 1e20)**.717 * field_t**.803
                   * surface_m2**.941 * 2 / fuel_mass_amu, 'L-H threshold', positive=True)


def greenwald_density(current_ma, minor_radius_m):
    _number(current_ma, 'current_ma', positive=True)
    _number(minor_radius_m, 'minor_radius_m', positive=True)
    return _number(1e20 * current_ma / (math.pi * minor_radius_m**2),
                   'Greenwald density', positive=True)


@dataclass(frozen=True)
class ReferenceTokamak:
    """Illustrative large tokamak geometry, not a proposed IGNIS device design."""

    major_radius_m: float = 6.2
    minor_radius_m: float = 2.0
    elongation: float = 1.7
    field_t: float = 5.3
    current_ma: float = 15.0

    def __post_init__(self):
        for name, value in self.__dict__.items():
            _number(value, name, positive=True)
        if self.major_radius_m / self.minor_radius_m <= 2.7:
            raise ValueError('this uncorrected Martin reference requires aspect ratio R/a > 2.7')
        if self.elongation < 1:
            raise ValueError('elliptical reference requires elongation >= 1')
        _number(self.volume_m3, 'volume', positive=True)
        _number(self.surface_m2, 'surface', positive=True)

    @property
    def volume_m3(self):
        return 2 * math.pi**2 * self.major_radius_m * self.minor_radius_m**2 * self.elongation

    @property
    def surface_m2(self):
        # Pappus's theorem, with Ramanujan II approximation for ellipse perimeter.
        a, b = self.minor_radius_m, self.minor_radius_m * self.elongation
        h = ((a - b) / (a + b))**2
        perimeter = math.pi * (a + b) * (1 + 3 * h / (10 + math.sqrt(4 - 3 * h)))
        return 2 * math.pi * self.major_radius_m * perimeter
