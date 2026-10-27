"""Uniform D–T fusion and steady-state power accounting, in SI units.

Reactivity: Bosch & Hale (1992), checked against UKAEA PROCESS and NRL 2023.
This is a prescribed plasma scenario, not a confinement or reactor-design solver.
"""

from dataclasses import asdict, dataclass
import math

EV_J = 1.602176634e-19
KEV_J = 1000 * EV_J
ALPHA_ENERGY_J = 3.5e6 * EV_J
NEUTRON_ENERGY_J = 14.1e6 * EV_J
MIN_REACTIVITY_KEV = 0.2
MAX_REACTIVITY_KEV = 100.0
MAX_BALANCE_KEV = 30.0  # Chosen conservative scope for the nonrelativistic loss model.
NRL_URL = "https://homepage.physics.uiowa.edu/~ghowes/teach/phys4731/resources/NRL_Plasma_Formulary_2023.pdf"
BOSCH_HALE_DOI = "https://doi.org/10.1088/0029-5515/32/4/i07"


def _number(value, name, minimum=0.0, maximum=None, *, positive=False):
    if (not math.isfinite(value) or value < minimum
            or (positive and value == 0) or (maximum is not None and value > maximum)):
        interval = f"[{minimum:g}, {maximum:g}]" if maximum is not None else f">= {minimum:g}"
        if positive:
            interval = f"in (0, {maximum:g}]" if maximum is not None else "> 0"
        raise ValueError(f"{name} must be finite and {interval}")
    return value


def dt_reactivity(temperature_kev):
    """Maxwellian relative-velocity averaged <sigma v> in m³/s, 0.2–100 keV.

    Bosch–Hale thermal reactivity parametrization, D–T coefficients.
    Input is kT, not the single-ion mean energy 3kT/2. No extrapolation.
    """
    t = _number(temperature_kev, "D–T temperature (keV)", MIN_REACTIVITY_KEV,
                MAX_REACTIVITY_KEV)
    correction = t * (1.51361e-2 + t * (4.60643e-3 - t * 1.06750e-4)) / (
        1 + t * (7.51886e-2 + t * (1.35000e-2 + t * 1.36600e-5)))
    theta = t / (1 - correction)
    xi = (34.3827**2 / (4 * theta)) ** (1 / 3)
    return 1e-6 * 1.17302e-9 * theta * math.sqrt(xi / (1.124656e6 * t**3)) * math.exp(-3 * xi)


@dataclass(frozen=True)
class PlantConfig:
    """Pure, fully ionized D–T; Te=Ti. All efficiencies are assumed inputs.

    confinement_s is the TRANSPORT energy timescale; radiation is separate.
    other_loss_mw is total extra plasma loss, and plant_load_mw is electrical.
    """

    volume_m3: float = 1.0
    tritium_fraction: float = 0.5
    confinement_s: float = 1.0
    alpha_deposition: float = 1.0
    heater_efficiency: float = 0.5
    thermal_efficiency: float = 0.4
    heat_recovery: float = 0.9
    plant_load_mw: float = 1.0
    other_loss_mw: float = 0.0

    def __post_init__(self):
        for name in ("volume_m3", "confinement_s"):
            _number(getattr(self, name), name, positive=True)
        for name in ("tritium_fraction", "alpha_deposition", "thermal_efficiency", "heat_recovery"):
            _number(getattr(self, name), name, 0, 1)
        _number(self.heater_efficiency, "heater_efficiency", 0, 1, positive=True)
        for name in ("plant_load_mw", "other_loss_mw"):
            _number(getattr(self, name), name)


@dataclass(frozen=True)
class FusionCase:
    temperature_kev: float
    density_m3: float
    reactivity_m3_s: float
    deuterium_density_m3: float
    tritium_density_m3: float
    reaction_rate_m3_s: float
    reactions_s: float
    fusion_power_w: float
    alpha_power_w: float
    neutron_power_w: float
    alpha_heating_w: float
    stored_energy_j: float
    transport_loss_w: float
    bremsstrahlung_loss_w: float
    other_loss_w: float
    total_loss_w: float
    required_heating_w: float
    extra_exhaust_w: float
    self_heating_margin_w: float
    plasma_gain_q: float | None
    recovered_heat_w: float
    gross_electric_w: float
    heater_electric_w: float
    plant_load_w: float
    net_electric_w: float

    def to_dict(self):
        return asdict(self)


def evaluate_case(temperature_kev, density_m3, config=None):
    """Evaluate a steady target state; additional heat exhaust closes a surplus.

    Required auxiliary heating and additional exhaust are nonnegative. A surplus
    is not evidence of a stable ignition point. Net electricity assumes the
    prescribed state, heating, heat recovery and exhaust can actually be supplied.
    """
    config = config if config is not None else PlantConfig()
    t = _number(temperature_kev, "balance temperature (keV)", MIN_REACTIVITY_KEV,
                MAX_BALANCE_KEV)
    n = _number(density_m3, "ion density (m^-3)", positive=True)
    try:
        sv = dt_reactivity(t)
        nd = n * (1 - config.tritium_fraction)
        nt = n * config.tritium_fraction
        # D and T are distinct reactants: no identical-particle factor of 1/2.
        rate = nd * nt * sv
        reactions = rate * config.volume_m3
        alpha = reactions * ALPHA_ENERGY_J
        neutron = reactions * NEUTRON_ENERGY_J
        fusion = alpha + neutron
        alpha_heat = config.alpha_deposition * alpha
        # Ion plus electron thermal energy, ne=ni=n and Te=Ti.
        stored = 3 * n * t * KEV_J * config.volume_m3
        transport = stored / config.confinement_s
        # NRL 2023 p58 Eq30, converted from cm^-3 and W/cm³ to SI; Z=1.
        brem = 1.69e-38 * n * n * math.sqrt(1000 * t) * config.volume_m3
        other = config.other_loss_mw * 1e6
        losses = transport + brem + other
        required = max(0.0, losses - alpha_heat)
        extra_exhaust = max(0.0, alpha_heat - losses)
        margin = alpha_heat - losses
        q = fusion / required if required > 0 else None
        # Steady state: escaping neutrons + alphas + plasma losses + extra exhaust
        # sum to Pfusion + Paux. Recover deposited auxiliary heating ONCE.
        recovered = config.heat_recovery * (fusion + required)
        gross_electric = config.thermal_efficiency * recovered
        heater_electric = required / config.heater_efficiency
        plant_load = config.plant_load_mw * 1e6
        net = gross_electric - heater_electric - plant_load
        case = FusionCase(t, n, sv, nd, nt, rate, reactions, fusion, alpha, neutron,
                          alpha_heat, stored, transport, brem, other, losses, required,
                          extra_exhaust, margin, q, recovered, gross_electric,
                          heater_electric, plant_load, net)
    except OverflowError as exc:
        raise ValueError("scenario exceeds floating-point range; reduce its inputs") from exc
    if any(value is not None and not math.isfinite(value) for value in asdict(case).values()):
        raise ValueError("scenario exceeds floating-point range; reduce its inputs")
    if stored <= 0 or brem <= 0 or transport <= 0:
        raise ValueError("scenario is below floating-point resolution; increase its inputs")
    return case
