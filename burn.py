"""Closed-volume D–D/D–T burn with one thermal temperature and adaptive RK4.

State: three reaction extents per initial ion, thermal energy per initial ion
(keV), and integrated transport/bremsstrahlung energy in the same energy units.
Inventories follow stoichiometry; no species is silently clipped or replenished.
"""

from dataclasses import asdict, dataclass
import math

from fusion import KEV_J, MAX_BALANCE_KEV, MIN_REACTIVITY_KEV, _number, dd_reactivities, dt_reactivity

# Per-reaction birth energies from NRL 2023 p44, in keV.
CHARGED_KEV = (3500.0, 4030.0, 820.0)  # DT, DD->T+p, DD->He3+n
NEUTRON_KEV = (14100.0, 0.0, 2450.0)


@dataclass(frozen=True)
class BurnConfig:
    temperature_kev: float = 10.0
    density_m3: float = 1e20
    tritium_fraction: float = 0.0
    volume_m3: float = 1.0
    confinement_s: float = 1.0
    charged_deposition: float = 1.0
    heating_mw: float = 0.0
    heater_efficiency: float = 0.5
    thermal_efficiency: float = 0.4
    heat_recovery: float = 0.9
    plant_load_mw: float = 1.0
    other_loss_mw: float = 0.0
    duration_s: float = 1.0
    max_step_s: float = 0.02
    rtol: float = 1e-7
    samples: int = 201

    def __post_init__(self):
        _number(self.temperature_kev, 'temperature_kev', MIN_REACTIVITY_KEV, MAX_BALANCE_KEV)
        for name in ('density_m3', 'volume_m3', 'confinement_s', 'duration_s', 'max_step_s'):
            _number(getattr(self, name), name, positive=True)
        for name in ('tritium_fraction', 'charged_deposition', 'thermal_efficiency', 'heat_recovery'):
            _number(getattr(self, name), name, 0, 1)
        _number(self.heater_efficiency, 'heater_efficiency', 0, 1, positive=True)
        for name in ('heating_mw', 'plant_load_mw', 'other_loss_mw'):
            _number(getattr(self, name), name)
        _number(self.rtol, 'rtol', 1e-11, 1e-2)
        if not isinstance(self.samples, int) or not 2 <= self.samples <= 10001:
            raise ValueError('samples must be an integer between 2 and 10001')
        scale = self.density_m3 * self.volume_m3 * KEV_J
        if not math.isfinite(scale) or scale <= 0:
            raise ValueError('density × volume is outside floating-point range')


class _InvalidStage(ValueError):
    """An explicit RK stage left the physical/fit domain; reject that step."""


def inventories(state, config):
    """Species counts normalized by initial total fuel ion count.

    Neutrons are a cumulative escaped count, included only in nucleon accounting.
    """
    dt, ddp, ddn = state[:3]
    return {
        'deuterium': 1 - config.tritium_fraction - dt - 2 * ddp - 2 * ddn,
        'tritium': config.tritium_fraction - dt + ddp,
        'proton': ddp,
        'helium3': ddn,
        'helium4': dt,
        'escaped_neutron': dt + ddn,
    }


def _thermal_state(state, config):
    species = inventories(state, config)
    if any(value < 0 or not math.isfinite(value) for value in species.values()):
        raise _InvalidStage('negative or nonfinite reaction inventory')
    ions = math.fsum(species[name] for name in ('deuterium', 'tritium', 'proton', 'helium3', 'helium4'))
    electrons = (species['deuterium'] + species['tritium'] + species['proton']
                 + 2 * (species['helium3'] + species['helium4']))
    capacity = 1.5 * (ions + electrons)
    temperature = state[3] / capacity
    if not MIN_REACTIVITY_KEV <= temperature <= MAX_BALANCE_KEV:
        raise _InvalidStage('temperature left the 0.2–30 keV model domain')
    return species, ions, electrons, temperature


def _terms(state, config):
    species, ions, electrons, temperature = _thermal_state(state, config)
    n = config.density_m3
    d, t = species['deuterium'], species['tritium']
    sv_ddp, sv_ddn = dd_reactivities(temperature)
    rates = (n * d * t * dt_reactivity(temperature),
             0.5 * n * d * d * sv_ddp, 0.5 * n * d * d * sv_ddn)
    charged = math.fsum(rate * energy for rate, energy in zip(rates, CHARGED_KEV))
    neutron = math.fsum(rate * energy for rate, energy in zip(rates, NEUTRON_KEV))
    transport = state[3] / config.confinement_s
    z2 = (species['deuterium'] + species['tritium'] + species['proton']
          + 4 * (species['helium3'] + species['helium4']))
    # W/m³ divided by initial ion density × J/keV -> keV per initial ion per s.
    brem = 1.69e-38 * n * electrons * z2 * math.sqrt(1000 * temperature) / KEV_J
    energy_scale = n * config.volume_m3 * KEV_J
    aux = config.heating_mw * 1e6 / energy_scale
    other = config.other_loss_mw * 1e6 / energy_scale
    derivative = (*rates, config.charged_deposition * charged + aux - transport - brem - other,
                  transport, brem)
    if not all(math.isfinite(value) for value in (*derivative, neutron, charged)):
        raise ValueError('burn rates exceed floating-point range; reduce the inputs')
    return derivative, species, ions, electrons, temperature, charged, neutron, other, aux


def _rhs(state, config):
    return _terms(state, config)[0]


def _rk4(state, step, config):
    k1 = _rhs(state, config)
    k2 = _rhs(tuple(v + step * k / 2 for v, k in zip(state, k1)), config)
    k3 = _rhs(tuple(v + step * k / 2 for v, k in zip(state, k2)), config)
    k4 = _rhs(tuple(v + step * k for v, k in zip(state, k3)), config)
    out = tuple(v + step * (a + 2 * b + 2 * c + d) / 6
                for v, a, b, c, d in zip(state, k1, k2, k3, k4))
    _thermal_state(out, config)
    return out


def snapshot(time_s, state, config):
    derivative, species, ions, electrons, temperature, charged, neutron, other, aux = _terms(state, config)
    scale = config.density_m3 * config.volume_m3 * KEV_J
    charged_energy = math.fsum(extent * energy for extent, energy in zip(state[:3], CHARGED_KEV))
    neutron_energy = math.fsum(extent * energy for extent, energy in zip(state[:3], NEUTRON_KEV))
    initial_energy = 3 * config.temperature_kev
    # Energy leaving the plasma now, rather than the steady-state Pfusion+Paux shortcut.
    outgoing_power = neutron + (1 - config.charged_deposition) * charged + derivative[4] + derivative[5] + other
    outgoing_energy = (neutron_energy + (1 - config.charged_deposition) * charged_energy
                       + state[4] + state[5] + other * time_s)
    gross_electric = config.thermal_efficiency * config.heat_recovery * outgoing_power * scale
    heater_electric = config.heating_mw * 1e6 / config.heater_efficiency
    operating_net_energy = (config.thermal_efficiency * config.heat_recovery * outgoing_energy * scale
                            - (heater_electric + config.plant_load_mw * 1e6) * time_s)
    startup = initial_energy * scale / config.heater_efficiency
    energy_residual = (state[3] - initial_energy - config.charged_deposition * charged_energy
                       - aux * time_s + state[4] + state[5] + other * time_s)
    nucleons = (2 * species['deuterium'] + 3 * species['tritium'] + species['proton']
                + 3 * species['helium3'] + 4 * species['helium4'] + species['escaped_neutron'])
    row = {
        'time_s': time_s,
        'temperature_kev': temperature,
        **{name + '_density_m3': fraction * config.density_m3 for name, fraction in species.items()},
        'ion_density_m3': ions * config.density_m3,
        'electron_density_m3': electrons * config.density_m3,
        'dt_reactions': state[0] * config.density_m3 * config.volume_m3,
        'dd_tritium_reactions': state[1] * config.density_m3 * config.volume_m3,
        'dd_helium3_reactions': state[2] * config.density_m3 * config.volume_m3,
        'dt_reactions_s': derivative[0] * config.density_m3 * config.volume_m3,
        'dd_tritium_reactions_s': derivative[1] * config.density_m3 * config.volume_m3,
        'dd_helium3_reactions_s': derivative[2] * config.density_m3 * config.volume_m3,
        'stored_energy_j': state[3] * scale,
        'fusion_power_w': (charged + neutron) * scale,
        'charged_power_w': charged * scale,
        'neutron_power_w': neutron * scale,
        'deposited_fusion_heating_w': config.charged_deposition * charged * scale,
        'aux_heating_w': config.heating_mw * 1e6,
        'transport_loss_w': derivative[4] * scale,
        'bremsstrahlung_loss_w': derivative[5] * scale,
        'other_loss_w': other * scale,
        'gross_electric_w': gross_electric,
        'heater_electric_w': heater_electric,
        'plant_load_w': config.plant_load_mw * 1e6,
        'net_electric_w': gross_electric - heater_electric - config.plant_load_mw * 1e6,
        'fusion_energy_j': (charged_energy + neutron_energy) * scale,
        'neutron_energy_j': neutron_energy * scale,
        'deposited_fusion_energy_j': config.charged_deposition * charged_energy * scale,
        'aux_energy_j': aux * time_s * scale,
        'transport_energy_j': state[4] * scale,
        'bremsstrahlung_energy_j': state[5] * scale,
        'other_loss_energy_j': other * time_s * scale,
        'outgoing_heat_energy_j': outgoing_energy * scale,
        'startup_electric_energy_j': startup,
        'operating_net_electric_energy_j': operating_net_energy,
        'net_electric_energy_j': operating_net_energy - startup,
        'energy_residual_j': energy_residual * scale,
        'charge_relative_error': electrons - 1,
        'nucleon_relative_error': nucleons / (2 + config.tritium_fraction) - 1,
    }
    if not all(math.isfinite(value) for value in row.values()):
        raise ValueError('burn output exceeds floating-point range')
    return row


@dataclass
class BurnResult:
    config: BurnConfig
    rows: list
    status: str
    accepted_steps: int
    rejected_steps: int
    max_energy_relative_error: float
    max_charge_relative_error: float
    max_nucleon_relative_error: float

    def to_dict(self):
        return asdict(self)


def simulate(config=None):
    """Integrate with RK4 step doubling, positive inventories and range termination.

    The accepted state is the two-half-step solution. /15 estimates its local
    error; no clipping, out-of-range reactivity or artificial heat control is used.
    Boundary termination records the last valid state within 1e-8 keV.
    """
    config = config if config is not None else BurnConfig()
    state = (0.0, 0.0, 0.0, 3 * config.temperature_kev, 0.0, 0.0)
    time = 0.0
    rows = [snapshot(time, state, config)]
    accepted = rejected = 0
    step = min(config.max_step_s, config.duration_s / (config.samples - 1))
    target_index = 1
    status = 'completed'
    max_energy_error = max_charge_error = max_nucleon_error = 0.0
    while time < config.duration_s:
        terms = _terms(state, config)
        species, ions, electrons, temperature = terms[1:5]
        capacity = 1.5 * (ions + electrons)
        # Only DT and the neutron-producing DD branch reduce retained ion count.
        temp_derivative = (terms[0][3] + 1.5 * temperature * (terms[0][0] + terms[0][2])) / capacity
        if temperature <= MIN_REACTIVITY_KEV + 1e-8 and temp_derivative < 0:
            status = 'temperature_lower_bound'
            break
        if temperature >= MAX_BALANCE_KEV - 1e-8 and temp_derivative > 0:
            status = 'temperature_upper_bound'
            break
        if 1 - config.tritium_fraction > 0 and species['deuterium'] <= 1e-12:
            status = 'deuterium_depleted'
            break
        if accepted + rejected >= 200000:
            raise ValueError('burn solver exceeded 200000 step attempts; shorten duration or increase max step')
        target = config.duration_s * target_index / (config.samples - 1)
        remaining = target - time
        step = min(step, config.max_step_s, remaining)
        # Land exactly on output times; a few ulps must not create a tiny
        # residual step followed by dozens of unnecessary growth steps.
        if remaining - step <= 8 * math.ulp(target):
            step = remaining
        if step <= 0 or time + step == time:
            raise ValueError('burn solver step underflow')
        try:
            coarse = _rk4(state, step, config)
            halfway = _rk4(state, step / 2, config)
            fine = _rk4(halfway, step / 2, config)
            error = max(abs(a - b) / (15 * (1e-12 + config.rtol * max(abs(v), abs(a))))
                        for v, a, b in zip(state, fine, coarse))
        except _InvalidStage:
            rejected += 1
            step *= 0.5
            continue
        if error > 1:
            rejected += 1
            step *= max(0.1, 0.9 * error ** -0.2)
            continue
        time = target if step == target - time else time + step
        state = fine
        accepted += 1
        row = snapshot(time, state, config)
        energy_scale = max(rows[0]['stored_energy_j'], row['fusion_energy_j'], row['aux_energy_j'])
        max_energy_error = max(max_energy_error, abs(row['energy_residual_j']) / energy_scale)
        max_charge_error = max(max_charge_error, abs(row['charge_relative_error']))
        max_nucleon_error = max(max_nucleon_error, abs(row['nucleon_relative_error']))
        if time == target:
            rows.append(row)
            target_index += 1
        step *= 2 if error == 0 else min(2.0, max(0.5, 0.9 * error ** -0.2))
    if rows[-1]['time_s'] != time:
        rows.append(snapshot(time, state, config))
    return BurnResult(config, rows, status, accepted, rejected, max_energy_error,
                      max_charge_error, max_nucleon_error)


def convergence_check(config):
    """Independent tighter run: halve max step and tighten relative tolerance 16×."""
    from dataclasses import replace
    coarse = simulate(config)
    fine_config = replace(config, max_step_s=min(config.max_step_s, config.duration_s / (config.samples - 1)) / 2,
                          rtol=max(1e-11, config.rtol / 16))
    fine = simulate(fine_config)
    a, b = coarse.rows[-1], fine.rows[-1]
    comparisons = {}
    for key in ('time_s', 'temperature_kev', 'deuterium_density_m3', 'tritium_density_m3',
                'fusion_energy_j', 'net_electric_energy_j'):
        absolute = abs(a[key] - b[key])
        scale = max(abs(a[key]), abs(b[key]), 1e-30)
        comparisons[key] = {'absolute_difference': absolute, 'relative_difference': absolute / scale}
    initial_energy = coarse.rows[0]['stored_energy_j']
    floors = {'time_s': 1e-8 * max(1.0, config.duration_s), 'temperature_kev': 1e-7,
              'deuterium_density_m3': 1e-10 * config.density_m3,
              'tritium_density_m3': 1e-10 * config.density_m3,
              'fusion_energy_j': 1e-9 * initial_energy,
              'net_electric_energy_j': 1e-9 * max(initial_energy, abs(b['net_electric_energy_j']))}
    for key, comparison in comparisons.items():
        comparison['tolerance'] = floors[key] + 20 * config.rtol * max(abs(a[key]), abs(b[key]))
        comparison['passed'] = comparison['absolute_difference'] <= comparison['tolerance']
    passed = coarse.status == fine.status and all(c['passed'] for c in comparisons.values())
    return coarse, {'passed': passed, 'same_termination': coarse.status == fine.status,
                    'refined_status': fine.status, 'refined_accepted_steps': fine.accepted_steps,
                    'refined_max_step_s': fine_config.max_step_s, 'refined_rtol': fine_config.rtol,
                    'comparisons': comparisons}
