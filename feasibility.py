"""Confinement-coupled D–T operating-point screening and deterministic sensitivity.

The study holds uniform target states and fuel mix fixed; it is separate from
closed-volume burn evolution. Every ranking and break-even search filters limits.
"""

from collections import Counter
from dataclasses import asdict, dataclass, field, replace
from itertools import product
import math

from confinement import MU0, ReferenceTokamak, greenwald_density, ipb98y2_time, martin_threshold_mw
from fusion import KEV_J, PlantConfig, _number, evaluate_case


@dataclass(frozen=True)
class StudyConfig:
    device: ReferenceTokamak = field(default_factory=ReferenceTokamak)
    radiation_convention: str = 'transport'
    alpha_deposition: float = 1.0
    heater_efficiency: float = .5
    thermal_efficiency: float = .4
    heat_recovery: float = .9
    plant_load_mw: float = 100.0
    extra_radiation_mw: float = 0.0
    heating_capacity_mw: float = 200.0
    greenwald_fraction_limit: float = 1.0
    beta_n_limit: float = 3.0
    hmode_margin: float = 1.3
    radiation_fraction_limit: float = .3
    pulse_duration_s: float = 10.0
    startup_multiplier: float = 1.0
    sensitivity_fraction: float = .2

    def __post_init__(self):
        if isinstance(self.device, dict):
            object.__setattr__(self, 'device', ReferenceTokamak(**self.device))
        if not isinstance(self.device, ReferenceTokamak):
            raise ValueError('device must be a ReferenceTokamak')
        if self.radiation_convention not in ('transport', 'total'):
            raise ValueError('radiation_convention must be transport or total')
        for name in ('alpha_deposition', 'thermal_efficiency', 'heat_recovery', 'radiation_fraction_limit'):
            _number(getattr(self, name), name, 0, 1)
        _number(self.heater_efficiency, 'heater_efficiency', 0, 1, positive=True)
        for name in ('plant_load_mw', 'extra_radiation_mw', 'heating_capacity_mw'):
            _number(getattr(self, name), name)
        for name in ('greenwald_fraction_limit', 'beta_n_limit', 'pulse_duration_s'):
            _number(getattr(self, name), name, positive=True)
        _number(self.hmode_margin, 'hmode_margin', 1)
        _number(self.startup_multiplier, 'startup_multiplier', 1)
        _number(self.sensitivity_fraction, 'sensitivity_fraction', 0, .5)


@dataclass(frozen=True)
class StudyGrid:
    temperatures_kev: tuple = (5, 7.5, 10, 12.5, 15, 20, 25, 30)
    densities_m3: tuple = (5e19, 7.5e19, 1e20, 1.15e20)
    h_factors: tuple = (.8, 1.0, 1.1, 1.2)
    tritium_fractions: tuple = (.1, .25, .5)

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            values = tuple(dict.fromkeys(getattr(self, name)))
            if not values:
                raise ValueError(f'{name} cannot be empty')
            for value in values:
                if name == 'temperatures_kev':
                    _number(value, name, .2, 30)
                elif name == 'tritium_fractions':
                    _number(value, name, 0, 1)
                else:
                    _number(value, name, positive=True)
            object.__setattr__(self, name, values)
        if math.prod(len(getattr(self, name)) for name in self.__dataclass_fields__) > 5000:
            raise ValueError('study grid is limited to 5000 points')


def evaluate_point(temperature_kev, density_m3, h_factor=1.0, tritium_fraction=.5, config=None):
    config = config if config is not None else StudyConfig()
    _number(h_factor, 'h_factor', positive=True)
    _number(tritium_fraction, 'tritium_fraction', 0, 1)
    device = config.device
    try:
        raw = evaluate_case(temperature_kev, density_m3, PlantConfig(
            volume_m3=device.volume_m3, tritium_fraction=tritium_fraction,
            alpha_deposition=config.alpha_deposition))
        mass = 2 + tritium_fraction
        coefficient = ipb98y2_time(device.current_ma, device.field_t, density_m3, 1.0,
                                  device.major_radius_m, device.elongation,
                                  device.minor_radius_m / device.major_radius_m, mass, h_factor)
        # U[MJ] = A P[MW]^(1-0.69). Solve the scaling/power feedback analytically.
        scaled_loss_mw = (raw.stored_energy_j / 1e6 / coefficient)**(1 / .31)
        scaling_time = ipb98y2_time(device.current_ma, device.field_t, density_m3, scaled_loss_mw,
                                   device.major_radius_m, device.elongation,
                                   device.minor_radius_m / device.major_radius_m, mass, h_factor)
        radiation = raw.bremsstrahlung_loss_w + config.extra_radiation_mw * 1e6
        if config.radiation_convention == 'transport':
            transport, total_loss = scaled_loss_mw * 1e6, scaled_loss_mw * 1e6 + radiation
        else:
            total_loss = scaled_loss_mw * 1e6
            transport = total_loss - radiation
        required = total_loss - raw.alpha_heating_w  # Signed, never silently clamped.
        threshold = martin_threshold_mw(density_m3, device.field_t, device.surface_m2, mass) * 1e6
        greenwald = density_m3 / greenwald_density(device.current_ma, device.minor_radius_m)
        pressure = 2 * density_m3 * temperature_kev * KEV_J
        beta_n = 100 * 2 * MU0 * pressure / device.field_t**2 * device.minor_radius_m * device.field_t / device.current_ma
        rad_fraction = radiation / total_loss
        reasons = []
        if transport <= 0:
            reasons.append('radiation_exceeds_total_loss')
        if required < 0:
            reasons.append('self_heating_surplus_no_fixed_point')
        if required > config.heating_capacity_mw * 1e6:
            reasons.append('heating_capacity')
        if greenwald > config.greenwald_fraction_limit:
            reasons.append('greenwald_density')
        if beta_n > config.beta_n_limit:
            reasons.append('thermal_beta')
        if transport < config.hmode_margin * threshold:
            reasons.append('hmode_threshold')
        if rad_fraction > config.radiation_fraction_limit:
            reasons.append('radiation_fraction')
        startup = config.startup_multiplier * raw.stored_energy_j / config.heater_efficiency
        conditional_net = recovered = gross = heater = gain = None
        if required >= 0 and transport > 0:
            recovered = config.heat_recovery * (raw.fusion_power_w + required)
            gross = config.thermal_efficiency * recovered
            heater = required / config.heater_efficiency
            gain = raw.fusion_power_w / required if required > 0 else None
            conditional_net = gross - heater - config.plant_load_mw * 1e6
        eligible = not reasons
        net = conditional_net if eligible else None
        net_energy = net * config.pulse_duration_s - startup if net is not None else None
        payback = startup / net if net is not None and net > 0 else None
        row = {
            'temperature_kev': temperature_kev, 'density_m3': density_m3, 'h_factor': h_factor,
            'tritium_fraction': tritium_fraction, 'radiation_convention': config.radiation_convention,
            'volume_m3': device.volume_m3, 'surface_m2': device.surface_m2,
            'stored_energy_j': raw.stored_energy_j, 'fusion_power_w': raw.fusion_power_w,
            'alpha_heating_w': raw.alpha_heating_w, 'radiation_loss_w': radiation,
            'reactions_s': raw.reactions_s, 'reactivity_m3_s': raw.reactivity_m3_s,
            'recovered_heat_w': recovered, 'gross_electric_w': gross, 'heater_electric_w': heater,
            'plant_load_w': config.plant_load_mw * 1e6, 'plasma_gain_q': gain,
            'scaled_loss_power_w': scaled_loss_mw * 1e6, 'scaling_confinement_s': scaling_time,
            'transport_loss_w': transport,
            'transport_confinement_s': raw.stored_energy_j / transport if transport > 0 else None,
            'total_loss_w': total_loss, 'required_heating_w': required,
            'hmode_threshold_w': threshold, 'greenwald_fraction': greenwald,
            'beta_n': beta_n, 'radiation_fraction': rad_fraction,
            'startup_electric_energy_j': startup, 'conditional_net_electric_w': conditional_net,
            'net_electric_w': net, 'net_electric_energy_j': net_energy,
            'startup_payback_s': payback, 'passes_screens': eligible, 'exclusion_reasons': reasons,
            'scaling_relative_residual': scaling_time * scaled_loss_mw * 1e6 / raw.stored_energy_j - 1,
            'plasma_balance_residual_w': raw.alpha_heating_w + required - total_loss,
        }
        if any(isinstance(v, (int, float)) and not math.isfinite(v) for v in row.values()):
            raise ValueError('study result exceeds floating-point range')
        return row
    except OverflowError as exc:
        raise ValueError('study exceeds floating-point range; reduce inputs') from exc


def _corner_config(config, signs):
    changes = {}
    fields = ('thermal_efficiency', 'heat_recovery', 'heater_efficiency', 'plant_load_mw', 'startup_multiplier')
    for name, sign in zip(fields, signs):
        value = getattr(config, name) * (1 + sign * config.sensitivity_fraction)
        if name in ('thermal_efficiency', 'heat_recovery', 'heater_efficiency'):
            value = min(1.0, value)
        if name == 'startup_multiplier':
            value = max(1.0, value)
        changes[name] = value
    return replace(config, **changes)


def sensitivity_corners(point, config):
    """64 bounded parameter corners × two radiation conventions; no probabilities."""
    energies, excluded = [], Counter()
    count = 0
    for signs in product((-1, 1), repeat=6):
        changed = _corner_config(config, signs[1:])
        h = point['h_factor'] * (1 + signs[0] * config.sensitivity_fraction)
        for convention in ('transport', 'total'):
            count += 1
            row = evaluate_point(point['temperature_kev'], point['density_m3'], h,
                                 point['tritium_fraction'], replace(changed, radiation_convention=convention))
            if row['passes_screens']:
                energies.append(row['net_electric_energy_j'])
            else:
                excluded.update(row['exclusion_reasons'])
    return {'scenario_count': count, 'eligible_scenarios': len(energies),
            'positive_energy_scenarios': sum(e > 0 for e in energies),
            'worst_valid_energy_j': min(energies) if energies else None,
            'best_valid_energy_j': max(energies) if energies else None,
            'all_corners_positive': len(energies) == count and min(energies) > 0,
            'corner_exclusion_counts': dict(excluded)}


def one_at_a_time(point, config):
    """Separate numerical effect from transitions out of the screened regime."""
    output = []
    names = ('h_factor', 'thermal_efficiency', 'heat_recovery', 'heater_efficiency', 'plant_load_mw', 'startup_multiplier')
    for name in names:
        for sign in (-1, 1):
            base_value = point['h_factor'] if name == 'h_factor' else getattr(config, name)
            value = base_value * (1 + sign * config.sensitivity_fraction)
            if name in ('thermal_efficiency', 'heat_recovery', 'heater_efficiency'):
                value = min(1.0, value)
            if name == 'startup_multiplier':
                value = max(1.0, value)
            changed = config if name == 'h_factor' else replace(config, **{name: value})
            row = evaluate_point(point['temperature_kev'], point['density_m3'],
                                 value if name == 'h_factor' else point['h_factor'], point['tritium_fraction'], changed)
            output.append({'parameter': name, 'value': value, 'direction': sign,
                           'net_energy_j': row['net_electric_energy_j'],
                           'delta_energy_j': row['net_electric_energy_j'] - point['net_electric_energy_j']
                           if row['passes_screens'] else None,
                           'exclusion_reasons': row['exclusion_reasons']})
    alternate = replace(config, radiation_convention='total' if config.radiation_convention == 'transport' else 'transport')
    row = evaluate_point(point['temperature_kev'], point['density_m3'], point['h_factor'], point['tritium_fraction'], alternate)
    output.append({'parameter': 'radiation_convention', 'value': alternate.radiation_convention, 'direction': 0,
                   'net_energy_j': row['net_electric_energy_j'],
                   'delta_energy_j': row['net_electric_energy_j'] - point['net_electric_energy_j']
                   if row['passes_screens'] else None, 'exclusion_reasons': row['exclusion_reasons']})
    return output


def break_even_roots(rows, config):
    """Find bracketed pulse-energy zeros only inside connected eligible intervals."""
    groups = {}
    for row in rows:
        key = row['density_m3'], row['h_factor'], row['tritium_fraction']
        groups.setdefault(key, []).append(row)
    roots = []
    for (density, h, fraction), points in groups.items():
        points.sort(key=lambda row: row['temperature_kev'])
        for left, right in zip(points, points[1:]):
            if not left['passes_screens'] or not right['passes_screens']:
                continue
            if left['net_electric_energy_j'] * right['net_electric_energy_j'] > 0:
                continue
            a, b = left, right
            root = None
            for _ in range(60):
                middle = evaluate_point((a['temperature_kev'] + b['temperature_kev']) / 2, density, h, fraction, config)
                if not middle['passes_screens']:
                    root = None
                    break
                root = middle
                if (abs(middle['net_electric_energy_j']) <= 1e-9 * middle['startup_electric_energy_j']
                        or b['temperature_kev'] - a['temperature_kev'] < 1e-9):
                    break
                if a['net_electric_energy_j'] * middle['net_electric_energy_j'] <= 0:
                    b = middle
                else:
                    a = middle
            if root is not None and not any(abs(root['temperature_kev'] - r['temperature_kev']) < 1e-6
                                           and (density, h, fraction) == (r['density_m3'], r['h_factor'], r['tritium_fraction'])
                                           for r in roots):
                roots.append({k: root[k] for k in ('temperature_kev', 'density_m3', 'h_factor', 'tritium_fraction',
                                                   'net_electric_energy_j', 'net_electric_w')})
    return roots


def run_study(config=None, grid=None):
    config = config if config is not None else StudyConfig()
    grid = grid if grid is not None else StudyGrid()
    rows, exclusions = [], Counter()
    for t, n, h, f in product(grid.temperatures_kev, grid.densities_m3, grid.h_factors, grid.tritium_fractions):
        row = evaluate_point(t, n, h, f, config)
        exclusions.update(row['exclusion_reasons'])
        row['sensitivity'] = sensitivity_corners(row, config) if row['passes_screens'] else None
        rows.append(row)
    eligible = [row for row in rows if row['passes_screens']]
    best = max(eligible, key=lambda row: row['net_electric_energy_j']) if eligible else None
    return {'model': 'reference_tokamak_ipb98_screening', 'schema_version': 1,
            'config': asdict(config), 'grid': asdict(grid), 'rows': rows,
            'summary': {'sampled_points': len(rows), 'eligible_points': len(eligible),
                        'positive_power_points': sum(row['net_electric_w'] > 0 for row in eligible),
                        'positive_pulse_energy_points': sum(row['net_electric_energy_j'] > 0 for row in eligible),
                        'all_corner_positive_points': sum(row['sensitivity']['all_corners_positive'] for row in eligible),
                        'exclusion_counts': dict(exclusions)},
            'best_eligible_point': best, 'break_even_roots': break_even_roots(rows, config),
            'best_point_sensitivity': one_at_a_time(best, config) if best else []}
