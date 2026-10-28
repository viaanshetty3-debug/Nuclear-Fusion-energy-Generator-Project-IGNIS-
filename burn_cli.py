"""Run IGNIS's time-dependent D–D/D–T model with fixed auxiliary heating."""

import argparse
from dataclasses import asdict
from pathlib import Path
import webbrowser

from burn import BurnConfig, convergence_check
from burn_report import render_burn_report, write_burn_csv, write_burn_json
from fusion import KEV_J

K_KEV_K = 1.380649e-23 / KEV_J


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, prog='ignis.py burn')
    parser.add_argument('--temperature-k', type=float, help='initial ion/electron temperature (K); default kT=10 keV')
    defaults = asdict(BurnConfig())
    help_text = {
        'density_m3': 'initial total fuel ion density (m^-3)', 'tritium_fraction': 'initial T/(D+T) number fraction',
        'volume_m3': 'fixed plasma volume (m^3)', 'confinement_s': 'transport energy timescale excluding radiation (s)',
        'charged_deposition': 'fraction of charged-product birth energy deposited in plasma',
        'heating_mw': 'constant deposited auxiliary heating (MW)', 'heater_efficiency': 'electricity-to-heating efficiency',
        'thermal_efficiency': 'heat-to-electricity efficiency', 'heat_recovery': 'outgoing heat recovery fraction',
        'plant_load_mw': 'other electrical plant demand (MW)', 'other_loss_mw': 'additional plasma heat loss (MW)',
        'duration_s': 'requested burn duration (s)', 'max_step_s': 'maximum solver step (s)',
        'rtol': 'relative solver tolerance; convergence run tightens it up to 16 times', 'samples': 'output time samples',
    }
    for name, default in defaults.items():
        if name == 'temperature_kev':
            continue
        parser.add_argument('--' + name.replace('_', '-'), type=int if name == 'samples' else float,
                            default=default, help=f'{help_text[name]} (default: {default:g})')
    parser.add_argument('--output', type=Path, default=Path('burn_results.html'))
    parser.add_argument('--csv', type=Path, help='export time series and inputs')
    parser.add_argument('--json', type=Path, help='export time series, inputs, diagnostics and convergence')
    parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args(argv)
    try:
        paths = [p.resolve() for p in (args.output, args.csv, args.json) if p is not None]
        if len(paths) != len(set(paths)):
            raise ValueError('all output paths must be different')
        inputs = {name: getattr(args, name) for name in defaults if name != 'temperature_kev'}
        if args.temperature_k is not None:
            inputs['temperature_kev'] = args.temperature_k * K_KEV_K
        config = BurnConfig(**inputs)
        result, convergence = convergence_check(config)
        report = render_burn_report(result, convergence)
    except (ValueError, OverflowError) as exc:
        parser.error(str(exc))
    try:
        for path in (args.output, args.csv, args.json):
            if path is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding='utf-8')
        if args.csv:
            write_burn_csv(args.csv, result)
        if args.json:
            write_burn_json(args.json, result, convergence)
    except OSError as exc:
        parser.exit(1, f'Cannot write results: {exc}\n')
    last = result.rows[-1]
    print(f"Burn status: {result.status}; simulated {last['time_s']:.8g} / {config.duration_s:g} s")
    print(f"Final kT: {last['temperature_kev']:.6g} keV; remaining T: {last['tritium_density_m3']:.6g} m^-3")
    print(f"Fusion energy: {last['fusion_energy_j'] / 1e6:.6g} MJ; net electricity including startup: "
          f"{last['net_electric_energy_j'] / 1e6:.6g} MJ")
    print(f"Convergence: {'PASS' if convergence['passed'] else 'FAIL'}")
    for label, path in (('HTML', args.output), ('CSV', args.csv), ('JSON', args.json)):
        if path is not None:
            print(f'{label}: {path.resolve()}')
    if not args.no_browser:
        try:
            opened = webbrowser.open(args.output.resolve().as_uri())
        except webbrowser.Error:
            opened = False
        if not opened:
            print('Open the HTML path above in a browser to view the report.')
    return 0 if convergence['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
