"""Run the IGNIS tokamak reference confinement and break-even screening study."""

import argparse
from dataclasses import asdict
from pathlib import Path
import webbrowser

from confinement import ReferenceTokamak
from feasibility import StudyConfig, StudyGrid, run_study
from feasibility_report import render_study_report, write_study_csv, write_study_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, prog='ignis.py study')
    device_defaults = asdict(ReferenceTokamak())
    config_defaults = asdict(StudyConfig())
    grid_defaults = asdict(StudyGrid())
    parser.add_argument('--temperatures-kev', type=float, nargs='+', default=grid_defaults['temperatures_kev'],
                        help='sampled kT values, 0.2–30 keV')
    parser.add_argument('--densities', type=float, nargs='+', default=grid_defaults['densities_m3'], help='total D+T density in m^-3')
    parser.add_argument('--h-factors', type=float, nargs='+', default=grid_defaults['h_factors'], help='multipliers on IPB98(y,2) confinement')
    parser.add_argument('--tritium-fractions', type=float, nargs='+', default=grid_defaults['tritium_fractions'])
    parser.add_argument('--radiation-convention', choices=('transport', 'total'), default='transport',
                        help='power used by confinement scaling: transported power or total loss')
    group = parser.add_argument_group('Reference device and screening inputs (not an IGNIS hardware design)')
    for name, value in {**device_defaults, **{k: v for k, v in config_defaults.items()
                                           if k not in ('device', 'radiation_convention')}}.items():
        group.add_argument('--' + name.replace('_', '-'), type=float, default=value,
                           help=f'{name.replace("_", " ")} (default: {value:g})')
    parser.add_argument('--output', type=Path, default=Path('study_results.html'))
    parser.add_argument('--csv', type=Path, help='complete case budgets and scenario inputs')
    parser.add_argument('--json', type=Path, help='complete study including sensitivities and roots')
    parser.add_argument('--no-browser', action='store_true')
    args = parser.parse_args(argv)
    try:
        paths = [p.resolve() for p in (args.output, args.csv, args.json) if p is not None]
        if len(paths) != len(set(paths)):
            raise ValueError('all output paths must be different')
        device = ReferenceTokamak(**{k: getattr(args, k) for k in device_defaults})
        config = StudyConfig(device=device, **{k: getattr(args, k) for k in config_defaults if k != 'device'})
        grid = StudyGrid(tuple(args.temperatures_kev), tuple(args.densities),
                         tuple(args.h_factors), tuple(args.tritium_fractions))
        study = run_study(config, grid)
        report = render_study_report(study)
    except (ValueError, OverflowError) as exc:
        parser.error(str(exc))
    try:
        for path in (args.output, args.csv, args.json):
            if path is not None:
                path.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding='utf-8')
        if args.csv:
            write_study_csv(args.csv, study)
        if args.json:
            write_study_json(args.json, study)
    except OSError as exc:
        parser.exit(1, f'Cannot write study: {exc}\n')
    summary = study['summary']
    print(f"Reference study: {summary['sampled_points']} sampled; {summary['eligible_points']} pass selected screens")
    print(f"Positive flat-top energy: {summary['positive_pulse_energy_points']}; "
          f"positive at every tested corner: {summary['all_corner_positive_points']}")
    print('Results are conditional screening outputs; device feasibility is unvalidated.')
    for label, path in (('HTML', args.output), ('CSV', args.csv), ('JSON', args.json)):
        if path is not None:
            print(f'{label}: {path.resolve()}')
    if not args.no_browser:
        try:
            opened = webbrowser.open(args.output.resolve().as_uri())
        except webbrowser.Error:
            opened = False
        if not opened:
            print('Open the HTML path above in a browser to view the study.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
