from __future__ import annotations

import argparse
from pathlib import Path

from . import __version__
from .pipeline import run_pipeline
from .resources import default_ref_yaml
from .toolchain import validate_toolchain


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('-i', '--input', required=True, help='Sample manifest file')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument(
        '-r', '--ref-yaml', dest='ref_yaml', default=default_ref_yaml(),
        help='Reference YAML file path. Default: packaged ref.yaml'
    )
    parser.add_argument('-g', '--genome', default='hg38', help='Genome key in ref yaml, e.g. hg38/mm10/rn6')
    parser.add_argument('--trim', type=int, default=5, help='Trim length, default=5')
    parser.add_argument('--qvalue', type=int, default=20, help='Base quality threshold, default=20')
    parser.add_argument('--thread', type=int, default=4, help='Threads, default=4')
    parser.add_argument('--concurrent', type=int, default=10, help='Max concurrent jobs')
    parser.add_argument('--refresh', type=int, default=30, help='Refresh time in seconds')
    parser.add_argument('--job-type', dest='job_type', choices=['sge', 'local'], default='local')
    parser.add_argument('--work-dir', dest='work_dir', default='.', help='Work directory')
    parser.add_argument('-cn', '--contract', help='Contract/project number for report')
    parser.add_argument('-pn', '--name', help='Project name for report')
    parser.add_argument('--diff-type', dest='diff_type', default='padj', help='padj or pval')
    parser.add_argument('--diff-type-num', dest='diff_type_num', default='0.1', help='Threshold value')
    parser.add_argument('--skip-tool-check', action='store_true', help='Skip external tool validation before run')


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='RE-seq pipeline')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    subparsers = parser.add_subparsers(dest='command')

    run_parser = subparsers.add_parser('run', help='Run the pipeline')
    add_run_arguments(run_parser)

    check_parser = subparsers.add_parser('check', help='Check required tools and packaged resources')
    check_parser.add_argument('-r', '--ref-yaml', dest='ref_yaml', default=default_ref_yaml(), help='Reference YAML path')

    return parser


def main() -> None:
    parser = build_parser()
    args, extras = parser.parse_known_args()

    # Backward compatibility: allow `reseq -i ... -o ...` without explicit subcommand.
    if args.command is None and extras == [] and hasattr(args, 'input'):
        # defensive; argparse won't reach here normally
        pass
    if args.command is None:
        legacy_parser = argparse.ArgumentParser(description='RE-seq pipeline')
        legacy_parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
        add_run_arguments(legacy_parser)
        args = legacy_parser.parse_args()
        command = 'run'
    else:
        command = args.command

    if command == 'check':
        statuses = validate_toolchain(strict=False)
        print('Tool check:')
        for status in statuses:
            flag = 'OK' if status.found else ('WARN' if not status.required else 'MISS')
            print(f'[{flag}] {status.name} -> {status.detail}')
        ref_yaml = Path(args.ref_yaml)
        print(f"[{'OK' if ref_yaml.exists() else 'MISS'}] REF_YAML -> {ref_yaml}")
        return

    if not args.skip_tool_check:
        validate_toolchain(strict=True)

    run_pipeline(
        input=args.input,
        out_dir=args.outdir,
        ref_yaml=args.ref_yaml,
        genome=args.genome,
        trim=args.trim,
        qvalue=args.qvalue,
        thread=args.thread,
        concurrent=args.concurrent,
        refresh=args.refresh,
        job_type=args.job_type,
        work_dir=args.work_dir,
        contract=args.contract,
        name=args.name,
        diff_type=args.diff_type,
        diff_type_num=args.diff_type_num,
    )
