from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .resources import scripts_dir


@dataclass(frozen=True)
class ToolStatus:
    name: str
    value: str
    found: bool
    required: bool
    detail: str


def env_or_default(env_name: str, default: str) -> str:
    return os.environ.get(env_name, default)


def resolve_command(cmd: str) -> tuple[bool, str]:
    if os.path.sep in cmd or cmd.startswith('.'):
        return Path(cmd).exists(), cmd
    resolved = shutil.which(cmd)
    return (resolved is not None), (resolved or cmd)


def packaged_script(name: str) -> str:
    path = scripts_dir() / name
    return str(path)


def required_tools_for_run() -> list[tuple[str, str, bool]]:
    return [
        ('PYTHON_BIN', env_or_default('PYTHON_BIN', 'python'), True),
        ('FASTP_BIN', env_or_default('FASTP_BIN', 'fastp'), True),
        ('FASTQC_BIN', env_or_default('FASTQC_BIN', 'fastqc'), True),
        ('HISAT2_BIN', env_or_default('HISAT2_BIN', 'hisat2'), True),
        ('SAMTOOLS_BIN', env_or_default('SAMTOOLS_BIN', 'samtools'), True),
        ('HTSEQ', env_or_default('HTSEQ', 'htseq-count'), True),
        ('GATK_BIN', env_or_default('GATK_BIN', 'gatk'), True),
        ('PERL_BIN', env_or_default('PERL_BIN', 'perl'), True),
        ('RSCRIPT', env_or_default('RSCRIPT', 'Rscript'), True),
        ('BGZIP', env_or_default('BGZIP', 'bgzip'), True),
        ('TABIX', env_or_default('TABIX', 'tabix'), True),
        ('ANNOVAR_BIN', env_or_default('ANNOVAR_BIN', 'annovar'), True),
        ('JAVA_BIN', env_or_default('JAVA_BIN', 'java'), False),
    ]


def collect_tool_statuses(items: Iterable[tuple[str, str, bool]] | None = None) -> list[ToolStatus]:
    items = list(items) if items is not None else required_tools_for_run()
    statuses: list[ToolStatus] = []
    for name, value, required in items:
        found, resolved = resolve_command(value)
        detail = resolved if found else f'not found: {value}'
        statuses.append(ToolStatus(name=name, value=value, found=found, required=required, detail=detail))
    return statuses


def validate_toolchain(strict: bool = True) -> list[ToolStatus]:
    statuses = collect_tool_statuses()
    missing = [s for s in statuses if s.required and not s.found]
    if strict and missing:
        summary = ', '.join(f'{m.name}={m.value}' for m in missing)
        raise RuntimeError(
            'Missing required external tools. Set them in PATH or export env vars first: '
            + summary
        )
    return statuses
