from __future__ import annotations

from importlib.resources import files
from pathlib import Path


def package_root() -> Path:
    return Path(str(files('reseq_pipeline')))


def scripts_dir() -> Path:
    return package_root() / 'scripts'


def config_dir() -> Path:
    return package_root() / 'config'


def default_ref_yaml() -> str:
    return str(config_dir() / 'ref.yaml')
