
import os
from importlib.resources import files

DEFAULT_REF_YAML = str(files("reseq_pipeline").joinpath("config/ref.yaml"))

def get_tool_bin(env_name: str, default: str) -> str:
    return os.environ.get(env_name, default)
