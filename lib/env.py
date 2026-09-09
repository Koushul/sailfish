from __future__ import annotations

import os
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
AF_TUTORIAL = PIPELINE_ROOT.parent
CONDA_BIN = AF_TUTORIAL / "conda_env" / "bin"
CONDA_LIB = AF_TUTORIAL / "conda_env" / "lib"
CARGO_BIN = AF_TUTORIAL / "cargo_tools" / "bin"


def setup_env(alevin_fry_home: Path) -> None:
    alevin_fry_home.mkdir(parents=True, exist_ok=True)
    os.environ["PYTHONNOUSERSITE"] = "1"
    os.environ["ALEVIN_FRY_HOME"] = str(alevin_fry_home)
    os.environ["PATH"] = f"{CONDA_BIN}:{CARGO_BIN}:{os.environ.get('PATH', '')}"
    conda_lib = str(CONDA_LIB)
    ld = os.environ.get("LD_LIBRARY_PATH", "")
    parts = [p for p in ld.split(":") if p and p != conda_lib]
    os.environ["LD_LIBRARY_PATH"] = ":".join([conda_lib] + parts)
    try:
        os.chdir(os.getcwd())
        import resource

        resource.setrlimit(resource.RLIMIT_NOFILE, (2048, 2048))
    except (ValueError, OSError):
        pass


def python_bin() -> str:
    return str(CONDA_BIN / "python")


def simpleaf_bin() -> str:
    return str(CONDA_BIN / "simpleaf")
