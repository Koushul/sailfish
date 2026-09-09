from __future__ import annotations

import os
import shutil
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
_SIMPLEAF = "simpleaf"


def _prepend_path(bin_dir: Path) -> None:
    os.environ["PATH"] = f"{bin_dir}:{os.environ.get('PATH', '')}"


def _prepend_ld(lib_dir: Path) -> None:
    if not lib_dir.is_dir():
        return
    lib = str(lib_dir)
    rest = [p for p in os.environ.get("LD_LIBRARY_PATH", "").split(":") if p and p != lib]
    os.environ["LD_LIBRARY_PATH"] = ":".join([lib] + rest)


def setup_env(alevin_fry_home: Path, tools: dict | None = None) -> None:
    alevin_fry_home.mkdir(parents=True, exist_ok=True)
    os.environ["ALEVIN_FRY_HOME"] = str(alevin_fry_home)
    os.environ.setdefault("PYTHONNOUSERSITE", "1")

    tools = tools or {}
    conda = tools.get("conda") or os.environ.get("SAILFISH_CONDA")
    cargo = tools.get("cargo") or os.environ.get("SAILFISH_CARGO")
    if not conda:
        sibling = PIPELINE_ROOT.parent / "conda_env"
        if sibling.is_dir():
            conda = str(sibling)
    if not cargo:
        sibling_cargo = PIPELINE_ROOT.parent / "cargo_tools"
        if sibling_cargo.is_dir():
            cargo = str(sibling_cargo)

    if conda:
        c = Path(conda)
        _prepend_path(c / "bin")
        _prepend_ld(c / "lib")
    if cargo:
        _prepend_path(Path(cargo) / "bin")

    try:
        import resource

        resource.setrlimit(resource.RLIMIT_NOFILE, (2048, 2048))
    except (ValueError, OSError):
        pass

    exe = shutil.which("simpleaf")
    if not exe:
        raise FileNotFoundError(
            "simpleaf not found. Put it on PATH, or set SAILFISH_CONDA / config tools.conda "
            "to a conda env that contains simpleaf, piscem, and alevin-fry."
        )


def simpleaf_bin() -> str:
    exe = shutil.which("simpleaf")
    if not exe:
        raise FileNotFoundError("simpleaf not found on PATH")
    return exe
