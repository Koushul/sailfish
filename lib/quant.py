from __future__ import annotations

import subprocess
from pathlib import Path

from .env import simpleaf_bin


def comma_paths(paths: list[str] | str) -> str:
    if isinstance(paths, str):
        return paths
    return ",".join(str(p) for p in paths)


def run_simpleaf_quant(
    *,
    reads1: list[str] | str,
    reads2: list[str] | str,
    index: Path,
    chemistry: str,
    output: Path,
    threads: int,
    min_reads: int = 10,
    resolution: str = "cr-like",
    unfiltered_pl: bool = True,
    anndata_out: bool = True,
    log_path: Path | None = None,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        simpleaf_bin(),
        "quant",
        "--reads1",
        comma_paths(reads1),
        "--reads2",
        comma_paths(reads2),
        "--threads",
        str(threads),
        "--index",
        str(index),
        "--chemistry",
        chemistry,
        "--resolution",
        resolution,
        "--output",
        str(output),
    ]
    if unfiltered_pl:
        cmd.append("--unfiltered-pl")
        cmd.extend(["--min-reads", str(min_reads)])
    if anndata_out:
        cmd.append("--anndata-out")
    log_path = log_path or (output.parent / f"{output.name}.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
    if proc.returncode != 0:
        tail = log_path.read_text(errors="replace")[-4000:]
        raise RuntimeError(f"simpleaf quant failed (rc={proc.returncode})\n{tail}")


def simpleaf_set_paths() -> None:
    subprocess.run([simpleaf_bin(), "set-paths"], check=True)


def ensure_chemistry(name: str, geometry: str, expected_ori: str = "fw") -> None:
    lookup = subprocess.run(
        [simpleaf_bin(), "chemistry", "lookup", "--name", name],
        check=False,
        capture_output=True,
        text=True,
    )
    if lookup.returncode == 0 and name in (lookup.stdout + lookup.stderr):
        return
    subprocess.run(
        [
            simpleaf_bin(),
            "chemistry",
            "add",
            "--name",
            name,
            "--geometry",
            geometry,
            "--expected-ori",
            expected_ori,
        ],
        check=True,
    )
