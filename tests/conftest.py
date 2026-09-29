"""Shared helpers: paths to the binaries (tests skip when one is not built, or when there is no GPU for the CUDA
ones) and a runner for simulator commands."""
import hashlib
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))


def binary(name, gpu=False):
    path = os.path.join(ROOT, "bin", name)
    if not os.path.exists(path):
        pytest.skip(f"bin/{name} is not built (make{' gpu' if gpu else ''})")
    if gpu and (shutil.which("nvidia-smi") is None or subprocess.call(["nvidia-smi", "-L"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)):
        pytest.skip("no NVIDIA GPU")
    return path


def run(cmd, cwd=None):
    """Run a simulator; returns its stderr (the log). Fails the test on a non-zero exit."""
    r = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, f"{' '.join(map(str, cmd))} failed:\n{r.stderr[-2000:]}"
    return r.stderr


def digest(prefix, exts=(".traj", ".moves", ".scores")):
    """md5 of a run's files, to check that two runs are bit-identical."""
    return {e: hashlib.md5(open(prefix + e, "rb").read()).hexdigest() for e in exts if os.path.exists(prefix + e)}
