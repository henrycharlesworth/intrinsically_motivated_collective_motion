"""Tools for FSM collective-motion runs: loading run files, flock measures, brute-force reference models (2D and
3D) and export for the browser viewer. Use with `PYTHONPATH=python` or `pip install -e .`."""
from .runs import load_run, dims_of

__all__ = ["load_run", "dims_of"]
