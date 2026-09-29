"""Export runs for the standalone browser viewer (viewer/): one script file per run, data/<id>.js, plus data/runs.js
listing them (see viewer/js/data.js for the format).

A run is stored as its initial state and every agent's move at every step. The viewer re-integrates the
trajectory, exactly as the simulators do, so nothing is lost; moves carry ~2.3 bits each and deflate to ~2.6 bits,
~40x smaller than the float32 positions. Runs with motor noise (which the moves do not determine) ship float32
positions instead. The per-step scores are left out (they do not compress): the 2D viewer recomputes any bird's
decision itself.
"""
import base64
import json
import os
import subprocess
import tempfile
import zlib

import numpy as np

from .runs import dims_of, load_run

VIEWER = os.path.join(os.path.dirname(__file__), "..", "..", "viewer")
CHUNK = 1 << 23  # base64 characters per chunk (the viewer decodes chunk by chunk)
MAX_FRAMES_BYTES = 1.2e9  # the browser holds a 2D run's frames as float32 [T, N, 3]; above this keep 1 in `stride`


def _array(arr, dtype):
    raw = np.ascontiguousarray(arr, dtype=dtype).tobytes()
    c = zlib.compressobj(9, zlib.DEFLATED, -15)  # raw deflate (DecompressionStream "deflate-raw")
    z = c.compress(raw) + c.flush()
    b64 = base64.b64encode(z).decode("ascii")
    kind = {np.uint8: "u8", np.float32: "f32", np.float64: "f64"}[dtype]
    return {"type": kind, "deflate": True, "chunks": [b64[i:i + CHUNK] for i in range(0, len(b64), CHUNK)]}


def replay2d(params, init, moves):
    """The 2D trajectory [T, N, 3] (x, y, theta) from the initial frame [N, 4] and the moves [T-1, N], as the
    simulators integrate it (no noise). The viewer's replay (fsm-core.js replay2d) does the same."""
    dth = np.radians(params["dtheta_deg"])
    dh = np.array([-1, 1, 0, 0, 0], float)[moves]
    ds = np.array([0, 0, -1, 0, 1], float)[moves]
    th = init[:, 2] + np.cumsum(dh * dth, axis=0)
    v = params["v0"] + ds * params["dv"]
    x = init[:, 0] + np.cumsum(v * params["dt"] * np.cos(th), axis=0)
    y = init[:, 1] + np.cumsum(v * params["dt"] * np.sin(th), axis=0)
    first = np.stack([init[:, 0], init[:, 1], init[:, 2]], -1)[None]
    return np.concatenate([first, np.stack([x, y, th], -1)])


def default_label(params, steps):
    head = ("3D · " if dims_of(params) == 3 else "") + f"N = {params['N']} · τ = {params['tau']}"
    if params.get("predator"):
        pr = params["predator"]
        head += f" · predator{'s' if pr.get('n', 1) > 1 else ''}"
    return f"{head} · {steps} steps"


def _analysis2d(params, init, moves, T, N, stride):
    """The 2D viewer's per-frame analysis (camera tracks, order, fragments) computed ahead with its own code
    (viewer/tools/analyse.js, needs node), so that opening a large run is fast. None when node is not available."""
    try:
        with tempfile.TemporaryDirectory() as d:
            np.asarray(init, np.float64).tofile(os.path.join(d, "init.f64"))
            np.asarray(moves, np.uint8).tofile(os.path.join(d, "moves.u8"))
            with open(os.path.join(d, "meta.json"), "w") as f:
                json.dump({"params": params, "T": T, "N": N, "stride": stride}, f)
            subprocess.run(["node", os.path.join(VIEWER, "tools", "analyse.js"), d], check=True, capture_output=True)
            return np.fromfile(os.path.join(d, "analysis.f32"), np.float32)
    except (OSError, subprocess.CalledProcessError):
        return None


def export_run(prefix, out_dir, run_id=None, label=None, stride=None, scores=False):
    """Write out_dir/<run_id>.js for the run at `prefix`; returns its entry for runs.js. stride (2D): keep one frame in
    this many in the browser (default: automatic for runs too large to hold every frame). scores (3D): include the
    per-step scores, which the 3D bird panel then shows (large)."""
    run = load_run(prefix, mmap=True)
    p, traj, moves = run["params"], run["traj"], run["moves"]
    dims = dims_of(p)
    T = min(len(traj), len(moves) + 1)
    N = traj.shape[1]
    run_id = run_id or os.path.basename(prefix)
    meta = {"params": p, "T": int(T), "N": int(N), "dims": dims}
    arrays = {"moves": _array(moves[: T - 1], np.uint8)}
    if dims == 2 and p.get("noise", 0) > 0:  # the moves do not determine the trajectory
        arrays["traj"] = _array(traj[:T, :, :3], np.float32)
    elif dims == 2:
        arrays["init"] = _array(traj[0], np.float64)
        if stride is None:
            stride = max(1, int(np.ceil(T * N * 12 / MAX_FRAMES_BYTES)))
        if stride > 1:
            meta["stride"] = int(stride)
        an = _analysis2d(p, traj[0], moves[: T - 1], T, N, stride)
        if an is not None:
            arrays["analysis"] = _array(an, np.float32)
    else:
        arrays["init"] = _array(traj[0], np.float64)
        if "pred" in run:
            pr = run["pred"][:T]
            arrays["pred"] = _array(pr, np.float32)
            meta["pred"], meta["pred_n"] = True, int(pr.shape[1])
        if "seen" in run:
            arrays["seen"] = _array(run["seen"][: T - 1], np.uint8)
            meta["seen"] = True
        if scores and len(run["scores"]) >= T - 1:
            arrays["scores"] = _array(run["scores"][: T - 1], np.float32)
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, run_id + ".js"), "w") as f:
        f.write("FSMData.register(" + json.dumps(run_id) + ", " + json.dumps({"meta": meta, "arrays": arrays}) + ");\n")
    return {"id": run_id, "dims": dims, "label": label or default_label(p, T - 1)}


def write_index(out_dir, entries):
    """data/runs.js: the bundled runs, in menu order (the first of each dimension opens by default)."""
    with open(os.path.join(out_dir, "runs.js"), "w") as f:
        f.write("FSMData.list(" + json.dumps(entries, ensure_ascii=False, indent=1) + ");\n")
