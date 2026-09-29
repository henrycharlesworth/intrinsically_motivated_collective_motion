"""Loading the files a simulation writes (bin/fsm2d, bin/fsm2d_gpu, bin/fsm3d_gpu --out PREFIX):

  PREFIX.json    parameters and progress
  PREFIX.traj    float64 [steps+1, N, 4] (x, y, theta, v), or [steps+1, N, 6] (x, y, z, theta, phi, v) in 3D
  PREFIX.moves   int8    [steps, N]    the move each agent took
  PREFIX.scores  float32 [steps, N, M] each agent's score for every first move (-1: every future collides)
  PREFIX.pred    float64 [steps+1, n, 7] 3D predators (x, y, z, vx, vy, vz, caught so far; NaN before they appear)
  PREFIX.seen    uint8   [steps, N]    3D predators with vision: base-3 digit per predator (0 unaware, 1 sees, 2 remembers)

Moves: 2D 0/1 turn -/+dtheta, 2/3/4 straight on at v0 - dv / v0 / v0 + dv; 3D 0/1 yaw -/+, 2/3 pitch -/+, 4/5/6
straight on at v0 - dv / v0 / v0 + dv.
"""
import json
import os

import numpy as np


def dims_of(params):
    return 3 if params.get("dims", 2) == 3 else 2


def load_run(prefix, mmap=False):
    """Load a run as a dict with "params", "traj", "moves", "scores" (and "pred", "seen" when present). Runs still
    going, or interrupted, are fine: incomplete trailing records are dropped. mmap: memory-map the files instead of
    reading them (large runs; slices are read on demand). Missing .scores (e.g. not copied) give an empty array."""
    with open(prefix + ".json") as f:
        p = json.load(f)
    N, dims = p["N"], dims_of(p)
    w, nm = (6, 7) if dims == 3 else (4, 5)
    read = (lambda f, dt: np.memmap(f, dtype=dt, mode="r")) if mmap else (lambda f, dt: np.fromfile(f, dtype=dt))
    traj = read(prefix + ".traj", np.float64)
    traj = traj[: traj.size // (w * N) * w * N].reshape(-1, N, w)
    moves = read(prefix + ".moves", np.int8)
    moves = moves[: moves.size // N * N].reshape(-1, N)
    scores = read(prefix + ".scores", np.float32) if os.path.exists(prefix + ".scores") else np.zeros(0, np.float32)
    scores = scores[: scores.size // (nm * N) * nm * N].reshape(-1, N, nm)
    run = {"params": p, "traj": traj, "moves": moves, "scores": scores}
    if os.path.exists(prefix + ".pred"):
        n = p.get("predator", {}).get("n", 1)
        pr = np.fromfile(prefix + ".pred", np.float64)
        run["pred"] = pr[: pr.size // (7 * n) * 7 * n].reshape(-1, n, 7)
    if os.path.exists(prefix + ".seen"):
        s = np.fromfile(prefix + ".seen", np.uint8)
        run["seen"] = s[: s.size // N * N].reshape(-1, N)
    return run


def write_init(path, frame):
    """Write a frame ([N, 4] in 2D, [N, 6] in 3D) as a text file for --init (e.g. to start from a flocked state)."""
    np.savetxt(path, frame, fmt="%.17g")
