"""The 3D FSM model by brute force: the reference for the 3D GPU simulator (bin/fsm3d_gpu), and the only CPU path in
3D (tiny runs). docs/NOTES.md describes the model.

Model (the 2D model with one more dimension):
  * Agents are spheres of radius R with a position, a heading given by azimuth theta and elevation phi (a
    gravity-aligned frame: no roll) and a speed. Moves change one of them per step: yaw by +-dtheta, pitch by
    +-dtheta, speed v0 +- dv (not cumulative, as in 2D), or nothing (7 moves). A pitch past +-phi_max is
    unavailable: those paths are pruned like collisions. Elevations stay on the lattice of multiples of dtheta
    (the initial ones are rounded to it), so the future tree depends only on the elevation index, and the GPU
    code needs one template per index.
  * Futures: every move sequence of length tau, the others ballistic at v0 along their headings; a path dies at
    its first collision (swept spheres, as in 2D).
  * Visual state of a future: the set of directions covered by the other agents' discs (half-angle atan(R/d), or
    asin), seen in the future's own gravity-aligned frame (rotated by its azimuth about the vertical). It is
    pixelised on a latitude-longitude grid (nlon x nlat, pixel centres); a disc covers a pixel when the pixel
    centre is within its half-angle. On a latitude row the covered longitudes are one interval, from
    cos(dlon) = (cos a - sin(lat) sin(lat_c)) / (cos(lat) cos(lat_c)).
  * Score (the paper's, "var"): W_b = (1 / (2 pi n_b)) sum_pixels Omega_pix S_b (n_b - S_b), with S_b the number of
    futures of move b that see the pixel covered and Omega_pix the pixel's solid angle. As in 2D it is n_b times
    the mean squared L2 distance between the futures' visual fields (up to a constant), the low-SNR limit of
    visual empowerment.

  python -m fsmtools.reference3d tree --tau 7                               # future-tree growth
  python -m fsmtools.reference3d resolution --N 200 --tau 3 --agents 30     # decisions vs grid resolution
  python -m fsmtools.reference3d simulate --N 60 --tau 3 --steps 300 --out runs/r3d_n60_t3   # a small run on the CPU
"""
import argparse
import json
import math
import os
from multiprocessing import Pool

import numpy as np

MOVES7 = [(-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 0), (0, 0, 1)]  # (d yaw, d pitch, d speed)
DEFAULTS = dict(v0=10.0, dv=2.0, dtheta_deg=15.0, R=1.0, dt=1.0, phi_max_deg=75.0, half_angle="atan")


def headings(theta, phi):
    return np.stack([np.cos(phi) * np.cos(theta), np.cos(phi) * np.sin(theta), np.sin(phi)], -1)


def init_state(N, box, heading_sd_deg=15.0, seed=1, p=DEFAULTS):
    """Random positions in a cube, headings around +x: azimuth ~ N(0, sd), elevation ~ N(0, sd) rounded to the
    lattice of multiples of dtheta (and within +-phi_max). [N, 6]: x, y, z, theta, phi, v."""
    rng = np.random.default_rng(seed)
    sd, dth = math.radians(heading_sd_deg), math.radians(p["dtheta_deg"])
    emax = int(math.floor(p["phi_max_deg"] / p["dtheta_deg"] + 1e-9))
    X = np.zeros((N, 6))
    X[:, :3] = rng.uniform(0, box, (N, 3))
    X[:, 3] = sd * rng.standard_normal(N)
    X[:, 4] = dth * np.clip(np.rint(sd * rng.standard_normal(N) / dth), -emax, emax)
    X[:, 5] = p["v0"]
    return X


def _hits(P, P0, Qk, Qk1, D2):
    """Swept-sphere test of the steps P0 -> P against the others' steps Qk1 -> Qk (contact distance^2 D2): a hit at
    the end of the step, or inside it unless the two overlap at its start."""
    r1 = P[:, None] - Qk[None]
    r0 = P0[:, None] - Qk1[None]
    end = (r1 ** 2).sum(-1) < D2
    d = r1 - r0
    dd = (d ** 2).sum(-1)
    sc = np.where(dd > 0, -(r0 * d).sum(-1) / np.where(dd > 0, dd, 1), -1)
    mid = (sc > 0) & (sc < 1) & (((r0 + sc[..., None] * d) ** 2).sum(-1) < D2)
    start = (r0 ** 2).sum(-1) < D2
    return (end | (mid & ~start)).any(axis=1)


def predator_seen(state, pred, p, blind_deg=60.0, max_range=np.inf):
    """Which agents see the predator (vision): the direction to it is outside the agent's rear blind cone (blind_deg
    wide) and within range, and the sight line to its centre passes no other agent (sphere of radius R) first."""
    P = np.asarray(pred[:3], float)
    d = P - state[:, :3]; L = np.linalg.norm(d, axis=1)
    u = headings(state[:, 3], state[:, 4])
    ok = (L <= max_range) & ((u * d).sum(1) / np.maximum(L, 1e-300) >= math.cos(math.pi - math.radians(blind_deg) / 2))
    e = d / np.maximum(L, 1e-300)[:, None]
    seen = np.zeros(len(state), bool)
    for i in np.flatnonzero(ok):
        r = state[:, :3] - state[i, :3]; t = r @ e[i]
        perp2 = (r ** 2).sum(1) - t ** 2
        blk = (t > 0) & (t < L[i]) & (perp2 < p["R"] ** 2); blk[i] = False
        seen[i] = not blk.any()
    return seen


def futures(state, i, p, tau, moves=MOVES7, pred=None):
    """Surviving futures of agent i, per level k = 1..tau: (first move, azimuth index, elevation, positions [n, 3],
    others at level k [N-1, 3]). Positions are in agent i's yaw frame (rotated by -theta_i about z, origin at i),
    so a future's azimuth is its index times dtheta. pred: a predator (x, y, z, vx, vy, vz), predicted straight on;
    futures passing within R + p["pred_radius"] + p["pred_grow"] k of it at level k die like collisions (it is not
    seen)."""
    R, dt, v0, dv = p["R"], p["dt"], p["v0"], p["dv"]
    dth, phmax = math.radians(p["dtheta_deg"]), math.radians(p["phi_max_deg"])
    th0 = state[i, 3]
    c, s = math.cos(th0), math.sin(th0)
    rel = state[:, :3] - state[i, :3]
    rel = np.stack([c * rel[:, 0] + s * rel[:, 1], -s * rel[:, 0] + c * rel[:, 1], rel[:, 2]], -1)
    others = np.arange(len(state)) != i
    u = headings(state[others, 3] - th0, state[others, 4])
    Q = [rel[others] + k * v0 * dt * u for k in range(tau + 1)]
    D2 = 4 * R * R
    preds = [] if pred is None else [pred] if np.ndim(pred) == 1 else list(pred)  # one predator or a list
    Qps = []
    for pd in preds:
        pr = np.asarray(pd[:3], float) - state[i, :3]
        pv = np.asarray(pd[3:6], float)
        rot = lambda v: np.array([c * v[0] + s * v[1], -s * v[0] + c * v[1], v[2]])
        Qps.append([(rot(pr) + k * dt * rot(pv))[None] for k in range(tau + 1)])
    if preds:
        Dp = [R + p["pred_radius"] + p.get("pred_grow", 0.0) * k for k in range(tau + 1)]
    m0 = np.zeros(1, int); h = np.zeros(1, int); ph = np.full(1, state[i, 4]); P = np.zeros((1, 3))
    levels = []
    for k in range(1, tau + 1):
        nm0, nh, nph, nP, oP = [], [], [], [], []
        for b, (dy, dp, ds) in enumerate(moves):
            ok = np.abs(ph + dp * dth) <= phmax + 1e-12
            hh, pp = h[ok] + dy, ph[ok] + dp * dth
            v = v0 + ds * dv
            nm0.append(m0[ok] if k > 1 else np.full(ok.sum(), b)); nh.append(hh); nph.append(pp)
            nP.append(P[ok] + v * dt * headings(hh * dth, pp)); oP.append(P[ok])
        m0, h, ph, P = (np.concatenate(a) for a in (nm0, nh, nph, nP))
        P0 = np.concatenate(oP)
        hit = _hits(P, P0, Q[k], Q[k - 1], D2)
        for Qp in Qps:
            hit |= _hits(P, P0, Qp[k], Qp[k - 1], Dp[k] ** 2)
        m0, h, ph, P = m0[~hit], h[~hit], ph[~hit], P[~hit]
        levels.append((m0, h, ph, P, Q[k]))
    return levels


class Grid:
    """Latitude-longitude pixel grid: nlon longitudes (a multiple of 360 / dtheta, so a turn is a whole shift)
    by nlat latitudes; pixel centres at (i + 1/2) 2pi / nlon and -pi/2 + (j + 1/2) pi / nlat."""

    def __init__(self, nlon, nlat):
        self.nlon, self.nlat = nlon, nlat
        self.lat = -np.pi / 2 + (np.arange(nlat) + 0.5) * np.pi / nlat
        edges = -np.pi / 2 + np.arange(nlat + 1) * np.pi / nlat
        self.omega = (2 * np.pi / nlon) * np.diff(np.sin(edges))  # solid angle of a pixel in each row


def paint(P, Q, h, p, g, chunk=None):
    """Visual masks [n, nlat, nlon] (bool) of futures at positions P [n, 3] with azimuth indices h, others at Q
    [m, 3], in each future's own yaw frame. Yields (slice, masks) in chunks of futures."""
    dth = math.radians(p["dtheta_deg"])
    per = g.nlon * dth / (2 * np.pi)
    assert abs(per - round(per)) < 1e-9, "nlon must be a multiple of 360 / dtheta"
    per = int(round(per))
    chunk = chunk or max(1, int(3e7 // (g.nlat * 3 * g.nlon)))  # bounds the difference arrays' memory
    for a in range(0, len(P), chunk):
        sl = slice(a, min(a + chunk, len(P)))
        D = Q[None] - P[sl, None]  # [n, m, 3]
        d = np.linalg.norm(D, axis=-1)
        x = p["R"] / d
        alpha = np.arctan(x) if p["half_angle"] == "atan" else np.arcsin(np.minimum(1, x))
        latc = np.arcsin(np.clip(D[..., 2] / d, -1, 1))
        lonc = np.arctan2(D[..., 1], D[..., 0])
        n, m = d.shape
        # rows whose centre latitude is within alpha of the disc centre's
        j0 = np.ceil((latc - alpha + np.pi / 2) * g.nlat / np.pi - 0.5).astype(int).clip(0, g.nlat - 1)
        j1 = np.floor((latc + alpha + np.pi / 2) * g.nlat / np.pi - 0.5).astype(int).clip(0, g.nlat - 1)
        cnt = np.maximum(j1 - j0 + 1, 0)
        f_idx = np.repeat(np.arange(n)[:, None], m, 1).ravel()
        c_idx = np.repeat(np.arange(m)[None], n, 0).ravel()
        rep = cnt.ravel()
        F, C = np.repeat(f_idx, rep), np.repeat(c_idx, rep)
        J = np.repeat(j0.ravel(), rep) + (np.arange(rep.sum()) - np.repeat(np.cumsum(rep) - rep, rep))
        lat = g.lat[J]
        lc, al = latc[F, C], alpha[F, C]
        num = np.cos(al) - np.sin(lat) * np.sin(lc)
        den = np.cos(lat) * np.cos(lc)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(den > 0, num / den, np.where(num <= 0, -2.0, 2.0))
        full = r <= -1
        keep = r <= 1
        dl = np.arccos(np.clip(r, -1, 1))
        lo = lonc[F, C]
        # first/last longitude pixel whose centre is inside [lo - dl, lo + dl], unwrapped over [-nlon, 2 nlon)
        k0 = np.ceil((lo - dl) * g.nlon / (2 * np.pi) - 0.5).astype(int)
        k1 = np.floor((lo + dl) * g.nlon / (2 * np.pi) - 0.5).astype(int)
        k0 = np.where(full, 0, k0); k1 = np.where(full, g.nlon - 1, k1)
        ok = keep & (k1 >= k0)
        diff = np.zeros((n, g.nlat, 3 * g.nlon + 1), np.int32)
        np.add.at(diff, (F[ok], J[ok], k0[ok] + g.nlon), 1)
        np.add.at(diff, (F[ok], J[ok], k1[ok] + g.nlon + 1), -1)
        cov = np.cumsum(diff[..., :-1], axis=-1)
        mask = (cov[..., :g.nlon] + cov[..., g.nlon:2 * g.nlon] + cov[..., 2 * g.nlon:]) > 0
        # into each future's own frame: its longitude q is the root frame's q + h * per
        sh = (h[sl] * per) % g.nlon
        idx = (np.arange(g.nlon)[None] + sh[:, None]) % g.nlon
        yield sl, np.take_along_axis(mask, idx[:, None, :], axis=2)


def branch_scores(levels, p, g, nm=len(MOVES7), keep_masks=False):
    """var scores W_b per first move (-1 when no future survives). keep_masks: also return the masks and first
    moves of all futures."""
    S = np.zeros((nm, g.nlat, g.nlon)); n = np.zeros(nm)
    kept, firsts = [], []
    for m0, h, ph, P, Q in levels:
        for sl, M in paint(P, Q, h, p, g):
            for b in range(nm):
                sel = m0[sl] == b
                if sel.any():
                    S[b] += M[sel].sum(0); n[b] += sel.sum()
            if keep_masks:
                kept.append(M); firsts.append(m0[sl])
    W = np.array([(g.omega[:, None] * S[b] * (n[b] - S[b])).sum() / (2 * np.pi * n[b]) if n[b] else -1.0
                  for b in range(nm)])
    return (W, np.concatenate(kept), np.concatenate(firsts)) if keep_masks else W


def tree_growth(tau_max, moves, merge=0.0, phi0=0.0, p=DEFAULTS):
    """Futures per level without collisions: paths, distinct states (heading, speed, position: the DAG) and
    distinct positions (optionally merged on a grid of size merge), summed over levels 1..tau."""
    dth, phmax = math.radians(p["dtheta_deg"]), math.radians(p["phi_max_deg"])
    lvl = {(0, 0, (0.0, 0.0, 0.0))}  # (azimuth index, elevation index, position)
    rows, paths, dag, pos = [], 0, 0, 0
    npaths = 1
    for k in range(1, tau_max + 1):
        new = set()
        for (a, e, q) in lvl:
            for (dy, dp, ds) in moves:
                if abs(phi0 + (e + dp) * dth) > phmax + 1e-12:
                    continue
                v = p["v0"] + ds * p["dv"]
                u = headings(np.float64((a + dy) * dth), np.float64(phi0 + (e + dp) * dth))
                new.add((a + dy, e + dp, tuple(round(q[i] + v * p["dt"] * u[i], 6) for i in range(3))))
        npaths *= len(moves)
        paths += npaths
        dag += len(new)
        qs = {s[2] for s in new}
        if merge > 0:
            qs = {tuple(round(c / merge) for c in q) for q in qs}
        pos += len(qs)
        rows.append((k, paths, dag, pos))
        lvl = new
    return rows


def check_resolution(N, tau, n_agents, grids, spacing, seed, run=None, frames=(), p=DEFAULTS):
    """Decisions of the var score at several grids, against the finest: share of agents choosing the same move,
    the mean regret (score lost, relative to the best, when following the coarse grid's choice), and the median
    relative error of the scores. States: a random initial state, or frames of a run (run prefix, frames)."""
    rng = np.random.default_rng(seed)
    if run:
        from .runs import load_run
        r = load_run(run)
        p = {k: r["params"][k] for k in DEFAULTS}
        states = [r["traj"][f] for f in frames]
    else:
        states = [init_state(N, spacing * N ** (1 / 3), seed=seed)]
    W = {gs: [] for gs in grids}
    for X in states:
        for i in rng.choice(len(X), min(n_agents, len(X)), replace=False):
            lv = futures(X, i, p, tau)
            for gs in grids:
                W[gs].append(branch_scores(lv, p, Grid(*gs)))
    ref = np.array(W[grids[-1]])
    srt = np.sort(ref, 1)
    margin = (srt[:, -1] - srt[:, -2]) / np.maximum(srt[:, -1], 1e-30)
    src = f"{run} frames {list(frames)}" if run else f"random state N={N} spacing={spacing}"
    print(f"{src}, tau={tau}, {len(ref)} decisions vs the {grids[-1][0]}x{grids[-1][1]} grid "
          f"(best-vs-second margin: median {np.median(margin):.2%}, {np.mean(margin < 0.01):.0%} under 1%)")
    for gs in grids[:-1]:
        w = np.array(W[gs])
        pick = w.argmax(1)
        same = pick == ref.argmax(1)
        best = ref.max(1)
        regret = (best - ref[np.arange(len(ref)), pick]) / np.maximum(best, 1e-30)
        ok = ref > 0
        rel = np.abs(w[ok] - ref[ok]) / ref[ok]
        print(f"  {gs[0]:5d} x {gs[1]:4d}: same move {100 * same.mean():5.1f}%   mean regret {regret.mean():.2e}"
              f"   max {regret.max():.2e}   median score error {np.median(rel):.2e}")


def _agent_scores(args):
    X, i, p, tau, gs = args
    return branch_scores(futures(X, i, p, tau), p, Grid(*gs))


def step(X, scores, p, rng, moves=MOVES7):
    """Synchronous update: every agent takes its best move (uniform among ties, as in 2D)."""
    dth = math.radians(p["dtheta_deg"])
    Y = X.copy()
    mv = np.zeros(len(X), np.int8)
    for i, w in enumerate(scores):
        best = np.flatnonzero(w >= w.max() - 1e-9 * max(1.0, abs(w.max())))
        m = best[rng.integers(len(best))] if w.max() >= 0 else 5  # everything collides: stay
        dy, dp, ds = moves[m]
        Y[i, 3] += dy * dth
        Y[i, 4] += dp * dth
        Y[i, 5] = p["v0"] + ds * p["dv"]
        mv[i] = m
    Y[:, :3] += (Y[:, 5] * p["dt"])[:, None] * headings(Y[:, 3], Y[:, 4])
    return Y, mv


def simulate(N, tau, steps, out, grid, spacing, seed, workers, p=DEFAULTS):
    """A small run with the reference (slow; no DAG, no pruning of the tree before painting). Writes the same files
    as the simulators (fsmtools.load_run). Note: its initial state uses numpy's generator, not the simulators'
    random streams, so it does not reproduce a GPU run of the same seed."""
    X = init_state(N, spacing * N ** (1 / 3), seed=seed)
    rng = np.random.default_rng(seed + 1)
    meta = dict(p, dims=3, N=N, tau=tau, steps=0, grid=list(grid), box=spacing * N ** (1 / 3), seed=seed,
                backend="python-reference")
    ft, fm, fs = open(out + ".traj", "wb"), open(out + ".moves", "wb"), open(out + ".scores", "wb")
    ft.write(X.astype(np.float64).tobytes())
    with Pool(workers) as pool:
        for t in range(steps):
            W = pool.map(_agent_scores, [(X, i, p, tau, grid) for i in range(N)])
            X, mv = step(X, W, p, rng)
            for f, arr in ((ft, X.astype(np.float64)), (fm, mv), (fs, np.array(W, np.float32))):
                f.write(arr.tobytes()); f.flush()
            if t % 10 == 0 or t == steps - 1:
                u = headings(X[:, 3], X[:, 4]).mean(0)
                print(f"step {t}  order {np.linalg.norm(u):.3f}", flush=True)
                meta["steps"] = t + 1
                json.dump(meta, open(out + ".json", "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("check", choices=["tree", "resolution", "simulate"])
    ap.add_argument("--tau", type=int, default=3)
    ap.add_argument("--merge", type=float, default=0.0)
    ap.add_argument("--N", type=int, default=200)
    ap.add_argument("--agents", type=int, default=30)
    ap.add_argument("--spacing", type=float, default=11.0, help="initial mean spacing (box = spacing N^(1/3))")
    ap.add_argument("--grids", default="192x96,384x192,768x384,1536x768")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--grid", default="384x192")
    ap.add_argument("--out", default="runs/r3d")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--run", help="resolution: take the states from frames of this run")
    ap.add_argument("--frames", default="")
    a = ap.parse_args()
    if a.check == "tree":
        for k, paths, dag, pos in tree_growth(a.tau, MOVES7, a.merge):
            print(f"tau={k}: paths {paths:>10d}  dag {dag:>9d}  positions {pos:>9d}")
    elif a.check == "simulate":
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        simulate(a.N, a.tau, a.steps, a.out, tuple(int(v) for v in a.grid.split("x")), a.spacing, a.seed, a.workers)
    else:
        grids = [tuple(int(v) for v in s.split("x")) for s in a.grids.split(",")]
        check_resolution(a.N, a.tau, a.agents, grids, a.spacing, a.seed, a.run,
                         [int(f) for f in a.frames.split(",")] if a.frames else ())
