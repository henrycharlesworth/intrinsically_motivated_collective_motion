"""Numerical check of docs/NOTES.md: does the paper's score make the same decisions as visual empowerment?

  python scripts/empowerment_check.py runs/n500_t6 --tau 4 --frames 300,600,900 --agents 40

For a sample of agents, every first move b's futures (all 5^tau paths enumerated, pruned by collisions; visual fields
on 3072 bins) are scored four ways:
  empowerment   the target: (n_b / n_max) I(A; x_A + sigma Z), the information between the (uniformly chosen)
                surviving path and the future visual field seen through white noise sigma, by Monte Carlo
  var           the paper's score
  gaussian      the Gaussian upper bound (n_b / n_max) 1/2 log det(I + C / sigma^2) with all modes
  log n         survival counting, (n_b / n_max) log n_b (the no-noise limit of empowerment)
and for each noise level sigma the table gives the share of clear-cut decisions (best move ahead of the runner-up by
more than 2% of the target) on which each score picks the same move as empowerment. Brute force: keep tau at 4-5.
"""
import argparse
import os
import sys

import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from fsmtools import load_run  # noqa: E402
from fsmtools.reference2d import MOVES, futures, masks  # noqa: E402


def hamming2(X):
    """Squared distances between the futures' fields in the orthonormal basis: the fraction of bins that differ."""
    X = X.astype(np.float32)
    return ((X @ (1 - X).T + (1 - X) @ X.T) / X.shape[1]).astype(float)


def embed(D2):
    """Coordinates of the points (in the span of their differences) and their Gram matrix, from squared distances."""
    n = len(D2)
    J = np.eye(n) - 1 / n
    G = -0.5 * J @ D2 @ J
    w, U = np.linalg.eigh(G)
    keep = w > 1e-10 * max(w.max(), 1e-30)
    return U[:, keep] * np.sqrt(w[keep]), G


def mutual_information(Y, s, rng, reps=64):
    """I(A; Y_A + s Z) for A uniform over the rows of Y: mean_i E_z [log p(y | i) - log mean_j p(y | j)]."""
    n, r = Y.shape
    sq = (Y ** 2).sum(1)
    tot = 0.0
    for _ in range(reps):
        Yn = Y + s * rng.standard_normal((n, r))
        ll = -(sq[None, :] - 2 * Yn @ Y.T + (Yn ** 2).sum(1)[:, None]) / (2 * s * s)
        tot += np.mean(np.diag(ll) - (logsumexp(ll, axis=1) - np.log(n)))
    return tot / reps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--tau", type=int, default=4)
    ap.add_argument("--frames", default="300,600,900")
    ap.add_argument("--agents", type=int, default=40, help="agents per frame")
    ap.add_argument("--sigmas", default="0.05,0.1,0.2,0.3")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    run = load_run(a.run)
    rng = np.random.default_rng(a.seed)
    sigmas = [float(s) for s in a.sigmas.split(",")]
    names = ["var", "gaussian", "log n"]
    n_max = sum(len(MOVES) ** (k - 1) for k in range(1, a.tau + 1))
    agree = {(nm, s): [] for nm in names for s in sigmas}
    for fr in (int(f) for f in a.frames.split(",")):
        frame = run["traj"][fr]
        for i in rng.choice(len(frame), a.agents, replace=False):
            X, first = masks(futures(frame, i, run["params"], a.tau), run["params"], 3072)
            sc = {(nm, s): {} for nm in names + ["empowerment"] for s in sigmas}
            for b in range(len(MOVES)):
                sel = first == b
                n = int(sel.sum())
                if n == 0:
                    continue
                D2 = hamming2(X[sel])
                Y, G = embed(D2)
                S = X[sel].sum(0).astype(float)
                w = n / n_max
                for s in sigmas:
                    sc[("var", s)][b] = (S * (n - S)).sum() / n
                    sc[("log n", s)][b] = w * np.log(n)
                    sc[("gaussian", s)][b] = w * 0.5 * np.linalg.slogdet(np.eye(n) + G / (n * s * s))[1]
                    sc[("empowerment", s)][b] = w * (mutual_information(Y, s, rng) if n > 1 else 0.0)
            for s in sigmas:
                t = sc[("empowerment", s)]
                v = np.sort(list(t.values()))
                if len(v) < 2 or v[-1] - v[-2] < 0.02 * abs(v[-1]):
                    continue  # a near-tie of the target (within its Monte Carlo error)
                best = max(t, key=t.get)
                for nm in names:
                    agree[(nm, s)].append(max(sc[(nm, s)], key=sc[(nm, s)].get) == best)
    print(f"{a.run}, tau = {a.tau}: share of clear-cut decisions of Monte Carlo empowerment reproduced")
    print(f"  {'sigma':>6s}" + "".join(f"{nm:>11s}" for nm in names) + "   decisions")
    for s in sigmas:
        print(f"  {s:6g}" + "".join(f"{100 * np.mean(agree[(nm, s)]):10.0f}%" for nm in names) + f"   {len(agree[(names[0], s)])}")


if __name__ == "__main__":
    main()
