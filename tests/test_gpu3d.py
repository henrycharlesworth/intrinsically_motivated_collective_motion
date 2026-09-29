"""The 3D GPU simulator against the brute-force reference (python/fsmtools/reference3d.py), with and without a
predator, and its invariances."""
import numpy as np

from conftest import binary, digest, run
from fsmtools.reference3d import DEFAULTS, Grid, branch_scores, futures, init_state

GRID = (192, 96)


def score_state(tmp_path, X, tau, extra=()):
    np.savetxt(tmp_path / "state.txt", X, fmt="%.17g")
    run([binary("fsm3d_gpu", gpu=True), "--score", tmp_path / "state.txt", "--tau", tau, "--grid", f"{GRID[0]}x{GRID[1]}",
         "--reps", 0, "--dump", tmp_path / "s.f32", *extra])
    return np.fromfile(tmp_path / "s.f32", np.float32).reshape(len(X), 7)


def check(g, X, tau, agents, pred=None, p=DEFAULTS):
    for i in agents:
        ref = branch_scores(futures(X, i, p, tau, pred=pred), p, Grid(*GRID))
        assert np.array_equal(ref < 0, g[i] < 0), f"agent {i}: different dead moves"
        ok = ref > 0
        assert np.allclose(g[i][ok], ref[ok], rtol=2e-3), f"agent {i}: {g[i]} vs {ref}"


def test_scores_match_the_reference(tmp_path):
    X = init_state(80, 4.0 * 80 ** (1 / 3), seed=4)
    check(score_state(tmp_path, X, 3), X, 3, range(0, 80, 9))


def test_predator_scores_match_the_reference(tmp_path):
    X = init_state(80, 4.0 * 80 ** (1 / 3), seed=6)
    c = X[:, :3].mean(0)
    pred = np.array([c[0] - 25, c[1] + 3, c[2] - 2, 15.0, 0.0, 1.0])
    p = dict(DEFAULTS, pred_radius=3.0, pred_grow=4.0)
    g = score_state(tmp_path, X, 3, ["--pred-state", ",".join(map(str, pred)), "--pred-radius", 3, "--pred-grow", 4])
    check(g, X, 3, range(0, 80, 7), pred=pred, p=p)


PRED = ["--N", 200, "--tau", 3, "--box", 64, "--seed", 3, "--grid", "192x96", "--predator", "com", "--pred-start", 4,
        "--pred-dist", 100, "--pred-grow", 4, "--pred2", "intercept,0,70", "--pred-vision", "on"]
FILES = (".traj", ".moves", ".scores", ".pred", ".seen")


def test_splitting_over_gpus_does_not_change_the_run(tmp_path):
    run([binary("fsm3d_gpu", gpu=True), *PRED, "--steps", 12, "--out", tmp_path / "one"])
    run([binary("fsm3d_gpu", gpu=True), *PRED, "--steps", 12, "--gpus", "0,0", "--out", tmp_path / "two"])
    assert digest(str(tmp_path / "one"), FILES) == digest(str(tmp_path / "two"), FILES)


def test_resume_is_bit_identical(tmp_path):
    run([binary("fsm3d_gpu", gpu=True), *PRED, "--steps", 12, "--out", tmp_path / "full"])
    run([binary("fsm3d_gpu", gpu=True), *PRED, "--steps", 7, "--out", tmp_path / "part"])
    run([binary("fsm3d_gpu", gpu=True), *PRED, "--steps", 12, "--resume", "--out", tmp_path / "part"])
    assert digest(str(tmp_path / "full"), FILES) == digest(str(tmp_path / "part"), FILES)
