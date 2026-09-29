"""The 2D GPU simulator against the exact CPU simulator and the brute-force reference, and its invariances."""
import json

import numpy as np

from conftest import binary, digest, run
from fsmtools import load_run
from fsmtools.reference2d import agent_scores

MODEL = ["--N", 300, "--tau", 5, "--box", 194, "--seed", 5]


def gpu_scores(dump):
    meta = json.load(open(str(dump) + ".json"))
    return np.fromfile(dump, np.float32).reshape(len(meta["frames"]), meta["N"], 5)


def test_gpu_makes_the_cpu_decisions(tmp_path):
    cpu = tmp_path / "cpu"
    run([binary("fsm2d"), *MODEL, "--steps", 40, "--out", cpu])
    frames = [10, 25, 39]
    run([binary("fsm2d_gpu", gpu=True), "--score", cpu, "--frames", ",".join(map(str, frames)), "--dump", tmp_path / "g.f32"])
    g = gpu_scores(tmp_path / "g.f32")
    c = load_run(str(cpu))["scores"][frames]  # scores used for the decision at frame f
    ok = c.max(-1) > 0
    top = np.sort(c, -1)
    clear = ok & ((top[..., -1] - top[..., -2]) > 1e-3 * top[..., -1])  # not a near-tie
    assert (g.argmax(-1) == c.argmax(-1))[clear].mean() >= 0.995
    valid = c > 0
    assert np.median(np.abs(g[valid] - c[valid]) / c[valid]) < 1e-3


def test_gpu_matches_the_brute_force_reference(tmp_path):
    cpu = tmp_path / "cpu"
    run([binary("fsm2d"), "--N", 120, "--tau", 3, "--box", 120, "--seed", 2, "--steps", 12, "--out", cpu])
    run([binary("fsm2d_gpu", gpu=True), "--score", cpu, "--frames", 11, "--dump", tmp_path / "g.f32"])
    g = gpu_scores(tmp_path / "g.f32")[0]
    r = load_run(str(cpu))
    for i in range(0, 120, 17):
        ref = agent_scores(r["traj"][11], i, r["params"], tau=3)
        assert np.array_equal(ref < 0, g[i] < 0)
        assert np.allclose(g[i][ref > 0], ref[ref > 0], rtol=2e-3)


def test_splitting_over_gpus_does_not_change_the_run(tmp_path):
    run([binary("fsm2d_gpu", gpu=True), *MODEL, "--steps", 15, "--out", tmp_path / "one"])
    run([binary("fsm2d_gpu", gpu=True), *MODEL, "--steps", 15, "--gpus", "0,0,0", "--out", tmp_path / "three"])
    assert digest(str(tmp_path / "one")) == digest(str(tmp_path / "three"))


def test_resume_is_bit_identical(tmp_path):
    run([binary("fsm2d_gpu", gpu=True), *MODEL, "--steps", 15, "--out", tmp_path / "full"])
    run([binary("fsm2d_gpu", gpu=True), *MODEL, "--steps", 6, "--out", tmp_path / "part"])
    run([binary("fsm2d_gpu", gpu=True), *MODEL, "--steps", 15, "--resume", "--out", tmp_path / "part"])
    assert digest(str(tmp_path / "full")) == digest(str(tmp_path / "part"))
