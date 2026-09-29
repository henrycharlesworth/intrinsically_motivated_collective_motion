"""The viewer's run bundles: the moves-only format reproduces the trajectory, and the bundle decodes."""
import base64
import json
import shutil
import subprocess
import zlib

import numpy as np
import pytest

from conftest import ROOT, binary, run
from fsmtools import load_run
from fsmtools.viewer_export import export_run, replay2d


@pytest.fixture(scope="module")
def run2d(tmp_path_factory):
    d = tmp_path_factory.mktemp("run")
    run([binary("fsm2d"), "--N", 80, "--tau", 4, "--box", 90, "--steps", 30, "--seed", 7, "--out", d / "r"])
    return str(d / "r")


def test_moves_replay_the_trajectory(run2d):
    r = load_run(run2d)
    traj = replay2d(r["params"], r["traj"][0], r["moves"])
    assert np.allclose(traj, r["traj"][:, :, :3], rtol=0, atol=1e-9)


def test_the_viewer_replay_matches_the_simulator(run2d):
    if shutil.which("node") is None:
        pytest.skip("node is not installed")
    r = load_run(run2d)
    js = ("const C = require(process.argv[1]); const fs = require('fs'); const p = JSON.parse(fs.readFileSync(process.argv[2] + '.json'));"
          "const tr = new Float64Array(fs.readFileSync(process.argv[2] + '.traj').buffer.slice(0, 8 * 4 * p.N));"
          "const mv = new Uint8Array(fs.readFileSync(process.argv[2] + '.moves'));"
          "const T = mv.length / p.N + 1, out = C.replay2d(p, tr, mv, T, p.N);"
          "process.stdout.write(Buffer.from(out.buffer));")
    out = subprocess.run(["node", "-e", js, f"{ROOT}/viewer/js/fsm-core.js", run2d], capture_output=True, check=True).stdout
    traj = np.frombuffer(out, np.float32).reshape(-1, r["traj"].shape[1], 3)
    assert np.allclose(traj, r["traj"][:, :, :3], rtol=0, atol=1e-3)  # float32 output


def test_bundle_decodes(run2d, tmp_path):
    e = export_run(run2d, str(tmp_path), run_id="r")
    text = open(tmp_path / "r.js").read()
    assert text.startswith('FSMData.register("r", ') and e["dims"] == 2
    arrays = json.loads(text[len('FSMData.register("r", '):-3])["arrays"]
    a = arrays["moves"]
    moves = np.frombuffer(zlib.decompress(base64.b64decode("".join(a["chunks"])), -15), np.uint8)
    assert np.array_equal(moves, load_run(run2d)["moves"].ravel())
