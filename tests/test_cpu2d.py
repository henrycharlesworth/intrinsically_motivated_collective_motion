"""The exact CPU simulator: its fast score is the definition, the DAG is the full tree, and results do not depend on
the thread count or on interrupting and resuming a run."""
import re

from conftest import binary, digest, run

MODEL = ["--N", 60, "--tau", 4, "--box", 80, "--seed", 3]


def test_fast_score_equals_the_definition(tmp_path):
    log = run([binary("fsm2d"), *MODEL, "--steps", 1, "--verify", "--out", tmp_path / "v"])
    err = float(re.search(r"verify: max relative \|W - W_brute_force\| = (\S+)", log).group(1))
    assert err < 1e-10


def test_dag_gives_the_same_run_as_the_full_tree(tmp_path):
    run([binary("fsm2d"), *MODEL, "--steps", 4, "--out", tmp_path / "dag"])
    run([binary("fsm2d"), *MODEL, "--steps", 4, "--full-tree", "--out", tmp_path / "tree"])
    assert digest(str(tmp_path / "dag")) == digest(str(tmp_path / "tree"))


def test_thread_count_does_not_change_the_run(tmp_path):
    run([binary("fsm2d"), *MODEL, "--steps", 6, "--threads", 1, "--out", tmp_path / "t1"])
    run([binary("fsm2d"), *MODEL, "--steps", 6, "--threads", 4, "--out", tmp_path / "t4"])
    assert digest(str(tmp_path / "t1")) == digest(str(tmp_path / "t4"))


def test_resume_is_bit_identical(tmp_path):
    run([binary("fsm2d"), *MODEL, "--steps", 8, "--out", tmp_path / "full"])
    run([binary("fsm2d"), *MODEL, "--steps", 3, "--out", tmp_path / "part"])
    run([binary("fsm2d"), *MODEL, "--steps", 8, "--resume", "--out", tmp_path / "part"])
    assert digest(str(tmp_path / "full")) == digest(str(tmp_path / "part"))
