"""Run a set of simulations on this machine's GPUs, several at a time, resumably.

  python3 scripts/run_matrix.py --dims 2 --N 500,1000,2500 --tau 6,7,8 --steps 5000 --gpus 0-7 --out runs/matrix
  python3 scripts/run_matrix.py --dims 3 --N 1000,10000 --tau 6 --gpus-per-run 1,4 --steps 2000 --gpus 0-7 --out runs/m3d

Every (N, tau) pair is one run, named n{N}_t{tau}. The initial box grows with N at a fixed density: 250 sqrt(N/500)
in 2D (the paper's N = 500 in a 250 x 250 box), 11 N^(1/3) in 3D. Runs start as soon as enough GPUs are free, in
the order N descending, then tau descending (longest first); --gpus-per-run gives each N its number of GPUs (one
value, or one per N). Every run is started with --resume, so rerunning the same command after a crash or a reboot
continues where the runs stopped (bit-identical to uninterrupted runs); a finished run gets a .done marker and is
skipped; a failed run is retried once. Progress goes to OUT/queue.log. Extra simulator flags go after "--".
--poweroff shuts the machine down when everything has finished (for cloud machines). Standard library only.
"""
import argparse
import os
import subprocess
import sys
import threading
import time

BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin")


def parse_gpus(s):
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    return [int(v) for v in s.split(",") if v]


def build_jobs(a, extra):
    given = [int(v) for v in a.N.split(",")]
    per = [int(v) for v in a.gpus_per_run.split(",")]
    per = dict(zip(given, per if len(per) == len(given) else per * len(given)))
    Ns = sorted(given, reverse=True)
    taus = sorted((int(v) for v in a.tau.split(",")), reverse=True)
    exe = os.path.join(BIN, "fsm3d_gpu" if a.dims == 3 else "fsm2d_gpu")
    jobs = []
    for N in Ns:
        for tau in taus:
            name = f"n{N}_t{tau}"
            box = 11 * N ** (1 / 3) if a.dims == 3 else 250 * (N / 500) ** 0.5
            cmd = [exe, "--N", str(N), "--tau", str(tau), "--box", f"{box:.6g}", "--steps", str(a.steps), "--seed", str(a.seed),
                   "--log-every", str(a.log_every), "--out", os.path.join(a.out, name), "--resume"]
            if a.dims == 3:
                cmd += ["--grid", a.grid, "--merge", str(1.5 if tau <= 7 else 2.0)]
            jobs.append({"name": name, "cmd": cmd + extra, "gpus": per[N]})
    return jobs


def main():
    argv = sys.argv[1:]
    extra = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[: argv.index("--")] if "--" in argv else argv
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dims", type=int, choices=[2, 3], default=2)
    ap.add_argument("--N", required=True, help="comma-separated")
    ap.add_argument("--tau", required=True, help="comma-separated")
    ap.add_argument("--steps", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--gpus", default="0", help="GPUs to use: 0-7 or 0,1,3")
    ap.add_argument("--gpus-per-run", default="1", help="one value, or one per N (in the order given)")
    ap.add_argument("--grid", default="192x96", help="3D visual-field grid")
    ap.add_argument("--log-every", type=int, default=50)
    ap.add_argument("--out", default="runs/matrix")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--poweroff", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    jobs = [j for j in build_jobs(a, extra) if not os.path.exists(os.path.join(a.out, j["name"] + ".done"))]
    free = parse_gpus(a.gpus)
    if any(j["gpus"] > len(free) for j in jobs):
        sys.exit("a run needs more GPUs than --gpus gives")
    if a.dry_run:
        for j in jobs:
            print(f"{j['gpus']} GPU(s)  " + " ".join(j["cmd"]))
        return

    cv = threading.Condition()
    qlog = open(os.path.join(a.out, "queue.log"), "a")
    failed, threads = [], []

    def note(msg):
        with cv:
            qlog.write(time.strftime("%Y-%m-%d %H:%M:%S ") + msg + "\n")
            qlog.flush()

    def run(j, gpus):
        for attempt in (1, 2):
            note(f"start {j['name']} on GPUs {gpus} (attempt {attempt})")
            t0 = time.time()
            with open(os.path.join(a.out, j["name"] + ".log"), "a") as lg:
                rc = subprocess.call(j["cmd"] + ["--gpus", ",".join(map(str, gpus))], stdout=lg, stderr=subprocess.STDOUT)
            note(f"{'done' if rc == 0 else f'FAILED (exit {rc})'} {j['name']} in {(time.time() - t0) / 3600:.2f} h")
            if rc == 0:
                open(os.path.join(a.out, j["name"] + ".done"), "w").close()
                break
        else:
            with cv:
                failed.append(j["name"])
        with cv:
            free.extend(gpus)
            free.sort()
            cv.notify_all()

    note(f"queue: {len(jobs)} runs on GPUs {free}")
    for j in jobs:
        with cv:
            while len(free) < j["gpus"]:
                cv.wait()
            gpus, free[:] = free[: j["gpus"]], free[j["gpus"]:]
        t = threading.Thread(target=run, args=(j, gpus))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    note("all finished" + (f"; failed: {', '.join(failed)}" if failed else ""))
    if a.poweroff:
        qlog.close()
        subprocess.call(["sudo", "poweroff"])
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
