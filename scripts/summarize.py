"""One line of flock measures per run (after a burn-in), for comparing runs.

  python scripts/summarize.py runs/n500_t6 runs/n1000_t6 ... [--burn-in 200] [--every 10]

Columns: mean polar order; mean share of agents in the largest group (and its minimum) and mean number of groups
(groups are linked at 3 x the median nearest-neighbour distance); median nearest-neighbour distance; radius of
gyration; collisions per agent and step; in 2D, the mean fraction of the visual field that is occupied.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from fsmtools import load_run  # noqa: E402
from fsmtools.metrics import summary  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--burn-in", type=int, default=200)
    ap.add_argument("--every", type=int, default=10)
    a = ap.parse_args()
    print(f"{'run':28s} {'frames':>6s} {'order':>6s} {'largest (min)':>14s} {'groups':>7s} {'nn':>6s} {'Rg':>8s} {'collisions':>10s} {'opacity':>7s}")
    for r in a.runs:
        s = summary(load_run(r, mmap=True), a.burn_in, a.every)
        print(f"{os.path.basename(r):28s} {s['frames']:6d} {s['order']:6.3f} {s['largest_group']:7.3f} ({s['largest_group_min']:.2f}) "
              f"{s['groups']:7.1f} {s['nn_distance']:6.2f} {s['radius_of_gyration']:8.1f} {s['collisions_per_agent_step']:10.5f} "
              f"{s.get('opacity', float('nan')):7.3f}")


if __name__ == "__main__":
    main()
