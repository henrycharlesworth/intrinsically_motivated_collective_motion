"""Bundle runs with the standalone viewer, optionally as one zip to send around.

  python scripts/build_viewer.py runs/n2500_t6 "runs/pred/drive_above=3D · two predators, split and rejoin"
  python scripts/build_viewer.py runs/n2500_t6 runs/n1000_t6_3d --zip fsm_viewer.zip

Each argument is a run prefix, optionally with "=LABEL" (the menu label; default: N, tau and steps). The runs are
written to viewer/data (replacing what is there) in the order given; the first 2D and the first 3D run open by
default. With --zip, the viewer (pages, scripts and these runs) is packed into one file: unzip it anywhere and open
index.html (2D) or 3d.html (3D) in a browser. No server, installation or GPU is needed.
"""
import argparse
import glob
import os
import sys
import zipfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from fsmtools.viewer_export import export_run, write_index  # noqa: E402

VIEWER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "viewer")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="PREFIX or PREFIX=LABEL")
    ap.add_argument("--out", default=os.path.join(VIEWER, "data"))
    ap.add_argument("--zip", help="also pack the viewer and its data into this zip file")
    ap.add_argument("--scores", action="store_true", help="3D: include the per-step scores (for the bird panel; large)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for f in glob.glob(os.path.join(a.out, "*.js")):
        os.remove(f)
    entries = []
    for spec in a.runs:
        prefix, _, label = spec.partition("=")
        e = export_run(prefix, a.out, label=label or None, scores=a.scores)
        size = os.path.getsize(os.path.join(a.out, e["id"] + ".js"))
        print(f"{e['id']:32s} {size / 1e6:7.1f} MB  {e['label']}")
        entries.append(e)
    write_index(a.out, entries)
    if a.zip:
        with zipfile.ZipFile(a.zip, "w", zipfile.ZIP_STORED) as z:  # the data is already compressed
            for rel in ["index.html", "3d.html"] + sorted(os.path.relpath(f, VIEWER) for f in glob.glob(os.path.join(VIEWER, "js", "*.js"))):
                z.write(os.path.join(VIEWER, rel), os.path.join("fsm_viewer", rel))
            for f in sorted(glob.glob(os.path.join(a.out, "*.js"))):
                z.write(f, os.path.join("fsm_viewer", "data", os.path.basename(f)))
        print(f"{a.zip}: {os.path.getsize(a.zip) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
