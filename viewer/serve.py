"""Serve the viewer and a folder of simulator runs on this machine, so every run in it can be opened from the menu.

  python viewer/serve.py                    # runs/ (the default --out location), http://localhost:8000
  python viewer/serve.py --runs ~/results --port 8080

Any run written by the simulators under the folder (PREFIX.json + PREFIX.traj [+ .moves, .pred, .seen], in any
subfolder) appears in the menu of index.html (2D) and 3d.html (3D), listed by its path. Runs still being written are
fine: the menu and each run are read afresh when the page is reloaded. Needs only Python 3; nothing is sent anywhere
(the server listens on localhost only).
"""
import argparse
import http.server
import json
import os
import re
import urllib.parse
import webbrowser

VIEWER = os.path.dirname(os.path.abspath(__file__))
RUN_FILES = (".json", ".traj", ".moves", ".pred", ".seen")


def find_runs(root):
    """The runs under root, as the viewer's menu entries (id = path of the prefix relative to root)."""
    out = []
    for d, subdirs, files in os.walk(root):
        subdirs.sort()
        for f in sorted(files):
            if not f.endswith(".json") or f[:-5] + ".traj" not in files:
                continue
            prefix = os.path.join(d, f[:-5])
            try:
                with open(prefix + ".json") as fh:
                    p = json.load(fh)
                N, tau = int(p["N"]), p["tau"]
            except (OSError, ValueError, KeyError, TypeError):
                continue  # not a simulator's parameter file
            dims = 3 if p.get("dims", 2) == 3 else 2
            frames = os.path.getsize(prefix + ".traj") // (8 * N * (6 if dims == 3 else 4))
            steps = frames - 1
            if os.path.exists(prefix + ".moves"):
                steps = min(steps, os.path.getsize(prefix + ".moves") // N)
            if steps < 1:
                continue
            rel = os.path.relpath(prefix, root).replace(os.sep, "/")
            head = ("3D · " if dims == 3 else "") + f"{rel} · N = {N} · τ = {tau}"
            out.append({"id": rel, "dims": dims, "label": f"{head} · {steps} steps", "src": "server",
                        "N": N, "steps": steps, "files": [e[1:] for e in RUN_FILES if os.path.exists(prefix + e)]})
    return out


class Handler(http.server.SimpleHTTPRequestHandler):
    """The viewer's files at /, the runs folder at /runs/ (run files only), its listing at /runs/index.json. Byte
    ranges are honoured (the viewer reads only the first frame of a trajectory)."""

    runs_root = "."

    def translate_path(self, path):
        path = urllib.parse.unquote(urllib.parse.urlsplit(path).path)
        if path.startswith("/runs/"):
            rel = os.path.normpath(path[len("/runs/"):]).lstrip(os.sep)
            if rel.startswith("..") or not rel.endswith(RUN_FILES):
                return ""
            return os.path.join(self.runs_root, rel)
        return super().translate_path(path)

    def do_GET(self):
        if urllib.parse.urlsplit(self.path).path == "/runs/index.json":
            body = json.dumps(find_runs(self.runs_root), ensure_ascii=False).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)
            return
        m = re.fullmatch(r"bytes=(\d+)-(\d*)", self.headers.get("Range", ""))
        path = self.translate_path(self.path)
        if m and os.path.isfile(path):
            size = os.path.getsize(path)
            start = int(m.group(1))
            end = min(int(m.group(2)) if m.group(2) else size - 1, size - 1)
            if start > end:
                self.send_error(416)
                return
            self.send_response(206)
            self.send_header("Content-Type", "application/octet-stream")
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            self.send_header("Content-Length", str(end - start + 1))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            with open(path, "rb") as f:
                f.seek(start)
                self.wfile.write(f.read(end - start + 1))
            return
        super().do_GET()

    def end_headers(self):
        if not self.path.startswith("/runs/"):
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def log_message(self, fmt, *args):
        pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(VIEWER, "..", "runs"), help="folder of runs (default: runs/)")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--no-browser", action="store_true", help="do not open a browser tab")
    a = ap.parse_args()
    Handler.runs_root = os.path.abspath(os.path.expanduser(a.runs))
    if not os.path.isdir(Handler.runs_root):
        ap.error(f"no folder {Handler.runs_root}")
    handler = lambda *args, **kw: Handler(*args, directory=VIEWER, **kw)  # noqa: E731
    server = http.server.ThreadingHTTPServer(("127.0.0.1", a.port), handler)
    url = f"http://localhost:{a.port}/"
    print(f"serving {Handler.runs_root} ({len(find_runs(Handler.runs_root))} runs) at {url}  (Ctrl-C stops)", flush=True)
    if not a.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
