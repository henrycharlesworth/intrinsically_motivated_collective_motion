"""Render a video of a run from the viewer, frame by frame: headless Chrome draws each frame (driven through the
viewer's scripting hooks, window.FSM2D / window.FSM3D) and ffmpeg encodes them.

  python viewer/movie/movie.py viewer/index.html --run n10000_t6 viewer/movie/shots/2d_opening.js out.mp4
  python viewer/movie/movie.py viewer/3d.html --run pred_drive_above viewer/movie/shots/3d_predators.js out.mp4
  python viewer/movie/movie.py ... --stills 300,900      # check shots as PNGs (OUT_300.png, ...) before rendering

A shot file sets window.__SHOT = {...} (see directors/2d.js and directors/3d.js for its fields); the matching director
turns it into the camera, speed and overlays of every frame. --gl picks Chrome's WebGL: "gl" (the machine's GPU
through OpenGL; much faster) or "sw" (SwiftShader, CPU only; always works, identical pictures). Needs Google Chrome
or Chromium, and pip install websocket-client imageio-ffmpeg.
"""
import argparse
import base64
import json
import os
import shutil
import signal
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request

import imageio_ffmpeg
import websocket

HERE = os.path.dirname(os.path.abspath(__file__))
GL = {"sw": ["--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
      "gl": ["--use-angle=gl", "--ignore-gpu-blocklist", "--enable-gpu"]}


class Chrome:
    """A headless Chrome tab, driven over the DevTools protocol."""

    def __init__(self, size, gl, port):
        exe = next((shutil.which(c) for c in ("google-chrome", "chromium", "chromium-browser") if shutil.which(c)), None)
        if not exe:
            raise SystemExit("Google Chrome or Chromium is needed")
        self.profile = tempfile.mkdtemp(prefix="fsm-chrome-")
        self.proc = subprocess.Popen([exe, "--headless=new", f"--remote-debugging-port={port}", f"--remote-allow-origins=http://127.0.0.1:{port}",
                                      f"--user-data-dir={self.profile}", "--allow-file-access-from-files", *GL[gl], "--hide-scrollbars",
                                      f"--window-size={size[0]},{size[1]}", "about:blank"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(100):
            try:
                tabs = json.load(urllib.request.urlopen(f"http://127.0.0.1:{port}/json"))
                break
            except OSError:
                time.sleep(0.2)
        else:
            raise SystemExit("Chrome did not start")
        self.ws = websocket.create_connection([t for t in tabs if t["type"] == "page"][0]["webSocketDebuggerUrl"], timeout=600)
        self.id = 0
        self.call("Runtime.enable")
        self.call("Emulation.setDeviceMetricsOverride", width=size[0], height=size[1], deviceScaleFactor=1, mobile=False)

    def call(self, method, **params):
        self.id += 1
        self.ws.send(json.dumps({"id": self.id, "method": method, "params": params}))
        while True:
            m = json.loads(self.ws.recv())
            if m.get("method") == "Runtime.exceptionThrown":
                print("page error:", json.dumps(m["params"]["exceptionDetails"])[:400])
            if m.get("id") == self.id:
                if "error" in m:
                    raise RuntimeError(m["error"])
                return m.get("result", {})

    def js(self, expr):
        r = self.call("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
        if "exceptionDetails" in r:
            raise RuntimeError("in the page: " + json.dumps(r["exceptionDetails"])[:600])
        return r.get("result", {}).get("value")

    def screenshot(self, fmt="jpeg"):
        r = self.call("Page.captureScreenshot", format=fmt, **({"quality": 95} if fmt == "jpeg" else {}))
        return base64.b64decode(r["data"])

    def close(self):
        self.proc.send_signal(signal.SIGTERM)
        self.proc.wait()
        shutil.rmtree(self.profile, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("page", help="viewer/index.html (2D) or viewer/3d.html (3D), or any URL of them")
    ap.add_argument("shot", help="shot file (viewer/movie/shots/*.js)")
    ap.add_argument("out", help="output .mp4 (or the PNG name pattern with --stills)")
    ap.add_argument("--run", help="bundled run to show (the viewer's default otherwise)")
    ap.add_argument("--size", default="1920x1080")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--crf", type=int, default=18, help="H.264 quality (lower: better and larger)")
    ap.add_argument("--gl", choices=sorted(GL), default="gl")
    ap.add_argument("--port", type=int, default=9340)
    ap.add_argument("--stills", default="", help="comma-separated frame numbers: write PNGs instead of a video")
    a = ap.parse_args()
    size = tuple(int(v) for v in a.size.split("x"))
    url = a.page if "://" in a.page else "file://" + os.path.abspath(a.page)
    if a.run:
        url += ("&" if "?" in url else "?") + "run=" + urllib.parse.quote(a.run)
    shot = open(a.shot).read()
    director = open(os.path.join(HERE, "directors", "3d.js" if "3d" in os.path.basename(a.page) else "2d.js")).read()
    chrome = Chrome(size, a.gl, a.port)
    ff = None
    try:
        chrome.call("Page.navigate", url=url)
        time.sleep(1)
        chrome.js(shot + "\n" + director)
        t0 = time.time()
        while not chrome.js("window.__D.ready()"):
            if time.time() - t0 > 600:
                raise SystemExit("the viewer did not finish loading the run")
            time.sleep(0.5)
        info = chrome.js("Promise.resolve(window.__D.setup()).then((x) => JSON.stringify(x))")
        n = chrome.js("window.__D.n")
        print(f"{url}: {info}; {n} frames ({n / a.fps:.1f} s)")
        stills = [int(s) for s in a.stills.split(",") if s]
        if not stills:
            ff = subprocess.Popen([imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error", "-f", "image2pipe", "-framerate", str(a.fps),
                                   "-c:v", "mjpeg", "-i", "-", "-c:v", "libx264", "-preset", "slow", "-crf", str(a.crf), "-pix_fmt", "yuv420p",
                                   "-movflags", "+faststart", a.out], stdin=subprocess.PIPE)
        t0 = time.time()
        for q, i in enumerate(stills or range(n)):
            chrome.js(f"Promise.resolve(window.__D.frame({i}))")
            if stills:
                open(os.path.splitext(a.out)[0] + f"_{i}.png", "wb").write(chrome.screenshot("png"))
            else:
                ff.stdin.write(chrome.screenshot())
            if q % 100 == 0:
                print(f"frame {i}/{n}  {(time.time() - t0) / (q + 1):.2f} s/frame", flush=True)
        print(f"done: {len(stills) or n} frames in {time.time() - t0:.0f} s")
    finally:
        if ff:
            ff.stdin.close()
            ff.wait()
        chrome.close()


if __name__ == "__main__":
    main()
