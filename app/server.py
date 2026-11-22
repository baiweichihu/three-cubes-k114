#!/usr/bin/env python3
"""Local web console for the k114 engine.

Same shape as the Math Gambling console: a single page on 127.0.0.1, one child
process, live progress parsed from the child's stdout. Standard library only.

The engine is *dumb*: it scans `d in [dlo, dhi)` with `|z|/d in [rlo, rhi]` and
emits one JSON line per completed block. Everything stateful lives here:

  * Start  -> spawn the engine with `--resume --emit-json`
  * Stop   -> terminate it (at most one block of work is lost)
  * the set of scanned `d` intervals, the totals, and the plot

State is persisted to `<out>/console.json`, so closing the console and coming
back later keeps the picture.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import threading
import time
import webbrowser
from collections import deque
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STATIC = Path(__file__).resolve().parent / "static"
IS_WINDOWS = os.name == "nt"
EXE = ROOT / "target" / "release" / ("k114.exe" if IS_WINDOWS else "k114")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


class Console:
    """Owns the child process and the accumulated picture of the region."""

    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.lock = threading.Lock()
        self.proc: subprocess.Popen | None = None
        self.reader: threading.Thread | None = None
        self.started_at: float | None = None
        self.session_seconds = 0.0
        self.runs: list[list[int]] = []
        self.current: list[int] | None = None
        self.totals = {
            "points": 0,
            "curves": 0,
            "exact": 0,
            "hits": 0,
            "no_roots": 0,
            "div3": 0,
            "dense": 0,
            "blocks": 0,
        }
        self.hits: list[dict] = []
        self.last_error: str | None = None
        self.last_event = 0.0
        self.completed = False
        self.rate_window: deque = deque(maxlen=64)
        self.state_path = Path(cfg["out"]) / "console.json"
        self.load()

    # ---------------------------------------------------------------- state
    def load(self) -> None:
        try:
            data = json.loads(self.state_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            data = None
        if data and data.get("region") == self.region_signature():
            self.runs = [[int(a), int(b)] for a, b in data.get("runs", [])]
            self.totals.update(data.get("totals", {}))
            self.hits = data.get("hits", [])[:64]
            self.completed = bool(data.get("completed"))
        # Fall back to the engine's own frontier if we have no console state.
        if not self.runs:
            frontier = self.engine_frontier()
            if frontier and frontier > self.cfg["dlo"]:
                self.runs = [[int(self.cfg["dlo"]), int(min(frontier, self.cfg["dhi"]))]]

    def save(self) -> None:
        payload = {
            "region": self.region_signature(),
            "runs": self.runs,
            "totals": self.totals,
            "hits": self.hits[:64],
            "completed": self.completed,
            "saved_at": time.time(),
        }
        try:
            self.state_path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.state_path.with_suffix(".tmp")
            tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            tmp.replace(self.state_path)
        except OSError as exc:  # never let persistence kill the console
            log(f"could not save state: {exc}")

    def region_signature(self) -> dict:
        c = self.cfg
        return {
            "k": c["k"],
            "dlo": int(c["dlo"]),
            "dhi": int(c["dhi"]),
            "rlo": float(c["rlo"]),
            "rhi": float(c["rhi"]),
        }

    def engine_frontier(self) -> int | None:
        path = Path(self.cfg["out"]) / "progress.log"
        best = None
        try:
            with path.open("r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if not line.startswith("d<="):
                        continue
                    token = line[3:].split(None, 1)
                    if not token:
                        continue
                    try:
                        value = int(token[0])
                    except ValueError:
                        continue
                    best = value if best is None else max(best, value)
        except OSError:
            return None
        return best

    # ------------------------------------------------------------- coverage
    def coverage(self) -> dict | None:
        if not self.runs:
            return None
        d0 = self.runs[0][0]
        d1 = self.runs[-1][1]
        rlo = float(self.cfg["rlo"])
        rhi = float(self.cfg["rhi"])
        # H = max(|x|,|y|) = d * R^(3/2) / sqrt(3)
        h = lambda d, r: d * (r ** 1.5) / math.sqrt(3.0)
        return {
            "d": [d0, d1],
            "z": [rlo * d0, rhi * d1],
            "h_lo": h(d0, rlo),
            "h_hi": h(d1, rhi),
            "z_lo": rlo * d0,
            "z_hi": rhi * d1,
            "d_lo": d0,
            "d_hi": d1,
        }

    # -------------------------------------------------------------- process
    def start(self) -> str | None:
        with self.lock:
            if self.proc is not None and self.proc.poll() is None:
                return "already running"
            frontier = self.engine_frontier() or int(self.cfg["dlo"])
            if frontier >= int(self.cfg["dhi"]):
                self.completed = True
                return "region already complete"
            if not EXE.exists():
                return f"engine not built: {EXE}"

            out = Path(self.cfg["out"])
            out.mkdir(parents=True, exist_ok=True)
            err_log = (out / "engine.stderr.log").open("a", encoding="utf-8")

            args = [
                str(EXE),
                "--k", str(self.cfg["k"]),
                "--dlo", str(int(self.cfg["dlo"])),
                "--dhi", str(int(self.cfg["dhi"])),
                "--rlo", repr(float(self.cfg["rlo"])),
                "--rhi", repr(float(self.cfg["rhi"])),
                "--threads", str(int(self.cfg["threads"])),
                "--block", str(int(self.cfg["block"])),
                "--out", str(out),
                "--resume",
                "--emit-json",
            ]
            flags = 0
            if IS_WINDOWS and self.cfg.get("priority") == "below_normal":
                flags = getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0)
            try:
                self.proc = subprocess.Popen(
                    args,
                    stdout=subprocess.PIPE,
                    stderr=err_log,
                    stdin=subprocess.DEVNULL,
                    text=True,
                    encoding="utf-8",
                    bufsize=1,
                    creationflags=flags,
                )
            except OSError as exc:
                err_log.close()
                self.last_error = str(exc)
                return f"could not start engine: {exc}"

            self.started_at = time.time()
            self.session_seconds = 0.0
            self.last_error = None
            self.completed = False
            self.rate_window.clear()
            self.reader = threading.Thread(target=self._pump, args=(self.proc, err_log), daemon=True)
            self.reader.start()
            log(f"started pid={self.proc.pid} from d={frontier}")
            return None

    def stop(self) -> str | None:
        with self.lock:
            proc = self.proc
            if proc is None or proc.poll() is not None:
                return "not running"
        log(f"stopping pid={proc.pid}")
        try:
            proc.terminate()
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
        except OSError:
            pass
        with self.lock:
            if self.started_at is not None:
                self.session_seconds = time.time() - self.started_at
            self.started_at = None
            self.current = None
            self.proc = None
        self.save()
        return None

    def _pump(self, proc: subprocess.Popen, err_log) -> None:
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                line = line.strip()
                if not line.startswith("{"):
                    continue
                try:
                    event = json.loads(line)
                except ValueError:
                    continue
                self._apply(event)
        finally:
            try:
                err_log.close()
            except OSError:
                pass
            code = proc.poll()
            with self.lock:
                if self.started_at is not None:
                    self.session_seconds = time.time() - self.started_at
                self.started_at = None
                self.current = None
                self.proc = None
            if code not in (0, None):
                self.last_error = f"engine exited with code {code}"
            log(f"engine stopped (exit {code})")
            self.save()

    def _apply(self, event: dict) -> None:
        kind = event.get("event")
        now = time.time()
        with self.lock:
            self.last_event = now
            if kind == "start":
                self.current = None
            elif kind == "block":
                lo = int(event["d_lo"])
                hi = int(event["d_hi"])
                if self.runs and self.runs[-1][1] == lo:
                    self.runs[-1][1] = hi
                else:
                    self.runs.append([lo, hi])
                self.current = [lo, hi]
                self.totals["points"] += int(event.get("points", 0))
                self.totals["curves"] += int(event.get("curves", 0))
                self.totals["exact"] += int(event.get("exact", 0))
                self.totals["hits"] += int(event.get("hits", 0))
                self.totals["no_roots"] += int(event.get("no_roots", 0))
                self.totals["div3"] += int(event.get("div3", 0))
                self.totals["dense"] += int(event.get("dense", 0))
                self.totals["blocks"] += 1
                self.rate_window.append((now, hi - lo))
            elif kind == "hit":
                self.hits.insert(0, event)
                del self.hits[64:]
                self.totals["hits"] += 0  # counted via blocks; this is for display
            elif kind == "done":
                self.completed = bool(event.get("completed"))
                self.current = None
            elif kind == "complete":
                self.completed = True
                self.current = None

    # ----------------------------------------------------------------- view
    def state(self) -> dict:
        with self.lock:
            running = self.proc is not None and self.proc.poll() is None
            elapsed = (
                (time.time() - self.started_at)
                if (running and self.started_at)
                else self.session_seconds
            )
            # Throughput over the recent window.
            rate_d = 0.0
            if len(self.rate_window) >= 2:
                (t0, _), (t1, _) = self.rate_window[0], self.rate_window[-1]
                span = max(t1 - t0, 1e-9)
                rate_d = sum(d for _, d in self.rate_window) / span
            cov = self.coverage()
            return {
                "running": running,
                "pid": self.proc.pid if running else None,
                "region": self.region_signature(),
                "totals": dict(self.totals),
                "runs": [list(r) for r in self.runs],
                "current": self.current,
                "hits": self.hits[:20],
                "coverage": cov,
                "elapsed": elapsed,
                "rate_d_per_s": rate_d,
                "completed": self.completed,
                "last_error": self.last_error,
                "engine": str(EXE),
                "priority": self.cfg.get("priority", "normal"),
            }


class Handler(BaseHTTPRequestHandler):
    console: Console

    def log_message(self, fmt, *args):  # keep the terminal quiet
        pass

    def _send(self, code: int, body: bytes, ctype: str) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code: int = 200) -> None:
        self._send(code, json.dumps(obj).encode("utf-8"), "application/json; charset=utf-8")

    def do_GET(self) -> None:
        path = self.path.split("?", 1)[0]
        if path in ("/", "/index.html"):
            try:
                body = (STATIC / "index.html").read_bytes()
            except OSError:
                self._send(500, b"index.html missing", "text/plain; charset=utf-8")
                return
            self._send(200, body, "text/html; charset=utf-8")
            return
        if path == "/api/state":
            self._json(self.console.state())
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")

    def do_POST(self) -> None:
        path = self.path.split("?", 1)[0]
        if path == "/api/start":
            err = self.console.start()
            self._json({"ok": err is None, "error": err}, 200 if err is None else 409)
            return
        if path == "/api/stop":
            err = self.console.stop()
            self._json({"ok": err is None, "error": err})
            return
        self._send(404, b"not found", "text/plain; charset=utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="Local web console for the k114 engine")
    ap.add_argument("--k", type=int, default=114)
    ap.add_argument("--dlo", type=str, required=True, help="lower edge of the d interval, u64 or 1.2e18")
    ap.add_argument("--dhi", type=str, required=True, help="upper edge of the d interval")
    ap.add_argument("--rlo", type=float, default=3.4, help="lower edge of R = |z|/d")
    ap.add_argument("--rhi", type=float, default=64.0, help="upper edge of R = |z|/d")
    ap.add_argument("--threads", type=int, default=2, help="engine worker threads (keep low on a shared PC)")
    ap.add_argument("--block", type=int, default=1 << 16, help="d values per checkpointed block")
    ap.add_argument("--out", type=str, default=None, help="output directory (default runs/<k>-<dlo>)")
    ap.add_argument("--port", type=int, default=8731)
    ap.add_argument("--priority", choices=["normal", "below_normal"], default="below_normal")
    ap.add_argument("--open", action="store_true", help="open the browser on start")
    args = ap.parse_args()

    def as_u64(text: str) -> int:
        try:
            return int(text, 0)
        except ValueError:
            return int(float(text))

    dlo, dhi = as_u64(args.dlo), as_u64(args.dhi)
    if dlo <= 0 or dhi <= dlo:
        ap.error("need 0 < dlo < dhi")
    out = args.out or f"runs/k{args.k}-{dlo:.3e}".replace("+", "")

    cfg = {
        "k": args.k,
        "dlo": dlo,
        "dhi": dhi,
        "rlo": args.rlo,
        "rhi": args.rhi,
        "threads": args.threads,
        "block": args.block,
        "out": out,
        "priority": args.priority,
    }

    console = Console(cfg)
    Handler.console = console
    httpd = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    log(f"region k={args.k} d in [{dlo}, {dhi}) R in [{args.rlo}, {args.rhi}]")
    log(f"engine {EXE}")
    log(f"console ready at {url} (Ctrl+C to quit)")
    if args.open:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        if console.proc is not None and console.proc.poll() is None:
            console.stop()
        console.save()
        httpd.server_close()
        log("bye")
    return 0


if __name__ == "__main__":
    sys.exit(main())
