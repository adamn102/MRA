#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
serve_dev.py -- serve the portal locally for manual testing.

    python scripts/serve_dev.py            # http://localhost:8000/
    python scripts/serve_dev.py --port 8080
    python scripts/serve_dev.py --open     # also open a browser window

WHY A SCRIPT RATHER THAN `python -m http.server`
    Three things the stock server gets wrong for this portal:

    1. MIME types. A .js served as text/plain is refused by the browser's module
       loader, and .mjs is not in Python's default map at all. Both are fixed here.
    2. Caching. The stock server sends Last-Modified and the browser will happily
       serve a stale module after you edit it, which during manual testing looks like
       "my change did nothing". Code and pages go out with no-store.
    3. The service worker. sw.js caches aggressively on purpose -- localhost counts as
       a secure context, so it registers during development too, and then serves you
       yesterday's JavaScript. This server sends no-store for sw.js as well, and the
       page skips registration entirely when you add ?nosw=1 to the URL.

    Files under data/<version>/ are still sent cacheable: a release directory is
    immutable by design, so caching them is what production does and nothing is gained
    by fighting it. Rebuild with --force and restart if you regenerate them.

This is a development convenience. In production the portal is just a folder of static
files behind any web server; it needs nothing from this script.
"""
import argparse
import functools
import http.server
import pathlib
import socketserver
import sys
import webbrowser

ROOT = pathlib.Path(__file__).resolve().parent.parent

EXTRA_TYPES = {
    ".js": "text/javascript",
    ".mjs": "text/javascript",
    ".json": "application/json",
    ".css": "text/css",
    ".csv": "text/csv",
    ".fasta": "text/plain",
    ".map": "application/json",
}


class DevHandler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **EXTRA_TYPES}

    def end_headers(self):
        path = self.path.split("?", 1)[0]
        immutable = "/data/v" in path and not path.endswith("releases.json")
        if immutable:
            self.send_header("Cache-Control", "public, max-age=3600")
        else:
            # Code, pages and releases.json: never cached, so an edit is visible on reload.
            self.send_header("Cache-Control", "no-store, must-revalidate")
            self.send_header("Pragma", "no-cache")
            self.send_header("Expires", "0")
        # The production guarantee is that nothing a user submits leaves the machine.
        # Mirroring a restrictive CSP here means a stray off-origin reference fails
        # loudly during development rather than silently working and shipping.
        self.send_header("Content-Security-Policy",
                         "default-src 'self'; script-src 'self'; style-src 'self'; "
                         "img-src 'self' data:; connect-src 'self'; worker-src 'self'; "
                         "form-action 'none'; base-uri 'self'")
        self.send_header("Referrer-Policy", "no-referrer")
        super().end_headers()

    def log_message(self, fmt, *args):
        sys.stderr.write("  %s\n" % (fmt % args))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--host", default="127.0.0.1",
                    help="default 127.0.0.1; use 0.0.0.0 only if you really want it reachable "
                         "from the network, which for a tool that handles patient sequences "
                         "you probably do not")
    ap.add_argument("--open", action="store_true", help="open a browser window")
    args = ap.parse_args()

    if not (ROOT / "index.html").exists():
        sys.exit(f"ERROR: {ROOT} does not look like the portal folder (no index.html).")

    handler = functools.partial(DevHandler, directory=str(ROOT))
    socketserver.TCPServer.allow_reuse_address = True
    url = f"http://{'localhost' if args.host == '127.0.0.1' else args.host}:{args.port}/"
    try:
        httpd = socketserver.TCPServer((args.host, args.port), handler)
    except OSError as exc:
        sys.exit(f"ERROR: could not bind {args.host}:{args.port} -- {exc}\n"
                 f"Another server may already be running; try --port {args.port + 1}.")

    print(f"Serving {ROOT}")
    print(f"  {url}")
    print(f"  {url}?nosw=1      (same site, service worker not registered -- use this while editing)")
    print("Ctrl-C to stop.\n")
    if args.open:
        webbrowser.open(url + "?nosw=1")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        httpd.server_close()


if __name__ == "__main__":
    main()
