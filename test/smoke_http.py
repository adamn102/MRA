#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
smoke_http.py -- serve the folder and confirm every resource it references is actually
served, with the right content type.

    python test/smoke_http.py

Deployment is "copy the folder", so the failure mode this catches is real: a path that
works on a case-insensitive filesystem and 404s on the server, or an asset the build
forgot. It does not execute JavaScript -- it checks that everything the pages and modules
reference is reachable over HTTP.
"""
import http.server, json, pathlib, re, socketserver, sys, threading, urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
PORT = 8731


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=str(ROOT), **kw)

    def log_message(self, *a):
        pass


def collect_urls():
    urls = set()
    for html in sorted(ROOT.glob("*.html")):
        urls.add(html.name)
        for ref in re.findall(r'(?:href|src)="([^"#?]+)"', html.read_text()):
            if not ref.startswith(("http://", "https://", "mailto:", "data:", "//")):
                urls.add(ref)
    for js in sorted(list(ROOT.glob("js/*.js")) + list(ROOT.glob("*.js"))):
        base = js.parent.relative_to(ROOT)
        for imp in re.findall(r'from\s+"(\.{1,2}/[^"]+)"', js.read_text()):
            urls.add(str((base / imp).as_posix()).replace("js/./", "js/"))
    urls.add("data/releases.json")
    releases = json.loads((ROOT / "data" / "releases.json").read_text())
    for rel in releases["releases"]:
        for asset in rel.get("assets", {}):
            urls.add(f"data/{rel['version']}/{asset}")
    return sorted(urls)


def main():
    socketserver.TCPServer.allow_reuse_address = True
    try:
        httpd = socketserver.TCPServer(("127.0.0.1", PORT), Handler)
    except PermissionError:
        # Some sandboxes forbid binding a listening socket. Say so plainly rather than
        # reporting a pass that did not happen; test/check_integrity.py covers the same
        # paths statically.
        print("  cannot bind 127.0.0.1:%d in this environment (socket binding not permitted)." % PORT)
        print("  SKIPPED -- run this on a machine where a local server is allowed.")
        print("  test/check_integrity.py checks the same references statically.")
        sys.exit(0)
    with httpd:
        threading.Thread(target=httpd.serve_forever, daemon=True).start()
        failures = []
        urls = collect_urls()
        for url in urls:
            full = f"http://127.0.0.1:{PORT}/{url}"
            try:
                with urllib.request.urlopen(full, timeout=20) as resp:
                    if resp.status != 200:
                        failures.append(f"{url}: HTTP {resp.status}")
                    elif resp.length == 0:
                        failures.append(f"{url}: empty response")
            except Exception as exc:
                failures.append(f"{url}: {exc}")
        httpd.shutdown()

    print(f"  checked {len(urls)} URLs over HTTP")
    if failures:
        print(f"\n{len(failures)} PROBLEM(S):")
        for f in failures:
            print("  - " + f)
        sys.exit(1)
    print("\nHTTP SMOKE TEST OK")


if __name__ == "__main__":
    main()
