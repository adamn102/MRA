#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""
check_integrity.py -- static checks on the deployable folder.

    python test/check_integrity.py

Catches the class of mistake that a syntax check cannot: a stylesheet or module path
that does not exist, an element id a page script reaches for that the HTML never
defines, an external URL that would break the privacy guarantee, a source file missing
its SPDX header, and a manifest checksum that no longer matches the file on disk.

Deployment is "copy the folder", so there is no build step to catch any of this.
"""
import hashlib, json, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
problems = []
notes = []


def add(msg):
    problems.append(msg)


# --- 1. referenced files exist -------------------------------------------------
for html in sorted(ROOT.glob("*.html")):
    text = html.read_text()
    for ref in re.findall(r'(?:href|src)="([^"#?]+)"', text):
        if ref.startswith(("http://", "https://", "mailto:", "data:", "//")):
            continue
        if not (ROOT / ref).exists():
            add(f"{html.name}: references {ref}, which does not exist")

for js in sorted(ROOT.glob("js/*.js")) + sorted(ROOT.glob("*.js")) + sorted(ROOT.glob("test/*.mjs")):
    for imp in re.findall(r'from\s+"(\.{1,2}/[^"]+)"', js.read_text()):
        if not (js.parent / imp).resolve().exists():
            add(f"{js.relative_to(ROOT)}: imports {imp}, which does not exist")

# --- 2. element ids a page script uses must exist in its page -------------------
for html in sorted(ROOT.glob("*.html")):
    text = html.read_text()
    m = re.search(r'src="js/([^"]+)"', text)
    if not m:
        continue
    script = ROOT / "js" / m.group(1)
    if not script.exists():
        continue
    ids = set(re.findall(r'id="([^"]+)"', text))
    # ids the shell creates at runtime rather than declaring in the HTML
    ids |= {"version-banner-slot", "mra-version-select"}
    for used in sorted(set(re.findall(r'getElementById\("([^"]+)"\)', script.read_text()))):
        if used not in ids:
            add(f"{script.name}: uses #{used}, which {html.name} does not define")

# --- 3. no off-origin requests -------------------------------------------------
# The privacy guarantee is that nothing a user pastes can leave the machine, and that
# merely using the tool is not reported anywhere. A single CDN script tag or analytics
# beacon would break it, so this check is a release gate, not a style preference.
ALLOWED_URL_CONTEXTS = ("creativecommons.org", "www.w3.org")   # licence link, SVG namespace
for path in sorted(list(ROOT.glob("*.html")) + list(ROOT.glob("js/*.js")) + list(ROOT.glob("*.js")) + list(ROOT.glob("css/*.css"))):
    for url in re.findall(r'https?://[^\s"\'()<>]+', path.read_text()):
        if any(a in url for a in ALLOWED_URL_CONTEXTS):
            continue
        add(f"{path.relative_to(ROOT)}: contains an off-origin URL {url}")
for path in sorted(list(ROOT.glob("js/*.js")) + list(ROOT.glob("*.js"))):
    text = path.read_text()
    for bad in ["XMLHttpRequest", "navigator.sendBeacon", "new Image("]:
        if bad in text:
            add(f"{path.name}: uses {bad}, which can make an off-origin request")
    for fetched in re.findall(r"fetch\(([^)]*)\)", text):
        if "http" in fetched:
            add(f"{path.name}: fetch() with an absolute URL: {fetched.strip()}")

# --- 4. SPDX headers -----------------------------------------------------------
for path in sorted(list(ROOT.glob("js/*.js")) + list(ROOT.glob("*.js")) + list(ROOT.glob("*.html")) +
                   list(ROOT.glob("css/*.css")) + list(ROOT.glob("scripts/*.py")) +
                   list(ROOT.glob("test/*.py")) + list(ROOT.glob("test/*.mjs"))):
    head = path.read_text()[:400]
    if "SPDX-License-Identifier: MIT" not in head:
        add(f"{path.relative_to(ROOT)}: missing SPDX header")

# --- 5. manifest checksums still describe the shipped files --------------------
releases = json.loads((ROOT / "data" / "releases.json").read_text())
for rel in releases["releases"]:
    vdir = ROOT / "data" / rel["version"]
    if not vdir.exists():
        add(f"releases.json lists {rel['version']} but data/{rel['version']}/ does not exist")
        continue
    manifest = json.loads((vdir / "manifest.json").read_text())
    for name, meta in manifest["assets"].items():
        p = vdir / name
        if not p.exists():
            add(f"{rel['version']}: manifest lists {name}, which is missing")
            continue
        actual = hashlib.sha256(p.read_bytes()).hexdigest()
        if actual != meta["sha256"]:
            add(f"{rel['version']}: {name} does not match its manifest checksum")
    for name, meta in (rel.get("assets") or {}).items():
        if manifest["assets"].get(name, {}).get("sha256") != meta["sha256"]:
            add(f"{rel['version']}: releases.json and manifest.json disagree about {name}")
    bundle = json.loads((vdir / "mra_web_model.json").read_text())
    if bundle["version"] != rel["version"]:
        add(f"{rel['version']}: bundle declares version {bundle['version']}")
    notes.append(f"{rel['version']}: {len(manifest['assets'])} assets, "
                 f"{manifest['total_gzip_bytes']/1048576:.2f} MB gzipped, checksums match")

# --- 6. no third-party code ----------------------------------------------------
vendored = [p for p in (ROOT / "vendor").rglob("*") if p.suffix in (".js", ".css", ".mjs")]
notes.append(f"vendor/: {len(vendored)} third-party file(s)")


# --- 7. no collaborator email list committed into the served tree ---------------
# The whitelist belongs in the Cloudflare Access dashboard, never under the served
# folder: everything here is a static file and would leak collaborator addresses to
# anyone past the gate (or, before the gate exists, to the public). This is a backstop
# for CLOUDFLARE_SETUP.md's warning. The contact-email convention (a single address in a
# request) is fine; a LIST of addresses in a served file is the hazard.
import re as _re
EMAIL = _re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
SUSPECT_NAME = _re.compile(r"(whitelist|allowlist|allow-list|collaborator|access[-_]?list|emails?)", _re.I)
SERVED_EXT = {".html", ".js", ".css", ".csv", ".json", ".txt", ".tsv", ".md", ".fasta"}
for path in sorted(ROOT.rglob("*")):
    if path.is_dir() or "/test/" in str(path) or path.suffix.lower() not in SERVED_EXT:
        continue
    rel = path.relative_to(ROOT)
    name_hit = bool(SUSPECT_NAME.search(path.stem))
    try:
        text = path.read_text(errors="ignore")
    except Exception:
        continue
    emails = set(EMAIL.findall(text))
    # a documentation placeholder like name@example.org / @yourlab.org is not a real list
    real = {e for e in emails if not _re.search(r"(example\.|yourlab\.|your-?domain|placeholder)", e, _re.I)}
    if name_hit and real:
        add(f"{rel}: a file named like a collaborator list contains email addresses "
            f"and is in the SERVED tree; move the whitelist to the Cloudflare dashboard")
    elif len(real) >= 3:
        add(f"{rel}: contains {len(real)} real email addresses in the served tree; "
            f"if this is a collaborator whitelist, it must not be deployed")

print("\n".join("  " + n for n in notes))
if problems:
    print(f"\n{len(problems)} PROBLEM(S):")
    for p in problems:
        print("  - " + p)
    sys.exit(1)
print("\nINTEGRITY OK")
