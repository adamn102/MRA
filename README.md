<!-- SPDX-License-Identifier: MIT -->
# HIV-1 Mutation Rate Atlas — web portal

A static, client-side portal for the HIV-1 Mutation Rate Atlas (MRA). It does two things:

1. **Browse, query and download** the published NL4-3 atlas — position-, mutation-type- and
   context-resolved de novo mutation rates, and the amino-acid mutabilities derived from them.
2. **Generate a custom atlas for any submitted HIV-1 sequence**, by applying the fitted MRA
   model to that sequence's own codons and local sequence contexts.

**Submitted sequences never leave the user's machine.** There is no server, no database, no
analytics, no telemetry, no error reporter and no third-party code. Alignment and every rate
calculation run in the browser. That is the reason for the static architecture, not a
side-effect of it.

## Status

Atlas **v0.1** — pre-publication, in-lab build. The NL4-3 round-trip test passes; see
[docs/validation_report.md](docs/validation_report.md). **v0.1 does not go to anyone outside
the lab until that report says PASSED for the assets being shipped.**

`doi` is `null` until one is minted. Exports say `[not yet assigned]` rather than inventing an
identifier.

## Deploy

Copy the folder. That is the whole deployment.

```bash
rsync -a mra-portal/ user@host:/var/www/mra/      # or: push to a gh-pages branch
```

No build step, no bundler, no `npm install`, no TypeScript. Plain HTML, CSS and ES modules.

## Run it locally

ES modules and service workers do not work from `file://`, so the folder has to be served
over HTTP. There is a development server in the repository that fixes the three things the
stock `python -m http.server` gets wrong here — JavaScript MIME types, caching of edited
modules, and the service worker caching your old code on localhost:

```bash
python scripts/serve_dev.py --open
# or: python scripts/serve_dev.py --port 8080
```

Then open **<http://localhost:8000/?nosw=1>**.

**Use `?nosw=1` while you are editing.** localhost counts as a secure context, so `sw.js`
registers during development too and will happily serve you yesterday's JavaScript. The
flag skips registration and unregisters any worker already installed for this origin.
Drop the flag when you specifically want to test offline behaviour: load the site once
without it, then turn off the network and reload.

Any static server works — `python -m http.server`, `npx serve`, nginx, GitHub Pages. The
dev script is a convenience, not a dependency.

### What to click first

1. `index.html` — the overview, the privacy statement and the coverage note.
2. `browse.html` — it opens on **PR 46** by default. Change the protein selector to
   **RT** and the position to **184**; the bar chart should put V at the top with a
   single-substitution flag. Paste-a-list is below, pre-filled with seven real
   drug-resistance mutations.
3. `drm.html` — the resistance & escape page. Pick a drug class (try **INI**), then
   click **codon sensitivity** on a glycine DRM such as **IN G140R** to see the 68× spread
   across synonymous codons with NL4-3's marked. The bNAb section below carries a
   provisional banner; its glycan rows show summed sequon-loss.
4. `sequence.html` — press **Load NL4-3 (worked example)** then **Compute custom atlas**.
   This is the round-trip test running in your browser: 2,984 codons, 62,518 amino-acid
   rows, 0 positions needing the spline fallback, 405 dropped for lying past 8,213.
   Then try `test/fixtures/JRFL_pol.fasta` via **Load FASTA file…** for a genuinely
   divergent sequence.
5. `compare.html` — it carries over the sequence from the previous page in the same tab.
   With JR-FL loaded, the first table is the codon-difference view; the genome track at
   the bottom defaults to G→A.
5. `download.html` — every asset with its SHA-256, verified before the download.

Open the browser console on each page. It should stay empty.

## Rebuild the data assets from the source tables

```bash
conda activate webportal     # nodejs + pandas, numpy, scipy, patsy, openpyxl
python scripts/build_web_assets.py \
    --version v0.1 \
    --release-date 2026-10-28 \
    --description "First release." \
    --src   /path/to/MutRatePaper/data/models/mra/production \
    --ref-root /path/to/MutRatePaper \
    --sdrd  /path/to/MutRatePaper/data/mutation_supply/SDRD
```

`--sdrd` adds the drug-resistance catalog and the provisional bNAb assets; omit it to build
the browse/sequence portal alone.

This runs `scripts/extract_mra_web_model.py` (the single authority on model extraction, whose
internal assertions are a hard stop), writes every asset into `data/<version>/`, computes a
SHA-256 manifest, checks the payload budget, and regenerates `data/releases.json`. It
**refuses** to write into a release directory that already carries a manifest — a published
release is never modified in place.

## Run the tests

```bash
node   test/roundtrip_nl43.mjs  --version v0.1   # the decisive test: engine vs published atlas
python test/compare_roundtrip.py --version v0.1
node   test/align_cases.mjs     --version v0.1   # alignment behaviour on real submission shapes
node   test/mutation_cases.mjs  --version v0.1   # mutation parsing and HXB2/subunit coordinates
node   test/codon_sensitivity_cases.mjs --version v0.1  # codon sensitivity + glycan sequon loss
python test/check_integrity.py                   # paths, element ids, off-origin URLs, checksums
python test/smoke_http.py                        # serve the folder, fetch every referenced URL
node   test/align_timing.mjs                     # performance, sets the batch cap
```

`smoke_http.py` needs permission to bind a local socket; it skips with a message where
that is not allowed, and `check_integrity.py` checks the same references statically.

Re-run all of them after any change to the engine, the aligner, the extractor or the source
tables. Treat an assertion failure as a hard stop, not a tolerance to widen.

## Cut a new release

See [docs/VERSIONING.md](docs/VERSIONING.md). In short: a new model is a **new directory**,
never an overwrite, because a link printed in a paper must resolve to the exact atlas version
that paper used, forever.

## Layout

```
index.html  browse.html  sequence.html  compare.html  download.html
css/mra.css                 one stylesheet
js/                         see docs/ARCHITECTURE.md
sw.js                       offline support after first load
data/releases.json          stable path, lists every release
data/v0.1/                  immutable release directory
scripts/                    asset build + model extraction
test/                       the regression suite
docs/                       validation report, caveats, data dictionary, versioning
figures/                    QC and demonstration figures
vendor/                     deliberately empty; see vendor/README.md
```

## Licence

- **Code** (`js/`, `css/`, `scripts/`, `test/`, `*.html`, `sw.js`): MIT — see `LICENSE`.
- **Atlas data** (`data/**`) and every file the portal exports: **CC-BY-4.0** — see
  `LICENSE-DATA`. Attribution is a licence condition.

Every page and every exported file carries a citation block naming the paper, the atlas
version and its DOI.

## Read these before trusting a number

- [docs/validation_report.md](docs/validation_report.md) — what was tested, the measured
  tolerances, and what the round-trip does *not* establish.
- [docs/CAVEATS.md](docs/CAVEATS.md) — the scientific limits of a custom atlas.
- [docs/data_dictionary.md](docs/data_dictionary.md) — every column, and the
  `bio_rate_median` versus `bio_rate_mean` trap.
