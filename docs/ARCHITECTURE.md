<!-- SPDX-License-Identifier: MIT -->
# Architecture

Written for the maintainer opening this a year from now. Plain HTML, CSS and ES modules;
no framework, no build step, no third-party code. Deployment is "copy the folder", so
there is nothing between what is in the repository and what runs.

## The one-paragraph version

A page loads `js/page.js`, which builds the header, resolves which atlas version to show
from `?v=`, and renders the citation footer. `js/data.js` fetches the shipped tables,
lazily. Browse views read published numbers out of those tables and never recompute them.
The sequence views hand the submitted sequence to `js/align.worker.js`, which uses
`js/align.js` to place it on NL4-3 and `js/mra_engine.js` to evaluate the model on the
result. Nothing leaves the origin.

## File by file

| File | What it is responsible for |
|---|---|
| `js/version.js` | The version spine. Parses `?v=`, loads `data/releases.json`, **refuses to start if a bundle's embedded version disagrees with its directory**, verifies SHA-256 where the browser allows, renders the not-latest banner and the version selector, formats the citation and the export provenance header. |
| `js/page.js` | The shell every page shares: header, nav (every link carries `?v=`), banner, footer, privacy panel, standing caveats, service-worker registration. |
| `js/data.js` | Fetching and indexing the shipped assets, each at most once per page. `loadSmallAssets()` is the sequence tool's entire requirement (~360 KB gzipped); the multi-megabyte browse tables load only when a browse view asks. Also a quote-aware CSV parser. |
| `js/mra_engine.js` | **The engine.** No DOM, no I/O, no dependencies. Evaluates the model for one site, enumerates all 63 alternative codons per codon, sums path rates into `mu_aa`, and resolves coordinates (HXB2 labels, Gag/Pol subunits, both directions). Pinned by `test/roundtrip_nl43.mjs`. |
| `js/align.js` | 12-mer seed voting plus a banded affine-gap alignment, query-global and free at both reference ends. Tries both orientations. Widens the band and re-runs if the optimal path touches a band edge that is not also a reference edge. Owns `MAX_BATCH_SEQUENCES`. |
| `js/align.worker.js` | Runs alignment *and* model evaluation off the main thread. Keeps results in worker memory and posts summaries, serving one sequence's full tables on request — a batch of 100 genomes would otherwise push millions of rows across the thread boundary at once. |
| `js/analysis.js` | Drives the worker for both sequence views, with a main-thread fallback for browsers without module workers, and a progress bar. |
| `js/mutations.js` | Parsing and resolving `RT K65R`-style strings. DOM-free so `test/mutation_cases.mjs` can test it. |
| `js/table.js` | Sortable, filterable tables with CSV export and a provenance header. Caps the DOM at 2,000 rows; downloads are complete. |
| `js/chart.js` | The three charts, in plain SVG: `barChart`, `genomeTrack`, `dotPlot`. |
| `js/{index,browse,sequence,compare,download}.js` | One file per page, containing only what is specific to that view. |
| `sw.js` | Offline after first load. Cache-first for `data/<version>/` (immutable by design), network-first for `releases.json` and everything else. |

## Decisions worth knowing about

**No vendored chart or table library.** The brief allowed one. Not using one means the
no-telemetry guarantee is checkable by reading four small files instead of auditing a
minified bundle, and the maintainer can read every line. See `vendor/README.md`.

**Codon frames are read, not derived.** `data/<version>/nl43_reference.json` carries each
protein's codon coordinates taken from the published path table. HIV-1 genes overlap and
*tat*/*rev* are spliced — one *rev* codon spans a 2,325 nt junction. Re-deriving reading
frames would mean special-casing all of that; reading the published assignment means it is
correct by construction. The build script asserts every codon string matches the reference
sequence at its stated coordinates.

**The engine computes rates; published numbers are read.** For the NL4-3 browse views
every displayed number comes out of a shipped file. The engine runs only for submitted
sequences — where there is nothing to read.

**Intervals for custom sequences come from `interval_ratios.json`.** The model bundle
carries posterior means of coefficients and nothing about posterior width, so it cannot
produce an interval on its own. Intervals are carried as multiplicative ratios to the
median: the published per-site ratio where the exact (position, substitution) pair exists,
a (type, 3′ base) mean where it does not. Rows using the approximation are flagged.

**Pooled context cells.** Cells fitted type-only still carry a small non-zero
`beta_sub_context_delta`, because the fit multiplied it by zero rather than removing it.
The engine must not apply it. `sub_context_pooled_cells` in the bundle is what tells it
which cells those are.

## Adding a view

1. New `page.html` with a `<main>` and empty slot `<div id="...">` elements.
2. New `js/page.js` that calls `bootstrapPage("page.html")` first.
3. Add it to `PAGES` in `js/page.js` and to `SHELL` in `sw.js`.
4. Run `python test/check_integrity.py` — it checks that every referenced path exists,
   that every element id a script reaches for is declared, that no off-origin URL has crept
   in, and that every source file carries its SPDX header.
