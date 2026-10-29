<!-- SPDX-License-Identifier: MIT -->
# vendor/ — deliberately empty

There is no third-party code in this portal.

The brief allowed a small vendored charting library and a small table library. They were
not used, for two reasons:

1. **Privacy is the point of the architecture.** Users paste patient-derived HIV-1
   sequences into this tool. Every byte of third-party code is something a reviewer
   would have to audit before believing the no-telemetry guarantee. With nothing
   vendored, the guarantee is checkable by reading the four files in `js/`.

2. **The maintainer has to read this in a year.** A 60 KB minified bundle is not
   readable. The charts this portal needs are a horizontal bar chart, a genome-position
   line track and a scatter; `js/chart.js` draws all three in plain SVG in about 300
   lines of commented code, and `js/table.js` does sorting, filtering and CSV export in
   about 200.

If a future view genuinely needs a plotting library, vendor it *here* as a static file
with its licence alongside — never as a CDN `<script>` tag, which would leak the fact
that a user is using the tool and could carry sequence data in a URL.
