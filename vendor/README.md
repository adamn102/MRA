<!-- SPDX-License-Identifier: MIT -->
# vendor/

Charts and tables are drawn by this portal, in `js/chart.js` and `js/table.js`. They are
not vendored.

The one third-party file is `html2canvas/html2canvas.min.js`. The PNG buttons use it to
rasterise a `<figure>`. It is served from this origin, the privacy note names it, and it
makes no network requests of its own. Do not replace it with a CDN `<script>` tag: that
would leak the fact that a user is using the tool and could carry sequence data in a URL.

If a future view genuinely needs another library, vendor it here as a static file with
its licence alongside.
