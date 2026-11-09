// SPDX-License-Identifier: MIT
//
// Per-figure PNG export. Run:
//   node test/figure_export_cases.mjs
//
// Two things are worth pinning down. The filename rule is pure and is what a user sees in
// their downloads folder, so it is checked directly. The attach/sync rule is the part that can
// go quietly wrong: a figure must offer a download exactly when it holds a plot, the button
// must not end up inside the image it produces, and attaching twice must not leave two
// buttons. Those need a DOM, so they run under jsdom and SKIP if it is absent, matching
// page_boot_cases.mjs.

import path from "node:path";
import fs from "node:fs";

const REPO = process.cwd();
let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

const MOD = path.join(REPO, "js/figure_export.js");
const { exportFileName, figureStamp } = await import(MOD);

console.log("Filenames");
check("a declared export name becomes the filename",
      exportFileName("batch_genome_track", "2026-11-09")
        === "mra_batch_genome_track_2026-11-09.png",
      exportFileName("batch_genome_track", "2026-11-09"));
check("spaces and punctuation collapse to single underscores",
      exportFileName("RT 103 \u2014 supply flow", "2026-11-09")
        === "mra_rt_103_supply_flow_2026-11-09.png",
      exportFileName("RT 103 \u2014 supply flow", "2026-11-09"));
check("a missing name still yields a usable filename",
      exportFileName("", "2026-11-09") === "mra_figure_2026-11-09.png",
      exportFileName("", "2026-11-09"));
check("a malformed date falls back to today, not to the literal string",
      /^mra_x_\d{4}-\d{2}-\d{2}\.png$/.test(exportFileName("x", "not-a-date")),
      exportFileName("x", "not-a-date"));
check("path separators cannot escape into the filename",
      !/[/\\]/.test(exportFileName("a/b\\c", "2026-11-09")),
      exportFileName("a/b\\c", "2026-11-09"));

console.log("Provenance stamp");
const stamp = figureStamp({ version: "v0.1", doi: "10.5281/zenodo.0" }, "2026-11-09");
check("the stamp names the atlas version", /v0\.1/.test(stamp), stamp);
check("the stamp carries the DOI when the release has one", /10\.5281/.test(stamp));
check("the stamp carries the retrieval date", /2026-11-09/.test(stamp));
check("the stamp carries the citation", /Nitido et al\. 2027/.test(stamp), stamp);
const cited = figureStamp({ version: "v0.1", paper: "Nitido et al. 2027" }, "2026-11-09");
check("a release without a DOI still names the DOI as unassigned",
      /DOI: \[not yet assigned\]/.test(cited) && /Nitido et al\. 2027/.test(cited), cited);
const bare = figureStamp(undefined, "2026-11-09");
check("a release-less stamp is still a sentence, not 'undefined'",
      !/undefined/.test(bare) && /Nitido et al\. 2027/.test(bare) &&
      /DOI: \[not yet assigned\]/.test(bare), bare);

console.log("The shipped pages declare a name for every figure");
// index.html is the landing page and carries the browse figures; browse.html is now a
// redirect stub with no figures of its own.
for (const page of ["index.html"]) {
  const html = fs.readFileSync(path.join(REPO, page), "utf8");
  // Match any figure whose class list contains "chart", not the exact attribute: the landing
  // page's track is class="chart plot-card" now that the figures are bordered cards. A
  // figure WITHOUT the chart class is a plain container -- index.html wraps the two
  // dashboard cards in one -- and offers no download, so it is not scanned here.
  const figs = (html.match(/<figure[^>]*class="[^"]*\bchart\b[^"]*"[^>]*>/g) || []);
  const named = figs.filter((f) => /data-export-name="[^"]+"/.test(f));
  check(`${page}: all ${figs.length} figures carry data-export-name`,
        figs.length > 0 && named.length === figs.length,
        `${named.length}/${figs.length}`);
  const names = [...html.matchAll(/data-export-name="([^"]+)"/g)].map((m) => m[1]);
  check(`${page}: export names are unique`,
        new Set(names).size === names.length, names.join(", "));
}

// --- attach / sync behaviour (needs a DOM) ----------------------------------------------
let JSDOM;
for (const spec of [process.env.JSDOM_PATH, "jsdom",
                    path.join(REPO, "node_modules/jsdom/lib/api.js"),
                    "/tmp/domtest/node_modules/jsdom/lib/api.js"]) {
  if (!spec) continue;
  try { ({ JSDOM } = await import(spec)); break; } catch { /* try the next location */ }
}
if (!JSDOM) {
  console.log("jsdom not found - SKIPPING attach/sync checks.");
  console.log("  enable with:  npm install jsdom   (or set JSDOM_PATH to its lib/api.js)");
} else {
  const { attachFigureDownloads, attachFigureDownload, figureIsDrawn, downloadFigure } =
    await import(MOD);

  const dom = new JSDOM(`<!doctype html><body>
    <figure class="chart" data-export-name="drawn_one"><svg></svg>
      <div class="track-slider" data-export-skip><input type="range"></div>
      <figcaption>caption</figcaption></figure>
    <figure class="chart" data-export-name="empty_one"><div></div>
      <figcaption></figcaption></figure>
    <figure class="chart" data-export-name="drawn_two"><div id="slot"><canvas></canvas></div>
      <figcaption></figcaption></figure>
  </body>`);
  const doc = dom.window.document;
  const [drawnOne, emptyOne, drawnTwo] = [...doc.querySelectorAll("figure.chart")];

  console.log("Which figures offer a download");
  check("a figure holding an <svg> counts as drawn", figureIsDrawn(drawnOne));
  check("a figure holding a <canvas> counts as drawn", figureIsDrawn(drawnTwo));
  check("an empty figure does not", !figureIsDrawn(emptyOne));

  const first = attachFigureDownloads({ root: doc });
  check("one button per drawn figure, none on the empty one",
        first.attached === 2 && first.drawn === 2,
        `attached=${first.attached} drawn=${first.drawn}`);
  check("the empty figure has no button", !emptyOne.querySelector(".figure-tools"));
  check("each drawn figure has exactly one button",
        drawnOne.querySelectorAll(".figure-download").length === 1 &&
        drawnTwo.querySelectorAll(".figure-download").length === 1);

  const again = attachFigureDownloads({ root: doc });
  check("syncing again adds nothing and removes nothing",
        again.attached === 0 && again.removed === 0,
        `attached=${again.attached} removed=${again.removed}`);
  check("still exactly one button per drawn figure",
        doc.querySelectorAll(".figure-download").length === 2,
        String(doc.querySelectorAll(".figure-download").length));

  console.log("The button stays out of its own image");
  check("the button row is marked data-export-skip",
        drawnOne.querySelector(".plot-tools").hasAttribute("data-export-skip"));
  check("the slider inside a track figure is marked data-export-skip",
        drawnOne.querySelector(".track-slider").hasAttribute("data-export-skip"));
  check("attaching to a figure that already has a button is a no-op",
        attachFigureDownload(drawnOne, {}) === null);

  console.log("A figure that is emptied again withdraws the offer");
  // The pages empty the inner chart slot, never the <figure> itself -- drawSite() does
  // el("codon-chart").innerHTML = "" at an unmodelled position. So the figure survives with a
  // stale button on it, and the next sync is what has to take that button away.
  doc.getElementById("slot").innerHTML = "";
  const third = attachFigureDownloads({ root: doc });
  check("the button is removed when its figure is cleared",
        third.removed === 1 && !drawnTwo.querySelector(".figure-tools"),
        `removed=${third.removed}`);
  check("the other figure keeps its button", !!drawnOne.querySelector(".figure-download"));

  console.log("Failure is reported, not swallowed");
  let msg = "";
  try { await downloadFigure(emptyOne, {}); } catch (e) { msg = e.message; }
  check("exporting an undrawn figure throws a readable error",
        /not been drawn/.test(msg), msg);
  // html2canvas is absent in this harness, so a drawn figure must fail on that, and must not
  // leave its provenance stamp behind in the page when it does.
  msg = "";
  try { await downloadFigure(drawnOne, {}); } catch (e) { msg = e.message; }
  check("a missing PNG library is reported by name", /PNG library/.test(msg), msg);
  check("a failed capture leaves no stamp behind in the page",
        !drawnOne.querySelector(".figure-stamp"));

  // What actually reaches html2canvas. Without a stub the checks above stop at the missing
  // library, so the contract that matters -- one figure in, controls left out, provenance
  // drawn in, nothing left behind -- would go untested.
  console.log("What the capture is handed");
  let seen = null;
  globalThis.html2canvas = async (node, opts) => {
    seen = {
      node,
      stampDuringCapture: node.querySelector(".figure-stamp")?.textContent || "",
      skipsToolRow: opts.ignoreElements(node.querySelector(".plot-tools")),
      skipsSlider: opts.ignoreElements(node.querySelector(".track-slider")),
      keepsPlot: opts.ignoreElements(node.querySelector("svg")),
      keepsCaption: opts.ignoreElements(node.querySelector("figcaption")),
      scale: opts.scale,
      background: opts.backgroundColor,
    };
    return { toBlob: (cb) => cb({ type: "image/png" }) };
  };
  // jsdom's URL has no object-URL support; the browser's does.
  let revoked = 0;
  const realURL = globalThis.URL;
  globalThis.URL = { createObjectURL: () => "blob:test", revokeObjectURL: () => { revoked++; } };
  let downloaded = "";
  const realCreate = doc.createElement.bind(doc);
  doc.createElement = (tag) => {
    const node = realCreate(tag);
    if (tag === "a") node.click = () => { downloaded = node.download; };
    return node;
  };
  try {
    await downloadFigure(drawnOne, { release: { version: "v0.1" }, dateISO: "2026-11-09" });
    check("html2canvas is handed the figure itself, not a page region",
          seen && seen.node === drawnOne);
    check("the button row is excluded from the capture", seen && seen.skipsToolRow === true);
    check("the slider is excluded from the capture", seen && seen.skipsSlider === true);
    check("the plot is kept", seen && seen.keepsPlot === false);
    check("the caption is kept, so the PNG says what it was drawn under",
          seen && seen.keepsCaption === false);
    check("the provenance stamp is present while the image is taken",
          seen && /atlas v0\.1/.test(seen.stampDuringCapture), seen && seen.stampDuringCapture);
    check("the stamp is gone from the page afterwards",
          !drawnOne.querySelector(".figure-stamp"));
    check("the capture is at 2x on white",
          seen && seen.scale === 2 && seen.background === "#ffffff");
    check("the file is offered under the figure's own name",
          downloaded === "mra_drawn_one_2026-11-09.png", downloaded);
    check("the object URL is revoked, so the blob is not held for the session",
          revoked === 1, String(revoked));
  } finally {
    delete globalThis.html2canvas;
    globalThis.URL = realURL;
    doc.createElement = realCreate;
  }
}

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL FIGURE EXPORT CASES PASSED");
