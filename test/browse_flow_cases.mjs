// SPDX-License-Identifier: MIT
//
// The precomputed LANL batch behind the landing page's supply-flow chart. Run:
//   node test/browse_flow_cases.mjs --version v0.1
//
// The risk this file exists for is silent staleness. The chart on the landing page is drawn
// from a file built once by scripts/build_lanl_flow.mjs; the chart on the batch tab is drawn
// from a batch run in the browser. They are supposed to be the same numbers, and nothing
// about shipping a JSON file enforces that -- a model change that is not followed by a
// rebuild would leave the landing page quietly showing last release's rates.
//
// So the expensive check is the real one: re-run the batch pipeline over the same FASTA and
// require the shipped file to agree, position for position. That costs ~30 s, which is why
// it only runs under --deep; the cheap structural checks always run.

import fs from "node:fs";
import path from "node:path";

const argv = process.argv.slice(2);
const version = argv.includes("--version") ? argv[argv.indexOf("--version") + 1] : "v0.1";
const deep = argv.includes("--deep");
const REPO = process.cwd();
const DATA = path.join(REPO, "data", version);
let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

const assetPath = path.join(DATA, "lanl_flow.json");
check("the precomputed batch is shipped", fs.existsSync(assetPath));
if (!fs.existsSync(assetPath)) {
  console.log("=".repeat(60));
  console.log("1 CHECK(S) FAILED"); process.exit(1);
}
const flow = JSON.parse(fs.readFileSync(assetPath, "utf8"));

console.log("The shipped file is what it claims to be");
check("it names the FASTA it was built from",
      flow.source === "lanl_batch_demo.fasta", flow.source);
check("it covers the whole compendium subset", flow.n_sequences === 182,
      String(flow.n_sequences));
check("row count matches the declared position count",
      flow.rows.length === flow.n_positions, `${flow.rows.length} vs ${flow.n_positions}`);
check("the residue order is the 20 amino acids plus stop",
      flow.aa_order.length === 21 && flow.aa_order.endsWith("*"), flow.aa_order);

console.log("Every row is well formed");
let badCounts = 0, badRates = 0, badIndex = 0, emptyCodons = 0, overCount = 0;
const proteins = new Set();
for (const [protein, position, nCovered, codons] of flow.rows) {
  proteins.add(protein);
  if (!codons.length) emptyCodons++;
  let summed = 0;
  for (const [codon, count, flat] of codons) {
    if (!/^[ACGT]{3}$/.test(codon) || !(count > 0)) badCounts++;
    summed += count;
    if (flat.length % 2 !== 0) { badIndex++; continue; }
    for (let i = 0; i < flat.length; i += 2) {
      if (!Number.isInteger(flat[i]) || flat[i] < 0 || flat[i] >= 21) badIndex++;
      // A rate is a per-replication probability: strictly positive (zeros are not stored)
      // and far below 1. A value at or above 1 would mean a certainty per cycle.
      if (!(flat[i + 1] > 0) || !(flat[i + 1] < 1)) badRates++;
    }
  }
  // Codon counts are samples covering this position, so they cannot exceed the batch size,
  // and they cannot exceed the position's own coverage.
  if (summed > nCovered || nCovered > flow.n_sequences) overCount++;
}
check("every codon is a plain ACGT triplet with a positive count", badCounts === 0,
      String(badCounts));
check("every destination index is inside the residue alphabet", badIndex === 0,
      String(badIndex));
check("every stored rate is in (0, 1)", badRates === 0, String(badRates));
check("no position ships with an empty codon list", emptyCodons === 0, String(emptyCodons));
check("codon counts never exceed the position's coverage", overCount === 0, String(overCount));
check("all nine proteins are represented", proteins.size >= 9,
      [...proteins].sort().join(", "));

console.log("The adapter rebuilds what the chart expects");
globalThis.fetch = async () => ({ ok: true, json: async () => flow });
const { lanlFlowAt, loadLanlFlow } = await import(`${REPO}/js/lanl_flow.js`);
const index = await loadLanlFlow(version);
check("the loader indexes every shipped position",
      index && index.byKey.size === flow.rows.length,
      index ? String(index.byKey.size) : "null");

const [p0, pos0] = [flow.rows[0][0], flow.rows[0][1]];
const got = lanlFlowAt(index, p0, pos0);
check(`${p0} ${pos0} resolves to a stat and a matrix`,
      !!(got && got.stat && got.matrix));
check("the stat carries codon counts sorted most-common first",
      got.stat.codon_counts.every((c, i, a) => i === 0 || a[i - 1][1] >= c[1]));
check("every codon in the stat has a destination vector",
      got.stat.codon_counts.every(([c]) => got.matrix.destinations.has("0|" + c)));
check("each codon maps to the residue the genetic code gives it",
      got.stat.codon_aa.size === got.stat.codon_counts.length &&
      [...got.stat.codon_aa.values()].every((a) => /^[A-Z*]$/.test(a)));
check("a position outside the compendium returns null",
      lanlFlowAt(index, "Gag", 99999) === null);
check("a missing index returns null rather than throwing",
      lanlFlowAt(null, p0, pos0) === null);

// The chart must accept what the adapter produces. This is the integration the adapter
// exists for, so it is checked against the real drawing function, not a stand-in.
let JSDOM;
for (const spec of [process.env.JSDOM_PATH, "jsdom",
                    "/tmp/domtest/node_modules/jsdom/lib/api.js"]) {
  if (!spec) continue;
  try { ({ JSDOM } = await import(spec)); break; } catch { /* next */ }
}
if (!JSDOM) {
  console.log("jsdom not found - SKIPPING the draw check.");
} else {
  const dom = new JSDOM("<!doctype html><body></body>");
  globalThis.window = dom.window;
  globalThis.document = dom.window.document;
  globalThis.Node = dom.window.Node;
  const { drawSupplyFlow, residueColour } = await import(`${REPO}/js/batch_view.js`);
  const box = dom.window.document.createElement("div");
  const cap = dom.window.document.createElement("div");
  const drew = drawSupplyFlow(box, cap, got.stat, got.matrix, { weight: "supply" });
  check("drawSupplyFlow renders from the precomputed data", drew !== false);
  check("it produced an SVG", !!box.querySelector("svg"));
  check("the caption is filled in", cap.textContent.trim().length > 20);

  // The shared colour key is the user-visible promise: the same residue must be the same
  // colour in the pie and in the flow chart, and two codons for one residue must agree.
  check("a residue's colour does not depend on which codon encodes it",
        residueColour("K") === residueColour("K") && residueColour("K") !== residueColour("E"));
  check("a stop is held apart from the residues",
        residueColour("*") !== residueColour("W"));

  // The user-visible promise is not that a palette function exists, but that the two figures
  // standing side by side agree. So render both for one position and compare the fills that
  // actually reach the SVG. A residue the flow chart paints must be painted the same colour
  // in the pie, or a reader carrying a colour from one figure to the other is misled.
  const { pieChart } = await import(`${REPO}/js/chart.js`);
  const aaOrder = flow.aa_order.split("");
  const reached = new Set();
  for (const [, , flat] of got.stat.codon_counts.map(([c]) =>
       [c, null, Array.from(got.matrix.destinations.get("0|" + c))])) {
    flat.forEach((v, i) => { if (v > 0 && aaOrder[i] !== "*") reached.add(aaOrder[i]); });
  }
  const pieBox = dom.window.document.createElement("div");
  pieChart(pieBox, [...reached].map((aa) => ({ label: aa, value: 1, colour: residueColour(aa) })),
           { title: "colour check", width: 420, height: 320 });
  // Only the data marks: wedges and bands are paths, rects and circles. Text and the
  // background carry fills too, and they are page chrome rather than a residue's colour.
  const fillsIn = (node) => new Set([...node.querySelectorAll("path, rect, circle")]
    .map((n) => (n.getAttribute("fill") || "").toLowerCase())
    .filter((f) => f.startsWith("#")));
  const pieFills = fillsIn(pieBox);
  const flowFills = fillsIn(box);
  const shared = [...reached].filter((aa) =>
    pieFills.has(residueColour(aa).toLowerCase()) &&
    flowFills.has(residueColour(aa).toLowerCase()));
  check("residues reached here are painted the same colour in both SVGs",
        reached.size > 0 && shared.length > 0,
        `${shared.length} of ${reached.size} residues found in both`);
  // The failure this guards against is the pie falling back to its positional palette, which
  // is what it used before: that palette assigns colour by slice order, so the same residue
  // changes colour whenever the ordering does. Every fill the pie emits must now be a
  // residue colour.
  const residuePalette = new Set(
    [...flow.aa_order].map((aa) => residueColour(aa).toLowerCase()));
  const strays = [...pieFills].filter((f) => !residuePalette.has(f));
  check("the pie uses no colour outside the residue palette", strays.length === 0,
        strays.join(", "));
  check("every residue the flow chart reaches has its residue colour in that SVG",
        [...reached].every((aa) => flowFills.has(residueColour(aa).toLowerCase())),
        [...reached].filter((aa) => !flowFills.has(residueColour(aa).toLowerCase())).join(","));

  // The bands between the two columns are the chart. Its label gutters used to be fixed at
  // 320 px left and 290 right, so below about 630 px they overran each other and the bands
  // were drawn with negative width -- two columns and nothing between them. Pin the span at
  // the width the dashboard uses and at widths narrow enough to have broken it.
  console.log("Bands are visible at every width the chart is drawn at");
  const bandSpan = (w) => {
    const b = dom.window.document.createElement("div");
    const c = dom.window.document.createElement("div");
    drawSupplyFlow(b, c, got.stat, got.matrix, { weight: "supply", width: w });
    let lo = Infinity, hi = -Infinity;
    for (const pth of b.querySelectorAll("path")) {
      for (const m of (pth.getAttribute("d") || "").matchAll(/[ML]\s*([-\d.]+),/g)) {
        const x = parseFloat(m[1]);
        if (x < lo) lo = x;
        if (x > hi) hi = x;
      }
    }
    return hi > lo ? hi - lo : 0;
  };
  const full = bandSpan(1140);
  check("at the dashboard's full width the bands span most of the figure",
        full > 400, `${full.toFixed(0)} px`);
  for (const w of [900, 700, 600, 480]) {
    const sp = bandSpan(w);
    check(`at ${w} px the bands are still drawn`, sp >= 150, `${sp.toFixed(0)} px`);
  }
  check("the gutters never invert, however narrow the figure",
        bandSpan(360) > 0 && bandSpan(300) > 0,
        `${bandSpan(360).toFixed(0)} / ${bandSpan(300).toFixed(0)} px`);

  // The dashboard panel itself: pie and flow chart in one row, and the panel still working
  // when no flow is available. The second case is the one that breaks quietly -- a position
  // the compendium does not cover, or an older release with no precomputed file.
  console.log("The position dashboard holds both figures");
  const { drawPositionPie } = await import(`${REPO}/js/track_view.js`);
  const rows = [
    { ref_aa: "K", alt_aa: "R", mu_aa: 3e-5, mu_aa_lo95: 1e-5, mu_aa_hi95: 9e-5,
      min_nt_changes: 1, n_paths: 1, synonymous: false },
    { ref_aa: "K", alt_aa: "E", mu_aa: 1e-5, mu_aa_lo95: 4e-6, mu_aa_hi95: 3e-5,
      min_nt_changes: 1, n_paths: 1, synonymous: false },
    { ref_aa: "K", alt_aa: "T", mu_aa: 2e-6, mu_aa_lo95: 8e-7, mu_aa_hi95: 6e-6,
      min_nt_changes: 2, n_paths: 2, synonymous: false },
  ];
  const panel = dom.window.document.createElement("div");
  const panelCap = dom.window.document.createElement("div");
  const drewPanel = drawPositionPie(panel, panelCap, rows, "RT 103",
                                    { custom: false, flow: got });
  check("the panel draws with a flow source", drewPanel === true);
  check("it holds both figures, the flow chart stacked under the pie",
        panel.querySelectorAll(".panel-row svg").length >= 2,
        String(panel.querySelectorAll(".panel-row svg").length));
  // The page draws the panel and then syncs the PNG buttons; do both, or the download check
  // below tests the harness rather than the app.
  const { attachFigureDownloads } = await import(`${REPO}/js/figure_export.js`);
  attachFigureDownloads({ root: panel });
  const cards = [...panel.querySelectorAll(".panel-row > .plot-card")];
  const vb = (c) => (c.querySelector("svg").getAttribute("viewBox") || "").split(" ");
  check("both plots sit in their own bordered card", cards.length === 2,
        cards.map((c) => c.getAttribute("data-export-name")).join(", "));
  check("the two cards are the same size",
        new Set(cards.map((c) => vb(c)[2] + "x" + vb(c)[3])).size === 1,
        cards.map((c) => vb(c)[2] + "x" + vb(c)[3]).join(" vs "));
  check("together they span the genome track's width",
        cards.reduce((t, c) => t + Number(vb(c)[2]), 0) + 16 === 1140,
        String(cards.reduce((t, c) => t + Number(vb(c)[2]), 0) + 16));
  check("each card carries its own PNG download",
        cards.every((c) => c.querySelector(".plot-tools .figure-download")));
  check("each card carries its own title inside the plot",
        cards.every((c) => /\S/.test((c.querySelector("svg text") || {}).textContent || "")),
        cards.map((c) => (c.querySelector("svg text") || {}).textContent).join(" | "));
  check("no emergence readout table is left in the panel",
        !panel.querySelector("table.track-stats"));
  // The explanation moved into each card's Details dialog; the shared caption under the
  // pair would otherwise repeat what the two plot titles already say.
  check("the shared caption is not another paragraph of prose",
        panelCap.textContent.trim() === "", panelCap.textContent.trim().slice(0, 40));
  check("each card explains itself behind a Details button",
        cards.every((c) => c.querySelector(".plot-tools .details-button")));
  const dialogs = [...panel.querySelectorAll("dialog.details-dialog")];
  check("the explanation is kept, not deleted", dialogs.length === 2 &&
        dialogs.every((x) => x.textContent.length > 120),
        dialogs.map((x) => x.textContent.length).join(", "));
  check("no details dialog still promises a probability curve",
        !dialogs.some((x) => /probabilit/i.test(x.textContent)));

  const solo = dom.window.document.createElement("div");
  const soloCap = dom.window.document.createElement("div");
  const drewSolo = drawPositionPie(solo, soloCap, rows, "RT 103", { custom: true });
  check("the panel still draws with no flow source", drewSolo === true);
  check("without a flow chart it shows one figure, not an empty slot",
        solo.querySelectorAll(".panel-row svg").length === 1,
        String(solo.querySelectorAll(".panel-row svg").length));
  check("and on its own the pie keeps a visible one-line caption",
        /\S/.test(soloCap.textContent));
}

if (deep) {
  console.log("The shipped file still matches a fresh batch (deep)");
  const { MraEngine } = await import(`${REPO}/js/mra_engine.js`);
  const { Analysis } = await import(`${REPO}/js/analysis.js`);
  const { parseFasta } = await import(`${REPO}/js/align.js`);
  const { buildAaIndex, aaPositionStats } = await import(`${REPO}/js/batch_core.js`);
  const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));
  const assets = { bundle: J("mra_web_model.json"), reference: J("nl43_reference.json"),
                   intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") };
  void new MraEngine(assets);
  const records = parseFasta(fs.readFileSync(path.join(DATA, "lanl_batch_demo.fasta"), "utf8"));
  const analysis = new Analysis(assets);
  await analysis.start();
  const batch = await analysis.batch(records, () => {});
  const stats = aaPositionStats(batch.aaMatrix, buildAaIndex(assets.reference.proteins),
                                { field: "nonsyn" });
  const live = new Map(stats.map((s) => [s.protein + "|" + s.position, s]));
  check("the fresh batch has the same number of covered positions",
        stats.length === flow.rows.length, `${stats.length} vs ${flow.rows.length}`);
  let mismatch = 0, rateOff = 0, compared = 0;
  for (const [protein, position, nCovered, codons] of flow.rows) {
    const s = live.get(protein + "|" + position);
    if (!s || s.n_covered !== nCovered) { mismatch++; continue; }
    const liveCounts = new Map(s.codon_counts);
    for (const [codon, count, flat] of codons) {
      if (liveCounts.get(codon) !== count) { mismatch++; continue; }
      const d = batch.aaMatrix.destinations.get(s.position_index + "|" + codon);
      for (let i = 0; i < flat.length; i += 2) {
        compared++;
        // Shipped at 3 significant figures, so agreement is relative, not exact.
        if (Math.abs(flat[i + 1] - d[flat[i]]) > 5e-3 * d[flat[i]]) rateOff++;
      }
    }
  }
  check("every position, codon and count matches the fresh batch", mismatch === 0,
        String(mismatch));
  check("every shipped rate matches to 3 significant figures", rateOff === 0,
        `${rateOff} off of ${compared}`);
} else {
  console.log("Deep rebuild comparison SKIPPED - pass --deep to run it (~30 s).");
}

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL BROWSE FLOW CASES PASSED");
