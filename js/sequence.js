// SPDX-License-Identifier: MIT
//
// sequence.js -- submit a sequence, get a custom atlas.
//
// The scientific framing is load-bearing here, not decoration. This view must never
// let a reader mistake a modelled rate for a measured one, so:
//   - the privacy note and the standing caveat sit at the end of the page, after the
//     results, which is where a reader meets them with the numbers already in view;
//   - the two tables under the charts match Browse: amino-acid changes, then codon paths,
//     with a sample column naming which submitted sequence each row came from;
//   - a sequence that aligns entirely outside the modelled region gets a refusal that
//     names the covered region, not an empty table.
//
// Up to MAX_DETAIL_SEQUENCES submissions share one genome track and one grid of pies, one
// pie per sequence at the marker. Checkboxes in the control bar choose which lines and pies
// are drawn; the table under the grid still lists every sequence. The two genome-wide
// tables match Browse, with a sample column, and each downloads in full.

import { bootstrapPage, privacyPanel, note, caveat, showError, exportProvenance,
         CUSTOM_CAVEAT, NEF_NOTE } from "./page.js";
import { loadSmallAssets } from "./data.js";
import { Analysis, progressRow, MAX_DETAIL_SEQUENCES } from "./analysis.js";
import { parseFasta, alignToReference } from "./align.js";
import { renderTable, formatRate } from "./table.js";
import { MraEngine, GENETIC_CODE } from "./mra_engine.js";
import { perSiteSupply, groupByPosition, buildFrameIndex, geneSegments, geneNtWindow,
         targetFor, smoothSeries, oneStepSlices } from "./explore_core.js";
import { drawSupplyTrack, bindAxisControls, CLASS_LABEL } from "./track_view.js";
import { pieChart, LAYOUT, bindZoomRange } from "./chart.js";
import { residueColour } from "./batch_view.js";
import { downloadFigure } from "./figure_export.js";
import { setDetails } from "./figure_details.js";

// Distinct from the dashed NL4-3 grey (#8a94a0). One colour per submitted sequence, in
// submission order, and stable when a line is toggled off.
const OVERLAY_COLOURS = [
  "#2f4b7c", "#c2561a", "#1f7a3f", "#6b3fa0", "#a4262c",
  "#0b6e8a", "#8a6d1f", "#3d6b4f", "#7a3e6d", "#3f6fbf",
];

const page = await bootstrapPage("sequence.html");
const region = page.release.modelled_region;

document.getElementById("privacy-slot").appendChild(privacyPanel());
document.getElementById("caveat-slot").appendChild(note(CUSTOM_CAVEAT));
document.getElementById("batch-note").textContent =
  `Up to ${MAX_DETAIL_SEQUENCES} sequences per batch.`;

const statusSlot = document.getElementById("run-status");
let assets, analysis, engine;
try {
  // The sequence tool runs on the small assets alone -- about 360 KB gzipped. The
  // multi-megabyte browse tables are never fetched on this page.
  assets = await loadSmallAssets(page.version);
  engine = new MraEngine({ bundle: assets.bundle, reference: assets.reference,
                           intervals: assets.intervals, hxb2: assets.hxb2 });
  analysis = new Analysis(assets);
  const inWorker = await analysis.start();
  statusSlot.appendChild(caveat(
    inWorker
      ? "Ready. Alignment and rate calculation run in a background thread, so the page stays responsive."
      : "Ready. This browser does not support module workers, so calculation runs on the main " +
        "thread and the page may pause briefly while a long sequence is processed."));
} catch (err) {
  showError(statusSlot, "Could not load the model: " + err.message);
  throw err;
}

// --- input -----------------------------------------------------------------------
const fasta = document.getElementById("fasta");

document.getElementById("load-file").addEventListener("click",
  () => document.getElementById("file-input").click());
document.getElementById("file-input").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (file) fasta.value = await file.text();        // read locally, never uploaded
});
for (const ev of ["dragover", "drop"]) {
  fasta.addEventListener(ev, (e) => { e.preventDefault(); });
}
fasta.addEventListener("drop", async (e) => {
  const file = e.dataTransfer.files[0];
  if (file) fasta.value = await file.text();
});
document.getElementById("load-example").addEventListener("click", async () => {
  // Shipped example: five reference genomes in HIV-references.fasta. Fetched from this
  // origin and written into the box; nothing is uploaded.
  try {
    const resp = await fetch("HIV-references.fasta");
    if (!resp.ok) throw new Error("HIV-references.fasta could not be loaded");
    fasta.value = await resp.text();
  } catch (err) {
    showError(document.getElementById("detail-slot"),
      "Could not load the example sequences. " + err.message);
  }
});
document.getElementById("clear").addEventListener("click", () => {
  fasta.value = "";
  document.getElementById("detail-slot").innerHTML = "";
  document.getElementById("results-heading").hidden = true;
  sessionStorage.removeItem("mra.sequences");
  analysis.dispose();
  analysis = new Analysis(assets);
  analysis.start();
});

// --- run -------------------------------------------------------------------------
document.getElementById("run").addEventListener("click", async () => {
  const text = fasta.value.trim();
  const detailSlot = document.getElementById("detail-slot");
  detailSlot.innerHTML = "";
  document.getElementById("results-heading").hidden = true;
  if (!text) { showError(detailSlot, "Paste a sequence first."); return; }

  const records = parseFasta(text).filter((r) => r.sequence.length >= 30);
  if (!records.length) {
    showError(detailSlot, "No sequence of at least 30 bases was found in that input.");
    return;
  }

  // Carried for the comparison view, in this tab only. sessionStorage is per-tab and is
  // discarded when the tab closes; nothing is written to disk by the site and nothing is
  // transmitted. "Clear everything" removes it immediately.
  try { sessionStorage.setItem("mra.sequences", text); } catch { /* quota: not essential */ }

  const bar = progressRow(document.getElementById("progress-slot"));
  const t0 = performance.now();
  let out;
  try {
    out = await analysis.run(records, (p) => bar.update(p));
  } catch (err) {
    showError(detailSlot, "Analysis failed: " + err.message);
    return;
  }
  const ms = Math.round(performance.now() - t0);
  bar.finish(`Done: ${out.summaries.length} sequence(s) in ${ms} ms.`);
  document.getElementById("results-heading").hidden = false;

  const rows = summariseRows(out, records);
  const analysed = rows.filter((r) => r.status === "analysed");
  if (!analysed.length) {
    if (out.truncated) detailSlot.appendChild(truncationNote());
    appendSummary(detailSlot, rows);
    return;
  }

  const loading = caveat("Building the tracks\u2026");
  detailSlot.appendChild(loading);
  const items = [];
  try {
    for (const r of analysed) {
      const detail = await analysis.detail(r.index);
      const item = prepareItem(r, detail, out.geometry[r.index], OVERLAY_COLOURS[items.length]);
      if (item.sites.length) items.push(item);
    }
  } catch (err) {
    showError(detailSlot, "Could not retrieve results: " + err.message);
    return;
  }
  loading.remove();
  if (!items.length) {
    if (out.truncated) detailSlot.appendChild(truncationNote());
    appendSummary(detailSlot, rows);
    detailSlot.appendChild(caveat("No modelled nucleotide position was covered, so there is nothing to plot."));
    return;
  }
  renderDashboard(detailSlot, items, rows, out.truncated);
});

function truncationNote() {
  return note(
    `Only the first ${MAX_DETAIL_SEQUENCES} sequences were analysed. The cap is a named ` +
    `constant (<code>MAX_DETAIL_SEQUENCES</code> in <code>js/align.js</code>).`);
}

function summariseRows(out, records) {
  return out.summaries.map((s, i) => ({
    index: i,
    name: s.name,
    length: records[i].sequence.length,
    status: s.unalignable ? "not placed"
      : s.codonsAnalysed === 0 ? "outside modelled region" : "analysed",
    identity: s.unalignable ? null : Math.round(s.identity * 1000) / 10,
    coverage: s.unalignable ? null : (out.geometry[i] ? Math.round(out.geometry[i].coverageOfQuery * 1000) / 10 : null),
    ref_span: s.unalignable || !out.geometry[i] ? "\u2014"
      : `${out.geometry[i].refStart.toLocaleString()}\u2013${out.geometry[i].refEnd.toLocaleString()}`,
    orientation: s.unalignable || !out.geometry[i] ? "\u2014"
      : (out.geometry[i].orientation === "reverse_complement" ? "reverse complement" : "forward"),
    codons: s.codonsAnalysed || 0,
    aa_rows: s.aaRows || 0,
    dropped: s.droppedOutsideRegion || 0,
    fallback: s.fallbackPositions || 0,
    summary: s,
  }));
}

const SUMMARY_COLUMNS = [
  { key: "name", label: "Sequence" },
  { key: "length", label: "bp", type: "num" },
  { key: "status", label: "Status" },
  { key: "identity", label: "% identity", type: "num" },
  { key: "coverage", label: "% of query placed", type: "num" },
  { key: "ref_span", label: "NL4-3 span" },
  { key: "orientation", label: "Orientation" },
  { key: "codons", label: "Codons", type: "num" },
  { key: "dropped", label: "Dropped > 8,213", type: "num",
    title: "Positions that aligned outside the modelled region and were removed" },
  { key: "fallback", label: "Spline fallback", type: "num",
    title: "Positions where the genome-mean spline was substituted because the site could not be aligned" },
];

function appendSummary(slot, rows) {
  const heading = document.createElement("h3");
  heading.textContent = "Sequences";
  slot.appendChild(heading);
  const tableSlot = document.createElement("div");
  slot.appendChild(tableSlot);
  renderTable(tableSlot, { columns: SUMMARY_COLUMNS, rows, filterPlaceholder: false });
  for (const r of rows.filter((x) => x.status !== "analysed")) {
    const div = document.createElement("div");
    div.className = "error-box";
    div.innerHTML = r.status === "not placed"
      ? `<strong>${escapeHtml(r.name)} was not placed on the reference.</strong> ` +
        escapeHtml(r.summary.reason)
      : `<strong>${escapeHtml(r.name)} lies entirely outside the modelled region.</strong> ` +
        `The atlas covers NL4-3 reference positions ${region.first_position.toLocaleString()}` +
        `&ndash;${region.last_position.toLocaleString()}; this sequence aligned beyond it, so ` +
        `no position-resolved rates can be given. ${NEF_NOTE}`;
    slot.appendChild(div);
  }
}

function prepareItem(row, detail, geom, colour) {
  const perSite = perSiteSupply(groupByPosition(detail.nt), "rate");
  const sites = [...perSite.keys()].sort((a, b) => a - b);
  const refSeq = assets.reference.sequence;
  const baseAt = new Map();
  for (const r of detail.nt) if (!baseAt.has(r.ref_position)) baseAt.set(r.ref_position, r.ref_nucleotide);
  const variants = [...baseAt.entries()]
    .filter(([pos, base]) => base !== refSeq[pos - 1]).map(([pos]) => pos).sort((a, b) => a - b);
  const aaByPosition = new Map();
  for (const r of detail.aa) {
    const k = r.protein + "|" + r.position;
    let list = aaByPosition.get(k);
    if (!list) aaByPosition.set(k, (list = []));
    list.push(r);
  }
  const proteinSpans = new Map();
  for (const r of detail.aa) {
    const cur = proteinSpans.get(r.protein);
    if (!cur) proteinSpans.set(r.protein, { min: r.position, max: r.position });
    else { cur.min = Math.min(cur.min, r.position); cur.max = Math.max(cur.max, r.position); }
  }
  return {
    row, detail, geom, colour, perSite, sites, baseAt, variants, aaByPosition, proteinSpans,
    indels: indelRatesFor(sites),
    frameAt: buildFrameIndex(assets.reference.proteins, aaByPosition),
    targets: engine.proteinTargets(proteinSpans),
    variantIdx: 0,
  };
}

/** Nearest nucleotide that is inside a codon this sequence's atlas can chart. A partial
 *  amplicon covers bases at its edge that belong to a codon the alignment did not complete,
 *  and those bases have a rate but no pie. */
function nearestCodon(item, nt) {
  if (item.frameAt.has(nt)) return nt;
  let best = null, bestD = Infinity;
  for (const s of item.frameAt.keys()) {
    const d = Math.abs(s - nt);
    if (d < bestD) { best = s; bestD = d; }
  }
  return best;
}

/**
 * One genome track for the whole submission, then a grid of pies, one per sequence, all
 * at the marker. Checkboxes choose which lines are drawn and do not hide a pie. The
 * smoothing and coordinate controls sit between the Results heading and the plot, and
 * the position slider sits in the track between the gene boxes and the rate curve,
 * matching Browse.
 */
function renderDashboard(slot, items, allRows, truncated) {
  slot.innerHTML = "";
  const SEGMENTS = geneSegments(assets.reference.proteins);
  const shown = new Set(items.map((_, i) => i));
  let active = 0;
  let selected = null;
  // Open on RT K65 (middle base 1954) when the sequence covers it.
  const prefer = 1954;
  let markerNt = nearestCodon(items[0], items[0].frameAt.has(prefer) ? prefer : (items[0].variants[0] || prefer));
  if (markerNt === null) {
    const donor = items.find((it) => it.frameAt.size);
    if (donor) markerNt = nearestCodon(donor, prefer);
  }
  // One amino acid for every pie. A new position starts on the largest one-step change.
  let picked = null;
  let pieProtein = null;
  let pieGeneAa = null;

  const fig = document.createElement("figure");
  fig.className = "chart plot-card";
  fig.setAttribute("data-export-name", "genome_track");
  const stack = document.createElement("div");
  stack.className = "genome-stack";
  const trackSlot = document.createElement("div");
  const sliderWrap = document.createElement("div");
  sliderWrap.className = "track-slider";
  sliderWrap.setAttribute("data-export-skip", "");
  const elSlider = document.createElement("input");
  elSlider.type = "range";
  elSlider.step = "1";
  elSlider.min = String(region.first_position);
  elSlider.max = String(region.last_position);
  elSlider.setAttribute("aria-label", "Genome position (NL4-3 nucleotide)");
  sliderWrap.appendChild(elSlider);
  stack.append(trackSlot, sliderWrap);
  const figCap = document.createElement("figcaption");
  const toggles = document.createElement("div");
  toggles.className = "overlay-row";
  toggles.setAttribute("data-export-skip", "");
  const toggleLabel = document.createElement("span");
  toggleLabel.className = "overlay-label";
  toggleLabel.textContent = "On the track";
  toggles.appendChild(toggleLabel);
  items.forEach((item, i) => {
    const lab = document.createElement("label");
    const cb = document.createElement("input");
    cb.type = "checkbox";
    cb.checked = true;
    cb.addEventListener("change", () => {
      if (cb.checked) shown.add(i); else shown.delete(i);
      drawTrack();
      syncPieGrid();
    });
    const swatch = document.createElement("i");
    swatch.className = "swatch";
    swatch.style.background = item.colour;
    const name = document.createElement("span");
    name.textContent = item.row.name;
    lab.append(cb, swatch, name);
    toggles.appendChild(lab);
  });
  const nl43 = document.createElement("span");
  nl43.className = "track-fixed";
  const nlSwatch = document.createElement("i");
  nlSwatch.className = "swatch";
  nlSwatch.style.background = "#8a94a0";
  nl43.append(nlSwatch, document.createTextNode("NL4-3"));
  toggles.appendChild(nl43);
  toggles.classList.add("track-names");
  fig.append(stack, figCap);

  const ctrl = document.createElement("div");
  ctrl.className = "control-bar";
  ctrl.setAttribute("data-export-skip", "");
  ctrl.innerHTML =
    `<div class="control-main">` +
    `<div><label>Supply shown</label><select class="sq-class">` +
    `<option value="total">Total substitution supply</option>` +
    `<option value="nonga">non-G&rarr;A only</option>` +
    `<option value="ga">G&rarr;A only</option>` +
    `<option value="insertion">Insertions</option>` +
    `<option value="deletion">Deletions</option></select></div>` +
    `<div><label>Jump to</label><select class="sq-sys">` +
    `<option value="protein">Protein position</option>` +
    `<option value="nt">NL4-3 nucleotide</option>` +
    `<option value="variant">Difference from NL4-3</option></select></div>` +
    `<div><label>Protein or subunit</label><select class="sq-target"></select></div>` +
    `<div><label>Position</label><input type="number" class="sq-pos" min="1" size="8"></div>` +
    `<div><label>Reading frame here</label><select class="sq-frame"></select></div>` +
    `<div class="zoom-fields"><label>Zoom</label>` +
    `<div class="zoom-range"><span class="zoom-range-fill sq-zfill"></span>` +
    `<input type="range" class="sq-z0" step="1" aria-label="Zoom from">` +
    `<input type="range" class="sq-z1" step="1" aria-label="Zoom to"></div>` +
    `<output class="zoom-readout sq-zread"></output>` +
    `<button type="button" class="secondary sq-zoom-reset">Full range</button></div>` +
    `<div class="bar-actions"><button class="sq-go">Go</button></div>` +
    `</div>`;
  const foot = document.createElement("div");
  foot.className = "control-foot";
  const axisBar = document.createElement("div");
  axisBar.className = "axis-controls";
  axisBar.innerHTML =
    `<span class="overlay-label">Y-axis</span>` +
    `<button type="button" class="toggle y-autoscale is-on" aria-pressed="true">Auto-scale</button>` +
    `<div class="seg" role="group" aria-label="Y-axis scale">` +
    `<button type="button" class="y-log is-on" aria-pressed="true">Log</button>` +
    `<button type="button" class="y-linear" aria-pressed="false">Linear</button></div>`;
  const smooth = document.createElement("div");
  smooth.className = "smooth-window";
  smooth.innerHTML =
    `<label>Smoothing window (nt)</label>` +
    `<input type="range" class="sq-win" min="25" max="501" step="1" value="250">` +
    `<output class="sq-win-out">250</output>`;
  foot.append(toggles, axisBar, smooth);
  ctrl.append(foot);
  const axis = bindAxisControls(ctrl, () => drawTrack());

  const head = document.createElement("div");
  const pieBoard = document.createElement("section");
  pieBoard.className = "pie-board plot-card";
  pieBoard.id = "pie-board";
  const pieHead = document.createElement("div");
  pieHead.className = "pie-board-head";
  const pieTitle = document.createElement("h3");
  pieTitle.className = "plot-title";
  const aaBar = document.createElement("div");
  aaBar.className = "pie-aa-bar";
  aaBar.setAttribute("data-export-skip", "");
  const aaLabel = document.createElement("label");
  aaLabel.textContent = "Amino acid change ";
  const aaSelect = document.createElement("select");
  aaSelect.setAttribute("aria-label", "Amino acid change");
  aaLabel.appendChild(aaSelect);
  aaBar.appendChild(aaLabel);
  pieHead.append(pieTitle, aaBar);
  const pieGrid = document.createElement("div");
  pieGrid.className = "pie-grid";
  const pieGridObserver = new ResizeObserver(() => sizePieGrid());
  pieGridObserver.observe(pieGrid);
  const selSlot = document.createElement("div");
  selSlot.className = "pie-selection";
  pieBoard.append(pieHead, pieGrid, selSlot);
  const lower = document.createElement("div");
  slot.append(ctrl, fig, head, pieBoard, lower);

  const elWin = ctrl.querySelector(".sq-win"), elWinOut = ctrl.querySelector(".sq-win-out");
  const elClass = ctrl.querySelector(".sq-class"), elSys = ctrl.querySelector(".sq-sys");
  const elTarget = ctrl.querySelector(".sq-target"), elPos = ctrl.querySelector(".sq-pos");
  const elFrame = ctrl.querySelector(".sq-frame");
  const elZ0 = ctrl.querySelector(".sq-z0"), elZ1 = ctrl.querySelector(".sq-z1");
  const axis0 = region.first_position, axis1 = region.last_position;
  let zoom0 = axis0, zoom1 = axis1;
  elZ0.min = elZ1.min = String(axis0);
  elZ0.max = elZ1.max = String(axis1);
  elZ0.value = String(axis0);
  elZ1.value = String(axis1);
  const zoomCtl = bindZoomRange(elZ0, elZ1, {
    fill: ctrl.querySelector(".sq-zfill"),
    readout: ctrl.querySelector(".sq-zread"),
    onChange: (a, b) => setZoom(a, b),
  });

  function markSelected() {
    toggles.querySelectorAll("label").forEach((lab, i) => {
      lab.classList.toggle("is-selected", i === active);
    });
  }
  function fillTargets() {
    elTarget.innerHTML = "";
    for (const t of items[active].targets) {
      const o = document.createElement("option");
      o.value = t.key;
      o.textContent = t.label;
      elTarget.appendChild(o);
    }
  }

  function placeTrackSlider() {
    const svg = trackSlot.querySelector("svg");
    if (!svg || !svg.dataset.sliderTop) return;
    const H = Number(svg.dataset.sliderSvgHeight);
    const W = Number(svg.dataset.sliderSvgWidth);
    sliderWrap.style.top = (100 * Number(svg.dataset.sliderTop) / H).toFixed(3) + "%";
    sliderWrap.style.left = (100 * Number(svg.dataset.sliderLeft) / W).toFixed(3) + "%";
    sliderWrap.style.width = (100 * Number(svg.dataset.sliderWidth) / W).toFixed(3) + "%";
  }

  function setZoom(a, b) {
    let lo = Math.round(Math.min(a, b));
    let hi = Math.round(Math.max(a, b));
    lo = Math.max(axis0, Math.min(lo, axis1));
    hi = Math.max(axis0, Math.min(hi, axis1));
    if (!(hi > lo)) return;
    if (lo === zoom0 && hi === zoom1) {
      zoomCtl.set(lo, hi);
      return;
    }
    zoom0 = lo;
    zoom1 = hi;
    zoomCtl.set(lo, hi);
    drawTrack();
  }
  function resetZoom() {
    zoom0 = axis0;
    zoom1 = axis1;
    zoomCtl.set(axis0, axis1);
    drawTrack();
  }

  function primeZoom(lo, hi) {
    lo = Math.max(axis0, Math.min(Math.round(lo), axis1));
    hi = Math.max(axis0, Math.min(Math.round(hi), axis1));
    if (!(hi > lo)) return;
    zoom0 = lo;
    zoom1 = hi;
    zoomCtl.set(lo, hi);
  }

  function zoomToGene(protein, geneAa, target) {
    const win = geneNtWindow(assets.reference.proteins, protein, {
      offset: target ? target.offset || 0 : 0,
      maxOwn: target && target.maxOwn ? target.maxOwn : null,
      geneAa,
    });
    if (win) primeZoom(win.start, win.end);
  }

  function zoomToNt(item, nt) {
    const frames = item.frameAt.get(nt) || [];
    const keep = selected && frames.find((f) => f.protein === selected.protein);
    const frame = keep || frames[0];
    if (!frame) return;
    zoomToGene(frame.protein, frame.geneAa, targetFor(item.targets, frame.protein, frame.geneAa));
  }

  function drawTrack() {
    const item = items[active];
    const win = Number(elWin.value);
    elWinOut.textContent = String(win);
    const cls = elClass.value;
    const indel = cls === "insertion" || cls === "deletion";
    const series = [];
    items.forEach((it, i) => {
      if (!shown.has(i)) return;
      series.push({
        label: shortName(it.row.name),
        name: it.row.name,
        colour: it.colour,
        points: smoothSeries(it.sites, indel ? it.indels : it.perSite, cls, win),
        width: i === active ? 1.7 : 1.15,
      });
    });
    const ref = referenceAtlas();
    const trackText = drawSupplyTrack(trackSlot, {
      perSite: indel ? item.indels : item.perSite, sites: item.sites, series, win, cls,
      segments: SEGMENTS,
      first: zoom0, last: zoom1,
      marker: markerNt,
      markerLabel: selected ? engine.hxb2Label(selected.protein, selected.geneAa) : String(markerNt),
      variants: item.variants,
      title: `Smoothed ${CLASS_LABEL[cls]} along the submitted sequences, with NL4-3 overlaid`,
      source: indel
        ? "the NL4-3 rate at each covered position; an insertion or deletion has no sequence context"
        : "computed from the codons and 3\u2032 context in the sequences you submitted, not from NL4-3",
      refPerSite: indel ? ref.indels : ref.perSite,
      subjectName: shortName(item.row.name),
      refBase: item.baseAt.get(markerNt),
      onSelect: (pos) => selectNt(pos),
      onZoom: (a, b) => setZoom(a, b),
      onZoomReset: resetZoom,
      yScale: axis.yScale,
      autoScale: axis.autoScale,
      // Room between the gene boxes and the rate curve for the position slider.
      sliderGap: 28,
    });
    elSlider.min = String(zoom0);
    elSlider.max = String(zoom1);
    placeTrackSlider();
    const where = selected ? ` (${engine.hxb2Label(selected.protein, selected.geneAa)})` : "";
    const zoomed = zoom0 !== axis0 || zoom1 !== axis1;
    figCap.innerHTML =
      `Genome-wide ${CLASS_LABEL[cls]}, ${win} nt mean. Marker at ` +
      `${markerNt.toLocaleString()}${where}. Red ticks are where ` +
      `<strong>${escapeHtml(shortName(item.row.name))}</strong> differs from NL4-3. ` +
      (indel
        ? `The dashed grey line is NL4-3, and each submitted sequence follows it wherever that sequence covers the site. `
        : `The dashed grey line is NL4-3. `) +
      (zoomed
        ? ` Showing positions ${zoom0.toLocaleString()}\u2013${zoom1.toLocaleString()}.`
        : ``) +
      ` Drag across the curve or the Zoom slider to zoom, click a gene to zoom to it, or double-click the plot for the full range.`;
    setDetails(fig, trackText, { title: "The genome supply track" });
  }

  const PIE_W = 170, PIE_H = 140, PIE_R = 52;

  function rowsAt(item) {
    if (pieProtein === null) return [];
    return item.aaByPosition.get(pieProtein + "|" + pieGeneAa) || [];
  }
  function oneStepAt(item) {
    return oneStepSlices(rowsAt(item));
  }
  function sliceKey(row) {
    return row.alt_aa + (row.synonymous ? "*" : "");
  }
  function pickKey(pick) {
    return pick && pick.altAa ? "aa:" + pick.altAa : "";
  }
  function pickFromKey(key) {
    if (key && key.startsWith("aa:")) return { kind: "aa", altAa: key.slice(3) };
    return null;
  }
  function sameChange(row, pick) {
    if (!pick || !row) return false;
    // The destination letter is one row. When this sequence already carries that residue,
    // the row is the synonymous change, and its rate is the one to report.
    return row.alt_aa === pick.altAa;
  }
  /** Highest-rate path of two or three nucleotide changes to this amino acid. */
  function bestFarPath(item, altAa) {
    if (!item._farPaths) {
      const m = new Map();
      for (const p of item.detail.paths) {
        if (p.n_nt_changes < 2 || p.n_nt_changes > 3) continue;
        const k = p.protein + "|" + p.position + "|" + p.alt_aa;
        const cur = m.get(k);
        if (!cur || p.mu_path > cur.mu_path) m.set(k, p);
      }
      item._farPaths = m;
    }
    if (pieProtein === null) return null;
    return item._farPaths.get(pieProtein + "|" + pieGeneAa + "|" + altAa) || null;
  }
  function changeCatalog() {
    const byKey = new Map();
    for (const item of items) {
      for (const r of rowsAt(item)) {
        if (r.synonymous) continue;
        const key = "aa:" + r.alt_aa;
        let cur = byKey.get(key);
        if (!cur) {
          cur = { key, altAa: r.alt_aa, one: false, score: -1, minSteps: r.min_nt_changes };
          byKey.set(key, cur);
        }
        if (r.min_nt_changes === 1) cur.one = true;
        if (r.mu_aa > cur.score) cur.score = r.mu_aa;
        cur.minSteps = Math.min(cur.minSteps, r.min_nt_changes);
      }
    }
    const list = [...byKey.values()];
    list.sort((a, b) => (b.one - a.one) || (b.score - a.score) || a.key.localeCompare(b.key));
    for (const c of list) {
      const name = c.altAa === "*" ? "stop" : c.altAa;
      c.label = c.one ? name : `${name} (${c.minSteps}-step)`;
    }
    return list;
  }
  function ensurePick() {
    const catalog = changeCatalog();
    if (!picked || !catalog.some((c) => c.key === pickKey(picked))) {
      const first = catalog.find((c) => c.one) || catalog[0];
      picked = first ? pickFromKey(first.key) : null;
    }
    return catalog;
  }

  function paintCell(cell, i) {
    const body = cell.querySelector(".plot-body");
    const item = items[i];
    const site = siteWithRef(pieProtein, pieGeneAa);
    const covered = item.aaByPosition.has(pieProtein + "|" + pieGeneAa);
    const oneStep = covered ? oneStepAt(item) : [];
    if (!covered || !oneStep.length) {
      body.innerHTML = "";
      const p = document.createElement("p");
      p.className = "caveat pie-empty";
      p.textContent = !covered
        ? `${item.row.name} does not cover ${site}.`
        : `No one-step amino acid change from ${item.row.name} at ${site}.`;
      body.appendChild(p);
      return;
    }
    const match = oneStep.find((r) => sameChange(r, picked));
    pieChart(body, oneStep.map((r) => ({
      label: sliceKey(r),
      value: r.mu_aa,
      colour: residueColour(r.alt_aa),
    })), {
      legend: false,
      drawTitle: false,
      width: PIE_W,
      height: PIE_H,
      radius: PIE_R,
      selected: match ? sliceKey(match) : null,
    });
    // The chart caps itself at the viewBox width. In this grid the cell is the limit,
    // so a pie can grow to half or a quarter of the card.
    const svg = body.querySelector("svg");
    if (svg) svg.style.maxWidth = "none";
  }

  function outcomeFor(item) {
    const site = siteWithRef(pieProtein, pieGeneAa);
    const covered = item.aaByPosition.has(pieProtein + "|" + pieGeneAa);
    if (!covered) return { sample: item.row.name, covered: false, site };
    const rows = rowsAt(item);
    const refAa = rows.length ? rows[0].ref_aa : "";
    const refCodon = rows.length ? rows[0].codon : "";
    const row = rows.find((r) => sameChange(r, picked));
    if (!row) return { sample: item.row.name, covered: true, missing: true, refAa, refCodon };
    const change = row.hxb2_change ||
      engine.hxb2Mutation(pieProtein, pieGeneAa, row.ref_aa, row.alt_aa);
    if (row.min_nt_changes === 1) {
      const one = oneStepAt(item);
      const total = one.reduce((t, r) => t + r.mu_aa, 0);
      return {
        sample: item.row.name, covered: true, colour: residueColour(row.alt_aa),
        refAa: row.ref_aa, refCodon: row.codon, aa: sliceKey(row), change, steps: 1,
        pct: total > 0 ? (100 * row.mu_aa / total) : null,
        mu: row.mu_aa,
      };
    }
    const path = bestFarPath(item, row.alt_aa);
    return {
      sample: item.row.name, covered: true, colour: residueColour(row.alt_aa),
      refAa: row.ref_aa, refCodon: row.codon, aa: sliceKey(row), change,
      steps: path ? path.n_nt_changes : row.min_nt_changes,
      pct: null,
      mu: path ? path.mu_path : null,
    };
  }

  function renderSelection() {
    selSlot.innerHTML = "";
    if (pieProtein === null || !picked) return;
    const heading = document.createElement("h3");
    heading.className = "pie-selection-title";
    heading.textContent = "Selected amino acid change";
    const table = document.createElement("table");
    table.className = "data pie-legend pie-by-sequence";
    const thead = document.createElement("thead");
    const headRow = document.createElement("tr");
    const headers = [
      { html: "Sequence", label: "Sequence" },
      { html: "Reference amino acid", label: "Reference amino acid" },
      { html: "Reference codon", label: "Reference codon" },
      { html: "", label: "Color" },
      { html: "Amino acid", label: "Amino acid" },
      { html: "Percent of <em>de novo</em> diversification",
        label: "Percent of de novo diversification" },
      { html: "Mutability (amino acid <em>de novo</em> generation rate)",
        label: "Mutability (amino acid de novo generation rate)" },
    ];
    for (const h of headers) {
      const th = document.createElement("th");
      th.scope = "col";
      if (h.html) th.innerHTML = h.html;
      else th.setAttribute("aria-label", h.label);
      headRow.appendChild(th);
    }
    thead.appendChild(headRow);
    const tbody = document.createElement("tbody");
    for (const item of items) {
      const o = outcomeFor(item);
      const tr = document.createElement("tr");
      const seq = document.createElement("td");
      seq.textContent = o.sample;
      const refCell = document.createElement("td");
      const codonCell = document.createElement("td");
      const colourCell = document.createElement("td");
      const aaCell = document.createElement("td");
      const pctCell = document.createElement("td");
      const rateCell = document.createElement("td");
      rateCell.className = "num";
      const writeRef = () => {
        refCell.textContent = o.refAa || "";
        const code = document.createElement("code");
        code.textContent = o.refCodon || "";
        codonCell.appendChild(code);
      };
      if (!o.covered) {
        aaCell.textContent = `not covered at ${o.site}`;
        aaCell.colSpan = 6;
        tr.append(seq, aaCell);
      } else if (o.missing || !o.aa) {
        writeRef();
        aaCell.textContent = "not reached";
        aaCell.colSpan = 4;
        tr.append(seq, refCell, codonCell, aaCell);
      } else {
        writeRef();
        const swatch = document.createElement("span");
        swatch.className = "pie-swatch";
        swatch.style.background = o.colour;
        colourCell.appendChild(swatch);
        aaCell.textContent = o.steps > 1 ? `${o.aa} (${o.steps}-step)` : o.aa;
        aaCell.title = o.steps > 1
          ? `${o.change}, highest ${o.steps}-step path`
          : o.change;
        pctCell.textContent = o.pct === null ? "\u2014" : `${o.pct.toFixed(1)}%`;
        rateCell.textContent = formatRate(o.mu);
        if (o.steps > 1) rateCell.title = `Highest ${o.steps}-step path rate`;
        tr.append(seq, refCell, codonCell, colourCell, aaCell, pctCell, rateCell);
      }
      tbody.appendChild(tr);
    }
    table.append(thead, tbody);
    selSlot.append(heading, table);
  }

  function fillAaSelect() {
    const catalog = ensurePick();
    aaSelect.innerHTML = "";
    for (const c of catalog) {
      const o = document.createElement("option");
      o.value = c.key;
      o.textContent = c.label;
      aaSelect.appendChild(o);
    }
    aaSelect.value = pickKey(picked);
    aaSelect.disabled = catalog.length === 0;
    aaBar.hidden = catalog.length === 0;
  }

  function drawPies(protein, geneAa) {
    pieProtein = protein;
    pieGeneAa = geneAa;
    picked = null;
    const site = siteWithRef(protein, geneAa);
    const label = engine.hxb2Label(protein, geneAa);
    pieTitle.textContent =
      `Distribution of one-step amino acid changes at ${label}`;
    fillAaSelect();
    syncPieGrid();
    setDetails(pieBoard, pieDetailHtml(site), {
      title: `Distribution of one-step amino acid changes at ${label}`,
    });
    ensurePieDownload();
    renderSelection();
  }

  // Ten pies are five across and two high. The grid stays that tall, and the checked
  // pies divide it: one row when five or fewer fit across, two rows after that, with
  // the second row stretched so it still reaches both edges.
  function tenPieBlockHeight(width) {
    const gap = 8, name = 18, pad = 12;
    const col = (width - 4 * gap) / 5;
    const row = col * (PIE_H / PIE_W) + name + pad;
    return Math.round(2 * row + gap);
  }
  function pieRowCounts(n) {
    if (n <= 5) return [n];
    const cols = Math.ceil(n / 2);
    return [cols, n - cols];
  }
  function sizePieGrid() {
    const width = pieGrid.clientWidth;
    if (!width || !shown.size) return;
    const next = tenPieBlockHeight(width) + "px";
    if (pieGrid.style.height !== next) pieGrid.style.height = next;
  }
  function syncPieGrid() {
    if (pieProtein === null) return;
    const n = shown.size;
    pieGrid.innerHTML = "";
    if (!n) {
      pieGrid.style.height = "";
      const p = document.createElement("p");
      p.className = "caveat pie-empty";
      p.textContent = "No sequences are checked.";
      pieGrid.appendChild(p);
      return;
    }
    const seqs = [];
    items.forEach((item, i) => { if (shown.has(i)) seqs.push(i); });
    let cursor = 0;
    for (const count of pieRowCounts(n)) {
      const row = document.createElement("div");
      row.className = "pie-grid-row";
      for (let k = 0; k < count; k++) {
        const i = seqs[cursor++];
        const item = items[i];
        const cell = document.createElement("figure");
        cell.className = "chart plot-card pie-cell";
        cell.dataset.seq = String(i);
        cell.dataset.exportName = "pie_" + exportSlug(item.row.name, i);
        const name = document.createElement("div");
        name.className = "pie-cell-name";
        name.textContent = item.row.name;
        name.title = item.row.name;
        const body = document.createElement("div");
        body.className = "plot-body";
        cell.append(name, body);
        row.appendChild(cell);
        paintCell(cell, i);
      }
      pieGrid.appendChild(row);
    }
    sizePieGrid();
  }

  function clearPies() {
    pieProtein = null;
    pieGeneAa = null;
    picked = null;
    pieTitle.textContent = "No modelled codon at this nucleotide";
    pieGrid.style.height = "";
    pieGrid.innerHTML = "";
    aaSelect.innerHTML = "";
    aaBar.hidden = true;
    selSlot.innerHTML = "";
  }

  function ensurePieDownload() {
    const row = pieBoard.querySelector(":scope > .plot-tools");
    if (!row || row.querySelector(".pie-grid-download")) return;
    const oneBtn = document.createElement("button");
    oneBtn.type = "button";
    oneBtn.className = "secondary pie-one-download";
    oneBtn.textContent = "Download a specific pie";
    oneBtn.addEventListener("click", openPieChooser);
    const allBtn = document.createElement("button");
    allBtn.type = "button";
    allBtn.className = "secondary pie-grid-download";
    allBtn.textContent = "Download all pies (PNG)";
    allBtn.addEventListener("click", () => downloadAllPies(allBtn));
    row.append(oneBtn, allBtn, pieChooserDialog());
  }

  function pieChooserDialog() {
    const dialog = document.createElement("dialog");
    dialog.className = "details-dialog pie-chooser";
    dialog.setAttribute("aria-label", "Download a specific pie");
    const heading = document.createElement("h3");
    heading.textContent = "Download a specific pie";
    const text = document.createElement("p");
    text.textContent = "Choose the sequence. The image uses the Browse pie: the title, the circle, and the legend table.";
    const lab = document.createElement("label");
    lab.textContent = "Sequence ";
    const sel = document.createElement("select");
    sel.className = "pie-chooser-select";
    lab.appendChild(sel);
    const actions = document.createElement("div");
    actions.className = "pie-chooser-actions";
    const save = document.createElement("button");
    save.type = "button";
    save.textContent = "Download";
    const cancel = document.createElement("button");
    cancel.type = "button";
    cancel.className = "secondary";
    cancel.textContent = "Cancel";
    cancel.addEventListener("click", () => dialog.close());
    save.addEventListener("click", async () => {
      const item = items[Number(sel.value)];
      if (!item) return;
      save.disabled = true;
      save.textContent = "Saving\u2026";
      try {
        await downloadOnePie(item);
        dialog.close();
      } catch (err) {
        showError(head, "Could not download the pie. " + err.message);
      } finally {
        save.disabled = false;
        save.textContent = "Download";
      }
    });
    actions.append(save, cancel);
    dialog.append(heading, text, lab, actions);
    return dialog;
  }

  function openPieChooser() {
    const dialog = pieBoard.querySelector("dialog.pie-chooser");
    const sel = dialog && dialog.querySelector("select");
    if (!dialog || !sel) return;
    const choices = items.filter((item) => oneStepAt(item).length);
    sel.innerHTML = "";
    for (const item of choices) {
      const o = document.createElement("option");
      o.value = String(items.indexOf(item));
      o.textContent = item.row.name;
      sel.appendChild(o);
    }
    if (!choices.length) {
      showError(head, "None of the sequences has a one-step pie at this position.");
      return;
    }
    if (typeof dialog.showModal === "function") dialog.showModal();
    else dialog.setAttribute("open", "");
  }

  /** A Browse-format pie for one sequence: title, circle, codon line, and legend table. */
  async function downloadOnePie(item) {
    const one = oneStepAt(item);
    if (!one.length) throw new Error(`${item.row.name} has no one-step pie at this position`);
    const site = siteForSequence(item);
    const who = item.row.name;
    const fig = document.createElement("figure");
    fig.className = "chart plot-card pie-export";
    fig.dataset.exportName = "pie_" + exportSlug(who, items.indexOf(item));
    const headEl = document.createElement("div");
    headEl.className = "plot-head";
    const title = document.createElement("h3");
    title.className = "plot-title";
    title.textContent = `Distribution of one-step amino acid changes from ${who} at ${site}`;
    headEl.appendChild(title);
    const body = document.createElement("div");
    body.className = "plot-body";
    const note = document.createElement("p");
    note.className = "plot-note";
    note.innerHTML = sequenceNote(item);
    const table = legendTable(item, one);
    fig.append(headEl, body, note, table);
    const stage = document.createElement("div");
    stage.className = "pie-export-stage";
    stage.appendChild(fig);
    document.body.appendChild(stage);
    const match = one.find((r) => sameChange(r, picked));
    const pieW = Math.floor((LAYOUT.width - LAYOUT.gap) / 2);
    try {
      pieChart(body, one.map((r) => ({
        label: sliceKey(r),
        value: r.mu_aa,
        colour: residueColour(r.alt_aa),
      })), {
        legend: false,
        drawTitle: false,
        width: pieW,
        height: 500,
        title: title.textContent,
        selected: match ? sliceKey(match) : null,
      });
      await downloadFigure(fig, { release: page.release });
    } finally {
      stage.remove();
    }
  }

  function siteForSequence(item) {
    const base = engine.hxb2Label(pieProtein, pieGeneAa);
    const rows = rowsAt(item);
    const aa = rows.length ? rows[0].ref_aa : "";
    const m = /^(.*\s)(\d+)$/.exec(String(base).trim());
    return aa && m ? m[1] + aa + m[2] : base;
  }

  function sequenceNote(item) {
    const rows = rowsAt(item);
    const rec = assets.reference.proteins[pieProtein];
    const refCodon = rec.codon_seq[pieGeneAa - rec.first_aa];
    const codonNt = rec.codon_nt[pieGeneAa - rec.first_aa];
    const yourCodon = rows.length ? rows[0].codon : refCodon;
    const aa = rows.length ? rows[0].ref_aa : (GENETIC_CODE[refCodon] || "");
    const hxb2 = engine.hxb2Label(pieProtein, pieGeneAa);
    const codonBit = yourCodon === refCodon
      ? `codon <code>${escapeHtml(yourCodon)}</code>`
      : `codon <code>${escapeHtml(yourCodon)}</code> (NL4-3 has <code>${escapeHtml(refCodon)}</code>)`;
    return `<strong>${escapeHtml(hxb2)}</strong> &mdash; <strong>${escapeHtml(item.row.name)}</strong>, ` +
      `amino acid <strong>${escapeHtml(aa)}</strong>, ${codonBit}. ` +
      `NL4-3 ${escapeHtml(pieProtein)} position ${pieGeneAa}; ` +
      `nucleotides ${codonNt[0].toLocaleString()}\u2013${codonNt[2].toLocaleString()}.`;
  }

  function legendTable(item, one) {
    const total = one.reduce((t, r) => t + r.mu_aa, 0);
    const table = document.createElement("table");
    table.className = "data pie-legend is-selectable";
    const thead = document.createElement("thead");
    const headRow = document.createElement("tr");
    const headers = [
      { html: "", label: "Color" },
      { html: "Amino acid", label: "Amino acid" },
      { html: "Percent of <em>de novo</em> diversification",
        label: "Percent of de novo diversification" },
      { html: "Mutability (amino acid <em>de novo</em> generation rate)",
        label: "Mutability (amino acid de novo generation rate)" },
    ];
    for (const h of headers) {
      const th = document.createElement("th");
      th.scope = "col";
      if (h.html) th.innerHTML = h.html;
      else th.setAttribute("aria-label", h.label);
      headRow.appendChild(th);
    }
    thead.appendChild(headRow);
    const tbody = document.createElement("tbody");
    for (const r of one) {
      const tr = document.createElement("tr");
      if (sameChange(r, picked)) tr.classList.add("is-selected");
      const colourCell = document.createElement("td");
      const swatch = document.createElement("span");
      swatch.className = "pie-swatch";
      swatch.style.background = residueColour(r.alt_aa);
      colourCell.appendChild(swatch);
      const aaCell = document.createElement("td");
      aaCell.textContent = sliceKey(r);
      if (r.synonymous) aaCell.title = "Synonymous change";
      const pctCell = document.createElement("td");
      pctCell.textContent = total > 0 ? `${((r.mu_aa / total) * 100).toFixed(1)}%` : "\u2014";
      const rateCell = document.createElement("td");
      rateCell.className = "num";
      rateCell.textContent = formatRate(r.mu_aa);
      tr.append(colourCell, aaCell, pctCell, rateCell);
      tbody.appendChild(tr);
    }
    table.append(thead, tbody);
    return table;
  }

  async function downloadAllPies(button) {
    if (typeof html2canvas === "undefined") {
      showError(head, "The PNG library did not load. Reload the page and try again.");
      return;
    }
    const prev = button.textContent;
    button.disabled = true;
    button.textContent = "Generating\u2026";
    try {
      const canvas = await html2canvas(pieBoard, {
        backgroundColor: "#ffffff",
        scale: 2,
        logging: false,
        useCORS: true,
        ignoreElements: (node) => !!(node.hasAttribute && node.hasAttribute("data-export-skip")),
      });
      const blob = await new Promise((res) => canvas.toBlob(res, "image/png"));
      if (!blob) throw new Error("the canvas could not be converted to an image");
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = `mra_pies_${new Date().toISOString().split("T")[0]}.png`;
      a.click();
      URL.revokeObjectURL(url);
    } catch (err) {
      showError(head, "Could not download the pies. " + err.message);
    } finally {
      button.disabled = false;
      button.textContent = prev;
    }
  }

  function apply(protein, geneAa, preferNt) {
    const item = items[active];
    const rec = assets.reference.proteins[protein];
    const codonNt = rec.codon_nt[geneAa - rec.first_aa];
    selected = { protein, geneAa };
    markerNt = (preferNt && codonNt.includes(preferNt)) ? preferNt : codonNt[1];
    elSlider.min = String(zoom0);
    elSlider.max = String(zoom1);
    elSlider.value = String(markerNt);
    head.innerHTML = "";

    elFrame.innerHTML = "";
    for (const f of (item.frameAt.get(markerNt) || [])) {
      const o = document.createElement("option");
      o.value = f.protein + "|" + f.geneAa;
      o.textContent = `${engine.hxb2Label(f.protein, f.geneAa)}  (${f.protein} ${f.geneAa})`;
      elFrame.appendChild(o);
    }
    elFrame.value = protein + "|" + geneAa;
    elFrame.disabled = (item.frameAt.get(markerNt) || []).length < 2;

    const t = targetFor(item.targets, protein, geneAa);
    if (t && elSys.value === "protein") {
      elTarget.value = t.key;
      elPos.value = String(engine.ownPositionFor(t, geneAa));
    }

    drawTrack();
    drawPies(protein, geneAa);
  }

  function selectNt(nt) {
    const item = items[active];
    markerNt = Math.min(Math.max(Math.round(nt), region.first_position), region.last_position);
    elSlider.value = String(markerNt);
    if (elSys.value === "nt") elPos.value = String(markerNt);
    if (!item.perSite.has(markerNt) || !(item.frameAt.get(markerNt) || []).length) {
      elFrame.innerHTML = "";
      elFrame.disabled = true;
      showError(head,
        item.perSite.has(markerNt)
          ? `NL4-3 nucleotide ${markerNt.toLocaleString()} is covered by ${item.row.name} ` +
            `but is not inside a codon this atlas models.`
          : `NL4-3 nucleotide ${markerNt.toLocaleString()} is not covered by ${item.row.name}. ` +
            `It covers ${item.sites[0].toLocaleString()}\u2013${item.sites[item.sites.length - 1].toLocaleString()}.`);
      const donor = items.find((it) => (it.frameAt.get(markerNt) || []).length);
      const frame = donor ? donor.frameAt.get(markerNt)[0] : null;
      selected = frame ? { protein: frame.protein, geneAa: frame.geneAa } : null;
      drawTrack();
      if (frame) drawPies(frame.protein, frame.geneAa);
      else clearPies();
      return;
    }
    const frames = item.frameAt.get(markerNt) || [];
    const keep = selected && frames.find((f) => f.protein === selected.protein);
    apply((keep || frames[0]).protein, (keep || frames[0]).geneAa, markerNt);
  }

  function jump(zoomGene) {
    zoomGene = zoomGene === true;
    const item = items[active];
    head.innerHTML = "";
    if (elSys.value === "variant") {
      if (!item.variants.length) {
        showError(head, `${item.row.name} matches NL4-3 at every modelled position.`);
        return;
      }
      const want = Number(elPos.value);
      item.variantIdx = Number.isInteger(want) && want >= 1 && want <= item.variants.length
        ? want - 1 : (item.variantIdx + 1) % item.variants.length;
      elPos.value = String(item.variantIdx + 1);
      const nt = item.variants[item.variantIdx];
      if (zoomGene) zoomToNt(item, nt);
      selectNt(nt);
      return;
    }
    const raw = Number(elPos.value);
    if (!Number.isInteger(raw) || raw < 1) { showError(head, "Enter a positive whole number."); return; }
    if (elSys.value === "nt") {
      if (!item.perSite.has(raw)) {
        showError(head,
          `NL4-3 nucleotide ${raw.toLocaleString()} is not covered by ${item.row.name}. ` +
          `It covers ${item.sites[0].toLocaleString()}\u2013${item.sites[item.sites.length - 1].toLocaleString()}.`);
        return;
      }
      if (zoomGene) zoomToNt(item, raw);
      selectNt(raw);
      return;
    }
    const t = new Map(item.targets.map((x) => [x.key, x])).get(elTarget.value);
    if (!t) { showError(head, "Pick a protein or subunit."); return; }
    const geneAa = engine.geneAaFromOwnPosition(t, raw);
    if (geneAa === null || !item.aaByPosition.has(t.protein + "|" + geneAa)) {
      const span = item.proteinSpans.get(t.protein);
      showError(head, `${t.key} ${raw} is not covered by ${item.row.name}` +
        (span ? ` \u2014 ${t.protein} ${span.min}\u2013${span.max} was analysed.` : "."));
      return;
    }
    if (zoomGene) zoomToGene(t.protein, geneAa, t);
    apply(t.protein, geneAa, null);
  }

  function renderIndelTable(slot, kind) {
    const noun = kind === "insertion" ? "insertion" : "deletion";
    const heading = document.createElement("h2");
    heading.textContent = kind === "insertion" ? "Insertions" : "Deletions";
    slot.appendChild(heading);
    const note = document.createElement("p");
    note.className = "caveat";
    note.textContent =
      "One row per nucleotide a sequence covers. The rate is the NL4-3 " + noun +
      " rate at that position. It does not depend on the codon or the surrounding bases, " +
      "so every sequence that covers a position has the same rate there. These are " +
      "single-nucleotide events, not the amino-acid changes charted for substitutions.";
    slot.appendChild(note);
    const rows = [];
    for (const item of items) {
      for (const pos of item.sites) {
        const rec = item.indels.get(pos);
        if (!rec || !(rec[kind] > 0)) continue;
        rows.push({
          sample: item.row.name,
          ref_position: pos,
          nucleotide: item.baseAt.get(pos) || "",
          rate: rec[kind],
        });
      }
    }
    const box = document.createElement("div");
    slot.appendChild(box);
    renderTable(box, {
      columns: [
        { key: "sample", label: "Sample" },
        { key: "ref_position", label: "NL4-3 position", type: "num" },
        { key: "nucleotide", label: "Nucleotide",
          title: "The base in the submitted sequence at this aligned position" },
        { key: "rate", label: "Rate", type: "rate",
          title: "NL4-3 " + noun + " rate at this position" },
      ],
      rows,
      sort: { key: "ref_position", dir: "asc" },
      csvHeader: "label",
      filterPlaceholder: "filter by sample, position, or nucleotide",
      downloadName: `mra_${page.version}_${kind}s.csv`,
      provenance: exportProvenance(page.release, [
        `NL4-3 ${noun} rates at positions covered by the submitted sequences.`,
        "An insertion or deletion has no nucleotide context, so the rate is the reference rate at the aligned coordinate.",
      ]),
    });
  }

  let analysisMode = "aa";
  function renderLower() {
    lower.innerHTML = "";
    if (truncated) lower.appendChild(truncationNote());
    for (const r of allRows.filter((x) => x.status !== "analysed")) {
      const div = document.createElement("div");
      div.className = "error-box";
      div.innerHTML = r.status === "not placed"
        ? `<strong>${escapeHtml(r.name)} was not placed on the reference.</strong> ` +
          escapeHtml(r.summary.reason)
        : `<strong>${escapeHtml(r.name)} lies entirely outside the modelled region.</strong> ` +
          `The atlas covers NL4-3 reference positions ${region.first_position.toLocaleString()}` +
          `&ndash;${region.last_position.toLocaleString()}; this sequence aligned beyond it, so ` +
          `no position-resolved rates can be given. ${NEF_NOTE}`;
      lower.appendChild(div);
    }
    if (analysisMode === "insertion" || analysisMode === "deletion") {
      renderIndelTable(lower, analysisMode);
      return;
    }
    renderCompiledTables(lower, items);
  }

  aaSelect.addEventListener("change", () => {
    picked = pickFromKey(aaSelect.value);
    pieGrid.querySelectorAll(".pie-cell").forEach((cell) => paintCell(cell, Number(cell.dataset.seq)));
    renderSelection();
  });
  elWin.addEventListener("input", () => { elWinOut.textContent = elWin.value; });
  elWin.addEventListener("change", drawTrack);
  elClass.addEventListener("change", () => {
    drawTrack();
    const indel = elClass.value === "insertion" || elClass.value === "deletion";
    pieBoard.hidden = indel;
    const mode = indel ? elClass.value : "aa";
    if (mode === analysisMode) return;
    analysisMode = mode;
    renderLower();
  });
  ctrl.querySelector(".sq-zoom-reset").addEventListener("click", resetZoom);
  elSlider.addEventListener("input", () => selectNt(Number(elSlider.value)));
  ctrl.querySelector(".sq-go").addEventListener("click", () => jump(true));
  elPos.addEventListener("keydown", (e) => { if (e.key === "Enter") jump(true); });
  elTarget.addEventListener("change", jump);
  elFrame.addEventListener("change", (e) => {
    const [p, g] = e.target.value.split("|");
    apply(p, Number(g), markerNt);
  });
  elSys.addEventListener("change", () => {
    const item = items[active];
    const mode = elSys.value;
    elTarget.disabled = mode !== "protein";
    const t = selected ? targetFor(item.targets, selected.protein, selected.geneAa) : null;
    elPos.value = mode === "nt" ? String(markerNt)
      : mode === "variant" ? String(item.variantIdx + 1)
      : (t ? String(engine.ownPositionFor(t, selected.geneAa)) : "1");
  });

  fillTargets();
  markSelected();
  elTarget.disabled = false;
  if (markerNt !== null) selectNt(markerNt);
  renderLower();
}


/** HXB2 site with NL4-3's residue, e.g. "RT K65". */
function siteWithRef(protein, geneAa) {
  const base = engine.hxb2Label(protein, geneAa);
  const rec = assets.reference.proteins[protein];
  const codon = rec && rec.codon_seq ? rec.codon_seq[geneAa - rec.first_aa] : "";
  const refAa = GENETIC_CODE[codon] || "";
  const m = /^(.*\s)(\d+)$/.exec(String(base).trim());
  return refAa && m ? m[1] + refAa + m[2] : base;
}

function exportSlug(name, i) {
  const s = String(name || "sequence").toLowerCase()
    .replace(/[^a-z0-9]+/g, "_").replace(/^_|_$/g, "");
  return (s.slice(0, 40) || "sequence") + "_" + (i + 1);
}

function pieDetailHtml(site) {
  return `<p>Each circle is one submitted sequence at <strong>${escapeHtml(site)}</strong>. ` +
    `A slice is an amino acid reachable by one nucleotide substitution from that sequence's codon. ` +
    `Its size is that amino acid's share of the one-step <em>de novo</em> diversification: ` +
    `the amino-acid <em>de novo</em> generation rate divided by the sum of those rates for every ` +
    `one-step amino acid. Colours name the destination residue.</p>` +
    `<p>The amino acid list chooses one change for every pie. A pie outlines that slice when the ` +
    `sequence can reach it in one nucleotide substitution. The table under the grid is that change ` +
    `for every sequence, in the same columns as the Browse legend.</p>` +
    `<p><strong>Percent of <em>de novo</em> diversification</strong> is the change's share of the ` +
    `one-step amino-acid generation rates. <strong>Mutability</strong> for a one-step change is the ` +
    `amino-acid <em>de novo</em> generation rate: the sum, per replication, of every nucleotide path ` +
    `to that amino acid. When a sequence cannot reach the change in one step, the row is the ` +
    `highest-rate path of two or three nucleotide substitutions, the percent column is blank, and ` +
    `mutability is that path's rate.</p>` +
    `<p>These rates describe generation in one replication, not establishment in a population. ` +
    `A sequence that does not cover this codon is named in its cell. Unchecking a sequence ` +
    `hides its line and its pie; the table still lists every sequence.</p>`;
}

// The track overlays NL4-3 from the reference's own custom atlas, computed once from assets
// already loaded. Cached so a batch of sequences pays the cost once.
/** A legend-length name: FASTA headers are routinely 40 characters and sit at the line end. */
function shortName(name) {
  const n = String(name || "your sequence").trim();
  return n.length <= 18 ? n : n.slice(0, 17) + "\u2026";
}

let _refCache = null;
function indelRatesFor(positions) {
  const out = new Map();
  for (const pos of positions) {
    const ins = engine.indelRate(pos, "insertion");
    const del = engine.indelRate(pos, "deletion");
    out.set(pos, {
      insertion: ins ? ins.rate : 0,
      deletion: del ? del.rate : 0,
    });
  }
  return out;
}
function referenceAtlas() {
  if (_refCache) return _refCache;
  const seq = assets.reference.sequence;
  const refAtlas = engine.customAtlas({
    name: "NL4-3", sequence: seq, alignment: alignToReference(seq, seq) });
  const perSite = perSiteSupply(groupByPosition(refAtlas.nt), "rate");
  _refCache = { perSite, indels: indelRatesFor([...perSite.keys()]) };
  return _refCache;
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

/** One amino-acid table and one codon-path table, rows from every submitted sequence.
 *  The columns match Browse; Sample is the extra first column. */
function renderCompiledTables(slot, items) {
  const aaRows = [];
  const pathRows = [];
  for (const item of items) {
    const sample = item.row.name;
    for (const r of item.detail.aa) {
      aaRows.push({
        sample,
        hxb2_position: r.hxb2_label,
        ref_aa: r.ref_aa,
        alt_aa: r.alt_aa,
        aa_change_hxb2: r.hxb2_change,
        mu_aa: r.mu_aa,
        mu_aa_lo95: r.mu_aa_lo95,
        mu_aa_hi95: r.mu_aa_hi95,
        min_nt_changes: r.min_nt_changes,
        n_paths: r.n_paths,
        synonymous: r.synonymous ? "yes" : "no",
      });
    }
    for (const r of item.detail.paths) {
      pathRows.push({
        sample,
        hxb2_position: engine.hxb2Label(r.protein, r.position),
        codon: r.codon,
        mutant_codon: r.mutant_codon,
        n_nt_changes: r.n_nt_changes,
        mutation_type_path: r.mutation_type_path,
        mu_path: r.mu_path,
        mu_path_lo95: r.mu_path_lo95,
        mu_path_hi95: r.mu_path_hi95,
      });
    }
  }
  const provenance = compiledProvenance(items);

  const allData = document.createElement("h2");
  allData.textContent = "All Data";
  slot.appendChild(allData);

  const aaHeading = document.createElement("h3");
  aaHeading.textContent = "Amino-acid mutability";
  slot.appendChild(aaHeading);
  const aaSlot = document.createElement("div");
  slot.appendChild(aaSlot);
  renderTable(aaSlot, {
    columns: [
      { key: "sample", label: "Sample", title: "The submitted sequence this row was computed from" },
      { key: "hxb2_position", label: "HXB2 Position", title: "HXB2 amino-acid numbering" },
      { key: "ref_aa", label: "Reference (From)" },
      { key: "alt_aa", label: "Mutation (To)" },
      { key: "aa_change_hxb2", label: "Amino Acid Change",
        title: "Protein, reference residue, HXB2 position, and destination residue" },
      { key: "mu_aa", label: "Amino Acid Mutability (mu_aa)", type: "rate",
        title: "Summed rate over all paths to this amino acid" },
      { key: "mu_aa_lo95", label: "95% Lower", type: "rate" },
      { key: "mu_aa_hi95", label: "95% Upper", type: "rate" },
      { key: "min_nt_changes", label: "Minimum Nucleotide", type: "num" },
      { key: "n_paths", label: "Paths", type: "num" },
      { key: "synonymous", label: "Synonymous" },
    ],
    rows: aaRows,
    sort: { key: "mu_aa", dir: "desc" },
    filterPlaceholder: "filter by sample, position, change\u2026",
    csvHeader: "label",
    downloadName: `mra_${page.version}_aa.csv`,
    provenance,
  });

  const pathHeading = document.createElement("h3");
  pathHeading.textContent = "Codon paths";
  slot.appendChild(pathHeading);
  const pathSlot = document.createElement("div");
  slot.appendChild(pathSlot);
  renderTable(pathSlot, {
    columns: [
      { key: "sample", label: "Sample", title: "The submitted sequence this row was computed from" },
      { key: "hxb2_position", label: "HXB2 Position", title: "HXB2 amino-acid numbering" },
      { key: "codon", label: "Reference Codon" },
      { key: "mutant_codon", label: "Mutant Codon" },
      { key: "n_nt_changes", label: "Nucleotide Changes", type: "num" },
      { key: "mutation_type_path", label: "Substitutions" },
      { key: "mu_path", label: "Mutation Rate", type: "rate" },
      { key: "mu_path_lo95", label: "95% Lower", type: "rate" },
      { key: "mu_path_hi95", label: "95% Upper", type: "rate" },
    ],
    rows: pathRows,
    sort: { key: "mu_path", dir: "desc" },
    filterPlaceholder: "filter by sample, codon, substitution\u2026",
    csvHeader: "label",
    downloadName: `mra_${page.version}_paths.csv`,
    provenance,
  });
  pathSlot.appendChild(caveat(
    "A path rate is the product of the per-site rates at the changed positions; its interval " +
    "is the same product over the marginal bounds. Amino-acid mutability above is the sum over " +
    "every path to that amino acid, including 2- and 3-substitution routes."));
}

/** Provenance for a merged export: which model, which version, and what each sample rests on. */
function compiledProvenance(items) {
  const lines = [
    `Custom atlas for ${items.length} submitted sequence(s). Rows from every sequence are merged; ` +
      `the sample column names which one a row came from.`,
    `Retrieved from the MRA portal; computed entirely in the browser.`,
  ];
  for (const item of items) {
    const s = item.row.summary;
    const geom = item.geom;
    const identity = s.identity == null ? "n/a" : `${(s.identity * 100).toFixed(3)}%`;
    lines.push(
      `Sample ${item.row.name}: NL4-3 ${geom.refStart}-${geom.refEnd}, identity ${identity}, ` +
      `${geom.orientation}. Dropped outside ${region.first_position}-${region.last_position}: ` +
      `${s.droppedOutsideRegion}. Spline fallback at ${s.fallbackPositions} position(s).`);
  }
  lines.push(
    "THIS IS NOT A MEASUREMENT OF THESE VIRUSES. It is the MRA model, fitted on NL4-3,",
    "applied to each submitted sequence. Codon identity and local 3' context come from that",
    "sequence; position effects were measured on NL4-3 and transfer only where it aligns.",
    "The amino-acid and codon-path tables exclude insertions and deletions.");
  return exportProvenance(page.release, lines);
}

// Restore anything carried over from a previous run in this tab.
const carried = sessionStorage.getItem("mra.sequences");
if (carried && !fasta.value) fasta.value = carried;
