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
// Up to MAX_DETAIL_SEQUENCES submissions share one genome track. Checkboxes choose which
// lines are drawn; one selector chooses which sequence the pie and the flow chart use.
// Under those charts the tables match Browse. Every submitted sequence is in both,
// distinguished by a sample column, and each table downloads in full.

import { bootstrapPage, privacyPanel, note, caveat, showError, exportProvenance,
         CUSTOM_CAVEAT, NEF_NOTE } from "./page.js";
import { loadSmallAssets } from "./data.js";
import { Analysis, progressRow, MAX_DETAIL_SEQUENCES } from "./analysis.js";
import { parseFasta, alignToReference } from "./align.js";
import { renderTable } from "./table.js";
import { MraEngine, GENETIC_CODE } from "./mra_engine.js";
import { perSiteSupply, groupByPosition, buildFrameIndex, geneSegments, targetFor,
         smoothSeries } from "./explore_core.js";
import { drawSupplyTrack, drawPositionPie, CLASS_LABEL } from "./track_view.js";
import { loadLanlFlow, lanlFlowAt } from "./lanl_flow.js";
import { attachFigureDownloads } from "./figure_export.js";
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

function syncFigures() {
  attachFigureDownloads({
    release: page.release,
    onError: (err) => console.error("figure export failed:", err),
  });
}

let lanlPromise = null;
function lanlIndex() {
  if (!lanlPromise) lanlPromise = loadLanlFlow(page.version);
  return lanlPromise;
}

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
 * One genome track for the whole submission, then the browse cards. Checkboxes choose
 * which sequences are drawn. The sample selector sits inside the pie chart and chooses
 * which sequence the pie and the flow chart analyse. The smoothing and coordinate
 * controls sit between the Results heading and the plot, and the position slider sits
 * in the track between the gene
 * boxes and the rate curve, matching Browse.
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
  let flowToken = 0;

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
    });
    const swatch = document.createElement("i");
    swatch.className = "swatch";
    swatch.style.background = item.colour;
    const name = document.createElement("span");
    name.textContent = item.row.name;
    lab.append(cb, swatch, name);
    toggles.appendChild(lab);
  });
  fig.append(stack, figCap, toggles);

  const ctrl = document.createElement("div");
  ctrl.className = "control-bar";
  ctrl.setAttribute("data-export-skip", "");
  ctrl.innerHTML =
    `<div><label>Smoothing window (nt)</label>` +
    `<input type="range" class="sq-win" min="25" max="501" step="2" value="151">` +
    `<output class="sq-win-out">151</output></div>` +
    `<div><label>Supply shown</label><select class="sq-class">` +
    `<option value="total">Total substitution supply</option>` +
    `<option value="nonga">non-G&rarr;A only</option>` +
    `<option value="ga">G&rarr;A only</option></select></div>` +
    `<div><label>Jump to</label><select class="sq-sys">` +
    `<option value="protein">Protein position</option>` +
    `<option value="nt">NL4-3 nucleotide</option>` +
    `<option value="variant">Difference from NL4-3</option></select></div>` +
    `<div><label>Protein or subunit</label><select class="sq-target"></select></div>` +
    `<div><label>Position</label><input type="number" class="sq-pos" min="1" size="8"></div>` +
    `<div><label>Reading frame here</label><select class="sq-frame"></select></div>` +
    `<div class="bar-actions"><button class="sq-go">Go</button></div>`;

  const head = document.createElement("div");
  const dash = document.createElement("figure");
  dash.className = "dashboard-figures";
  dash.setAttribute("data-export-name", "position_dashboard_supply");
  const chartEl = document.createElement("div");
  const captionEl = document.createElement("figcaption");
  dash.append(chartEl, captionEl);
  const lower = document.createElement("div");
  slot.append(ctrl, fig, head, dash, lower);

  const elWin = ctrl.querySelector(".sq-win"), elWinOut = ctrl.querySelector(".sq-win-out");
  const elClass = ctrl.querySelector(".sq-class"), elSys = ctrl.querySelector(".sq-sys");
  const elTarget = ctrl.querySelector(".sq-target"), elPos = ctrl.querySelector(".sq-pos");
  const elFrame = ctrl.querySelector(".sq-frame");
  // Lives inside the pie, not the control bar. drawPositionPie reattaches this node after
  // each redraw, so the listener below survives a new position or a new sample.
  const sampleWrap = document.createElement("div");
  sampleWrap.className = "pie-sample";
  sampleWrap.setAttribute("data-export-skip", "");
  const elWhich = document.createElement("select");
  elWhich.className = "sq-which";
  elWhich.id = "pie-sample";
  elWhich.setAttribute("aria-label", "Pie Chart Sample");
  const sampleLab = document.createElement("label");
  sampleLab.htmlFor = "pie-sample";
  sampleLab.textContent = "Pie Chart Sample";
  sampleWrap.append(sampleLab, elWhich);
  items.forEach((item, i) => {
    const o = document.createElement("option");
    o.value = String(i);
    o.textContent = item.row.name;
    elWhich.appendChild(o);
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

  function drawTrack() {
    const item = items[active];
    const win = Number(elWin.value);
    elWinOut.textContent = String(win);
    const cls = elClass.value;
    const series = [];
    items.forEach((it, i) => {
      if (!shown.has(i)) return;
      series.push({
        label: shortName(it.row.name),
        colour: it.colour,
        points: smoothSeries(it.sites, it.perSite, cls, win),
        width: i === active ? 1.7 : 1.15,
      });
    });
    const trackText = drawSupplyTrack(trackSlot, {
      perSite: item.perSite, sites: item.sites, series, win, cls, segments: SEGMENTS,
      first: region.first_position, last: region.last_position,
      marker: markerNt,
      markerLabel: selected ? engine.hxb2Label(selected.protein, selected.geneAa) : String(markerNt),
      variants: item.variants,
      title: `Smoothed ${CLASS_LABEL[cls]} along the submitted sequences, with NL4-3 overlaid`,
      source: "computed from the codons and 3\u2032 context in the sequences you submitted, not from NL4-3",
      refPerSite: referenceAtlas().perSite,
      subjectName: shortName(item.row.name),
      refBase: item.baseAt.get(markerNt),
      onSelect: (pos) => selectNt(pos),
      // Room between the gene boxes and the rate curve for the position slider.
      sliderGap: 28,
    });
    placeTrackSlider();
    const where = selected ? ` (${engine.hxb2Label(selected.protein, selected.geneAa)})` : "";
    figCap.innerHTML =
      `Genome-wide ${CLASS_LABEL[cls]}, ${win} nt mean. Marker at ` +
      `${markerNt.toLocaleString()}${where}. Red ticks are where ` +
      `<strong>${escapeHtml(shortName(item.row.name))}</strong> differs from NL4-3. ` +
      `The dashed grey line is NL4-3.`;
    setDetails(fig, trackText, { title: "The genome supply track" });
    syncFigures();
  }

  function showPie(item, rows, label, noteHtml, yourCodon, protein, geneAa, refAa) {
    const token = ++flowToken;
    const opts = {
      custom: true, noteHtml, refCodon: yourCodon, refAa,
      sourceName: item.row.name,
      markerPhrase: "the codon in " + escapeHtml(item.row.name),
      sampleControl: sampleWrap,
    };
    drawPositionPie(chartEl, captionEl, rows, label, opts);
    syncFigures();
    lanlIndex().then((index) => {
      if (token !== flowToken) return;
      const flow = index ? lanlFlowAt(index, protein, geneAa) : null;
      if (!flow) return;
      drawPositionPie(chartEl, captionEl, rows, label, { ...opts, flow });
      syncFigures();
    }).catch((err) => console.error("flow chart failed:", err));
  }

  function apply(protein, geneAa, preferNt) {
    const item = items[active];
    const rec = assets.reference.proteins[protein];
    const codonNt = rec.codon_nt[geneAa - rec.first_aa];
    selected = { protein, geneAa };
    markerNt = (preferNt && codonNt.includes(preferNt)) ? preferNt : codonNt[1];
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

    const aaRows = item.aaByPosition.get(protein + "|" + geneAa) || [];
    const label = engine.hxb2Label(protein, geneAa);
    const yourCodon = aaRows.length ? aaRows[0].codon : "?";
    const refCodon = rec.codon_seq[geneAa - rec.first_aa];
    const aa = aaRows.length ? aaRows[0].ref_aa : "?";
    const refAa = GENETIC_CODE[refCodon] || "";
    const own = t ? engine.ownPositionFor(t, geneAa) : null;
    const noteHtml =
      `<strong>${escapeHtml(label)}</strong> &mdash; <strong>${escapeHtml(item.row.name)}</strong>, ` +
      `codon <code>${escapeHtml(yourCodon)}</code>` +
      (yourCodon !== refCodon ? ` (NL4-3 has <code>${escapeHtml(refCodon)}</code>)` : " (same as NL4-3)") +
      `, amino acid <strong>${escapeHtml(aa)}</strong>. ` +
      `${escapeHtml(protein)} position ${geneAa}` +
      (t && t.offset ? ` (${escapeHtml(t.key)} ${own})` : "") +
      `; NL4-3 reference nucleotides ${codonNt[0].toLocaleString()}\u2013${codonNt[2].toLocaleString()}.`;
    drawTrack();
    showPie(item, aaRows, label, noteHtml, yourCodon, protein, geneAa, refAa);
  }

  function selectNt(nt) {
    const item = items[active];
    markerNt = Math.min(Math.max(Math.round(nt), region.first_position), region.last_position);
    elSlider.value = String(markerNt);
    if (elSys.value === "nt") elPos.value = String(markerNt);
    if (!item.perSite.has(markerNt) || !(item.frameAt.get(markerNt) || []).length) {
      selected = null;
      elFrame.innerHTML = "";
      elFrame.disabled = true;
      showError(head,
        item.perSite.has(markerNt)
          ? `NL4-3 nucleotide ${markerNt.toLocaleString()} is covered by ${item.row.name} ` +
            `but is not inside a codon this atlas models.`
          : `NL4-3 nucleotide ${markerNt.toLocaleString()} is not covered by ${item.row.name}. ` +
            `It covers ${item.sites[0].toLocaleString()}\u2013${item.sites[item.sites.length - 1].toLocaleString()}.`);
      chartEl.innerHTML = "";
      captionEl.textContent = "";
      drawTrack();
      return;
    }
    const frames = item.frameAt.get(markerNt) || [];
    const keep = selected && frames.find((f) => f.protein === selected.protein);
    apply((keep || frames[0]).protein, (keep || frames[0]).geneAa, markerNt);
  }

  function jump() {
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
      selectNt(item.variants[item.variantIdx]);
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
    apply(t.protein, geneAa, null);
  }

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
    renderCompiledTables(lower, items);
  }

  elWin.addEventListener("input", () => { elWinOut.textContent = elWin.value; });
  elWin.addEventListener("change", drawTrack);
  elClass.addEventListener("change", drawTrack);
  elSlider.addEventListener("input", () => selectNt(Number(elSlider.value)));
  ctrl.querySelector(".sq-go").addEventListener("click", jump);
  elPos.addEventListener("keydown", (e) => { if (e.key === "Enter") jump(); });
  elTarget.addEventListener("change", jump);
  elFrame.addEventListener("change", (e) => {
    const [p, g] = e.target.value.split("|");
    apply(p, Number(g), markerNt);
  });
  elWhich.addEventListener("change", () => {
    active = Number(elWhich.value);
    markSelected();
    fillTargets();
    const item = items[active];
    const nt = nearestCodon(item, markerNt);
    if (nt === null) {
      showError(head, `${item.row.name} has no modelled codon to chart.`);
      chartEl.innerHTML = "";
      captionEl.textContent = "";
    } else {
      selectNt(nt);
    }
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


// The track overlays NL4-3 from the reference's own custom atlas, computed once from assets
// already loaded. Cached so a batch of sequences pays the cost once.
/** A legend-length name: FASTA headers are routinely 40 characters and sit at the line end. */
function shortName(name) {
  const n = String(name || "your sequence").trim();
  return n.length <= 18 ? n : n.slice(0, 17) + "\u2026";
}

let _refCache = null;
function referenceAtlas() {
  if (_refCache) return _refCache;
  const seq = assets.reference.sequence;
  const refAtlas = engine.customAtlas({
    name: "NL4-3", sequence: seq, alignment: alignToReference(seq, seq) });
  _refCache = { perSite: perSiteSupply(groupByPosition(refAtlas.nt), "rate") };
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
    "Indels are excluded: their predictor has no nucleotide context.");
  return exportProvenance(page.release, lines);
}

// Restore anything carried over from a previous run in this tab.
const carried = sessionStorage.getItem("mra.sequences");
if (carried && !fasta.value) fasta.value = carried;
