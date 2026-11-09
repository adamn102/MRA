// SPDX-License-Identifier: MIT
//
// sequence.js -- submit a sequence, get a custom atlas.
//
// The scientific framing is load-bearing here, not decoration. This view must never
// let a reader mistake a modelled rate for a measured one, so:
//   - the privacy note and the standing caveat sit at the end of the page, after the
//     results, which is where a reader meets them with the numbers already in view;
//   - every result row carries flags naming which components it rests on;
//   - the codon-identity component is named first, because it is the one the submitted
//     sequence genuinely determines;
//   - a sequence that aligns entirely outside the modelled region gets a refusal that
//     names the covered region, not an empty table.
//
// Up to MAX_DETAIL_SEQUENCES submissions share one genome track. Checkboxes choose which
// lines are drawn; one selector chooses which sequence the pie and the flow chart use.

import { bootstrapPage, privacyPanel, note, caveat, showError, exportProvenance,
         CUSTOM_CAVEAT, NEF_NOTE } from "./page.js";
import { loadSmallAssets } from "./data.js";
import { Analysis, progressRow, MAX_DETAIL_SEQUENCES } from "./analysis.js";
import { parseFasta, alignToReference } from "./align.js";
import { renderTable, formatRate, toCsv, downloadText } from "./table.js";
import { MraEngine, GENETIC_CODE } from "./mra_engine.js";
import { perSiteSupply, groupByPosition, buildFrameIndex, geneSegments, targetFor,
         codonDifferences, accessibilityRows, smoothSeries } from "./explore_core.js";
import { drawSupplyTrack, drawPositionPie, CLASS_LABEL } from "./track_view.js";
import { trackSliderGeometry, dotPlot } from "./chart.js";
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
document.getElementById("load-nl43").addEventListener("click", () => {
  fasta.value = ">NL4-3_reference\n" + assets.reference.sequence.replace(/(.{60})/g, "$1\n");
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
 * which sequences are drawn; the "Pie and flow" selector chooses which sequence those
 * two charts analyse. The plot is the first thing under Results.
 */
function renderDashboard(slot, items, allRows, truncated) {
  slot.innerHTML = "";
  const SEGMENTS = geneSegments(assets.reference.proteins);
  const shown = new Set(items.map((_, i) => i));
  let active = 0;
  let selected = null;
  let markerNt = nearestCodon(items[0], items[0].frameAt.has(1600) ? 1600 : (items[0].variants[0] || 1600));
  let flowToken = 0;
  let lowerFor = -1;

  const fig = document.createElement("figure");
  fig.className = "chart plot-card";
  fig.setAttribute("data-export-name", "genome_track");
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
  Object.assign(elSlider.style, trackSliderGeometry());
  sliderWrap.appendChild(elSlider);
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
  fig.append(trackSlot, sliderWrap, figCap, toggles);

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
    `<div><label>Pie and flow</label><select class="sq-which"></select></div>` +
    `<div class="bar-actions"><button class="sq-go">Go</button></div>`;

  const head = document.createElement("div");
  const dash = document.createElement("figure");
  dash.className = "dashboard-figures";
  dash.setAttribute("data-export-name", "position_dashboard_supply");
  const chartEl = document.createElement("div");
  const captionEl = document.createElement("figcaption");
  dash.append(chartEl, captionEl);
  const statsSlot = document.createElement("div");
  const lower = document.createElement("div");
  slot.append(fig, ctrl, head, dash, statsSlot, lower);

  const elWin = ctrl.querySelector(".sq-win"), elWinOut = ctrl.querySelector(".sq-win-out");
  const elClass = ctrl.querySelector(".sq-class"), elSys = ctrl.querySelector(".sq-sys");
  const elTarget = ctrl.querySelector(".sq-target"), elPos = ctrl.querySelector(".sq-pos");
  const elFrame = ctrl.querySelector(".sq-frame"), elWhich = ctrl.querySelector(".sq-which");
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
      statsEl: statsSlot,
      subjectName: shortName(item.row.name),
      refBase: item.baseAt.get(markerNt),
      onSelect: (pos) => selectNt(pos),
    });
    const where = selected ? ` (${engine.hxb2Label(selected.protein, selected.geneAa)})` : "";
    figCap.innerHTML =
      `Genome-wide ${CLASS_LABEL[cls]}, ${win} nt mean. Marker at ` +
      `${markerNt.toLocaleString()}${where}. Red ticks are where ` +
      `<strong>${escapeHtml(shortName(item.row.name))}</strong> differs from NL4-3. ` +
      `The dashed grey line is NL4-3.`;
    setDetails(fig, trackText, { title: "The genome supply track" });
    syncFigures();
  }

  function showPie(item, rows, label, noteHtml, yourCodon, protein, geneAa) {
    const token = ++flowToken;
    const opts = {
      custom: true, noteHtml, refCodon: yourCodon,
      markerPhrase: "the codon in " + escapeHtml(item.row.name),
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
    showPie(item, aaRows, label, noteHtml, yourCodon, protein, geneAa);
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
    if (lowerFor === active) return;
    lowerFor = active;
    lower.innerHTML = "";
    if (truncated) lower.appendChild(truncationNote());
    appendSummary(lower, allRows);
    renderSequenceTables(lower, items[active]);
    syncFigures();
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
    lowerFor = -1;
    const nt = nearestCodon(item, markerNt);
    if (nt === null) {
      showError(head, `${item.row.name} has no modelled codon to chart.`);
      chartEl.innerHTML = "";
      captionEl.textContent = "";
    } else {
      selectNt(nt);
    }
    renderLower();
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


// --- comparison to NL4-3 ----------------------------------------------------------
// The accessibility view needs the NL4-3 rate at each position. Rather than fetch the heavy
// published amino-acid layer (which would break this page's "small assets only" design), build
// it once from the reference's OWN custom atlas: the round-trip guarantee makes that identical
// to the published layer, and it is computed in-browser from assets already loaded. Cached so a
// batch of sequences pays the cost once.
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
  const byPosition = new Map();
  for (const r of refAtlas.aa) {
    const key = r.protein + "|" + r.position;
    let list = byPosition.get(key);
    if (!list) byPosition.set(key, (list = []));
    list.push(r);
  }
  // The nucleotide layer costs nothing extra here -- the atlas is already built -- and it is
  // what the genome track overlays as the NL4-3 comparison line.
  _refCache = { byPosition, perSite: perSiteSupply(groupByPosition(refAtlas.nt), "rate") };
  return _refCache;
}
function refByPosition() { return referenceAtlas().byPosition; }

/** "Where your codons differ from NL4-3" -- reachability reshaped by codon identity. */
function renderCodonDifferences(slot, detail, row, geom) {
  const rows = codonDifferences(detail.aa, assets.reference.proteins, GENETIC_CODE);
  const heading = document.createElement("h4");
  heading.textContent = "Where your codons differ from NL4-3";
  slot.appendChild(heading);
  if (!rows.length) {
    slot.appendChild(caveat("Every analysed codon matches NL4-3, so there is nothing to " +
      "compare. (The worked NL4-3 example lands here.)"));
    return;
  }
  const syn = rows.filter((r) => r.kind === "synonymous codon change").length;
  const reshaped = rows.filter((r) => r.kind === "synonymous codon change" && (r.n_gained || r.n_lost)).length;
  slot.appendChild(note(
    `<strong>${rows.length.toLocaleString()} analysed codons differ from NL4-3</strong>, ` +
    `${syn.toLocaleString()} of them encoding the <em>same amino acid</em> by a different codon. ` +
    `<strong>${reshaped.toLocaleString()} of those synonymous differences change which amino acids are ` +
    `reachable in one substitution</strong> \u2014 the protein is unchanged, the accessible next steps ` +
    `are not. Gained and lost come from the genetic code alone, not from any rate.`));
  const tableSlot = document.createElement("div");
  slot.appendChild(tableSlot);
  renderTable(tableSlot, {
    columns: [
      { key: "hxb2_label", label: "HXB2" },
      { key: "query_nt_start", label: "Your nt", type: "num" },
      { key: "ref_codon", label: "NL4-3 codon" },
      { key: "your_codon", label: "Your codon" },
      { key: "ref_aa", label: "NL4-3 aa" },
      { key: "your_aa", label: "Your aa" },
      { key: "kind", label: "Kind" },
      { key: "gained", label: "Newly reachable in 1 nt",
        title: "Amino acids reachable from your codon in one substitution but not from the NL4-3 codon" },
      { key: "lost", label: "No longer reachable in 1 nt" },
    ],
    rows,
    sort: { key: "n_gained", dir: "desc" },
    filterPlaceholder: "filter by protein, codon, amino acid\u2026",
    downloadName: `mra_${page.version}_${safeName(row.name)}_codon_differences.csv`,
    provenance: customProvenance(row, row.summary, geom),
  });
}

/** "What has become more or less accessible" -- log2 FC of each AA change vs NL4-3. */
function renderAccessibility(slot, detail, row, geom) {
  const rows = accessibilityRows(detail.aa, refByPosition(), assets.reference.proteins, GENETIC_CODE);
  const heading = document.createElement("h4");
  heading.textContent = "What has become more or less accessible";
  slot.appendChild(heading);
  const withFc = rows.filter((r) => r.log2fc !== null && Number.isFinite(r.log2fc));
  if (!withFc.length) {
    slot.appendChild(caveat("No amino-acid change at a covered position differs in accessibility " +
      "from NL4-3 \u2014 expected when the sequence is identical to the reference."));
    return;
  }
  const chartSlot = document.createElement("div");
  const fig = document.createElement("figure"); fig.className = "chart";
  fig.appendChild(chartSlot);
  const cap = document.createElement("figcaption");
  fig.appendChild(cap);
  slot.appendChild(fig);
  const top = withFc.slice().sort((a, b) => Math.abs(b.log2fc) - Math.abs(a.log2fc)).slice(0, 25);
  dotPlot(chartSlot,
    top.map((r) => ({ label: `${r.hxb2_label} ${r.aa_change}`, value: r.log2fc,
                      note: `${r.ref_codon}\u2192${r.your_codon}` })),
    { title: "The 25 amino-acid changes whose accessibility shifts most", width: 740 });
  const shifted = withFc.filter((r) => Math.abs(r.log2fc) > 1).length;
  const stepChanged = rows.filter((r) => r.change && r.change.includes("substitutions")).length;
  cap.innerHTML =
    `${withFc.length.toLocaleString()} changes are in both atlases. ${shifted.toLocaleString()} ` +
    `shift more than two-fold, and ${stepChanged.toLocaleString()} need a different <em>number</em> of ` +
    `substitutions than in NL4-3 \u2014 a categorical change, not a quantitative one. ` +
    `Codon at the right of each row.`;
  const tableSlot = document.createElement("div");
  slot.appendChild(tableSlot);
  renderTable(tableSlot, {
    columns: [
      { key: "hxb2_change", label: "Mutation (HXB2)" },
      { key: "ref_codon", label: "NL4-3 codon" },
      { key: "your_codon", label: "Your codon" },
      { key: "mu_nl43", label: "mu_aa NL4-3", type: "rate" },
      { key: "mu_custom", label: "mu_aa yours", type: "rate" },
      { key: "log2fc", label: "log2 FC", type: "num",
        render: (r) => r.log2fc === null ? "\u2014" : r.log2fc.toFixed(2) },
      { key: "min_nt_nl43", label: "min nt NL4-3", type: "num" },
      { key: "min_nt_custom", label: "min nt yours", type: "num" },
      { key: "change", label: "Note" },
    ],
    rows,
    sort: { key: "log2fc", dir: "desc" },
    filterPlaceholder: "filter\u2026",
    downloadName: `mra_${page.version}_${safeName(row.name)}_accessibility_vs_nl43.csv`,
    provenance: customProvenance(row, row.summary, geom),
  });
  tableSlot.appendChild(caveat(
    "A fold change reflects your codon and local context, not a different mutation rate in your " +
    "virus: the per-site model is the same one, fitted on NL4-3."));
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

const COMPONENT_FLAGS = (r) => {
  const flags = [];
  flags.push('<span class="flag codon" title="Determined entirely by your sequence\'s codon">codon</span>');
  if (r.context_source === "sequence") {
    flags.push('<span class="flag context" title="3\u2032 context read from your sequence">context</span>');
  } else if (r.context_source === "pooled") {
    flags.push('<span class="flag" title="The model fitted this context cell type-only (fewer than 20 WT counts), so no context delta applies">context pooled</span>');
  } else if (r.context_source === "unavailable") {
    flags.push('<span class="flag fallback" title="The 3\u2032 base ran off the end of your sequence or was ambiguous; the type-only rate was used">context missing</span>');
  }
  flags.push(r.position_source === "aligned"
    ? '<span class="flag position" title="Fitted NL4-3 position spline applied at the aligned position">position</span>'
    : '<span class="flag fallback" title="Could not be aligned; the genome mean of the spline was substituted">position fallback</span>');
  return flags.join(" ");
};

/** Tables for the sequence selected in "Pie and flow". They sit below the plots. */
function renderSequenceTables(slot, item) {
  const { detail, row, geom } = item;
  const s = row.summary;
  const heading = document.createElement("h3");
  heading.textContent = `Tables for ${row.name}`;
  slot.appendChild(heading);

  // A sequence far from NL4-3 still aligns and still returns numbers, so the only thing
  // standing between a reader and a misread is saying so. The threshold is empirical: across
  // the LANL group M compendium, subtypes and CRFs sit at 85-93% identity to NL4-3, while
  // groups N/O/P and the SIVcpz/SIVgor outgroups fall at 74-78%. Below 80% the position
  // effects measured on NL4-3 are being extrapolated further than this model supports.
  if (s.identity !== null && s.identity < 0.80) {
    slot.appendChild(caveat(
      `<strong>This sequence is only ${(s.identity * 100).toFixed(1)}% identical to NL4-3.</strong> ` +
      `Group M subtypes and CRFs sit at 85&ndash;93%; non-M groups and SIV outgroups fall at ` +
      `74&ndash;78%. The rates below are still computed, but position effects fitted on NL4-3 ` +
      `are being carried a long way from where they were measured &mdash; treat them as ` +
      `indicative, not quantitative. If you pasted an RNA, protein or non-HIV-1 sequence, that ` +
      `is the more likely explanation.`));
  }
  const spans = (geom && geom.unalignedSpans) || [];
  if (spans.length) {
    slot.appendChild(caveat(
      `${spans.length} span(s) of ${escapeHtml(row.name)} could not be aligned: ` +
      spans.slice(0, 6).map((x) => `${x.queryStart}\u2013${x.queryEnd}`).join(", ") +
      (spans.length > 6 ? ", \u2026" : "") + "."));
  }
  if (s.proteins && s.proteins.Nef && s.proteins.Nef.codonsAnalysed > 0) {
    slot.appendChild(caveat(NEF_NOTE));
  }
  if (s.contextUnavailable || s.gaContextUnavailable) {
    slot.appendChild(caveat(
      `${s.contextUnavailable} site(s) had no readable 3&prime; base and fell back to the ` +
      `type-only rate; ${s.gaContextUnavailable} G&rarr;A site(s) were omitted because their ` +
      `3&prime; 2-mer could not be read \u2014 a G&rarr;A rate is not defined by this model without it. ` +
      `This happens at the very end of a submitted sequence.`));
  }

  renderCodonDifferences(slot, detail, row, geom);
  renderAccessibility(slot, detail, row, geom);

  const aaHeading = document.createElement("h4");
  aaHeading.textContent = "Amino-acid mutability";
  slot.appendChild(aaHeading);
  const aaSlot = document.createElement("div");
  slot.appendChild(aaSlot);
  renderTable(aaSlot, {
    columns: [
      { key: "hxb2_change", label: "Mutation (HXB2)", title: "Clinical / publication numbering" },
      { key: "codon", label: "Your codon", title: "Read from your sequence. This is what differs between viruses." },
      { key: "query_nt_start", label: "Your nt", type: "num", title: "1-based position in the sequence you submitted" },
      { key: "query_codon_number", label: "Your codon #", type: "num" },
      { key: "mu_aa", label: "mu_aa", type: "rate" },
      { key: "interval", label: "95% interval" },
      { key: "min_nt_changes", label: "min nt", type: "num" },
      { key: "n_paths", label: "paths", type: "num" },
      { key: "aa_change", label: "NL4-3 aa change", title: "Same change in NL4-3 gene numbering" },
      { key: "position", label: "NL4-3 aa", type: "num", title: "Internal coordinate, for traceability" },
    ],
    rows: detail.aa.map((r) => ({ ...r,
      interval: `${formatRate(r.mu_aa_lo95)} \u2013 ${formatRate(r.mu_aa_hi95)}` })),
    sort: { key: "mu_aa", dir: "desc" },
    filterPlaceholder: "filter by protein, change, codon\u2026",
    downloadName: `mra_${page.version}_${safeName(row.name)}_aa.csv`,
    provenance: customProvenance(row, s, geom),
  });

  const ntHeading = document.createElement("h4");
  ntHeading.textContent = "Nucleotide rates";
  slot.appendChild(ntHeading);
  const ntSlot = document.createElement("div");
  slot.appendChild(ntSlot);
  renderTable(ntSlot, {
    columns: [
      { key: "query_position", label: "Your nt", type: "num" },
      { key: "ref_position", label: "NL4-3 nt", type: "num" },
      { key: "ref_nucleotide", label: "Base" },
      { key: "mut_type", label: "Substitution" },
      { key: "context_3p", label: "3\u2032 context", title: "Read from your sequence, not from NL4-3" },
      { key: "rate", label: "rate", type: "rate" },
      { key: "interval", label: "95% interval" },
      { key: "interval_source", label: "Interval from",
        title: "published = the atlas interval for this exact site and substitution; type-context = a mean over sites, an approximation" },
      { key: "components", label: "Rests on", render: COMPONENT_FLAGS },
    ],
    rows: detail.nt.map((r) => ({ ...r,
      interval: `${formatRate(r.lo95)} \u2013 ${formatRate(r.hi95)}` })),
    sort: { key: "rate", dir: "desc" },
    filterPlaceholder: "filter by substitution, context\u2026",
    downloadName: `mra_${page.version}_${safeName(row.name)}_nt.csv`,
    provenance: customProvenance(row, s, geom),
  });
  ntSlot.appendChild(caveat(
    "Insertions and deletions are deliberately absent. Their predictor has no nucleotide " +
    "context, so a \u201ccustom\u201d indel rate would just be the NL4-3 indel rate at that " +
    "coordinate and would tell you nothing about your virus. They are in the published " +
    "browse and download layer, where they are a real result."));

  const btn = document.createElement("button");
  btn.textContent = "Download this custom atlas (both tables + paths)";
  btn.addEventListener("click", () => {
    const prov = customProvenance(row, s, geom);
    downloadText(`mra_${page.version}_${safeName(row.name)}_paths.csv`,
      toCsv(detail.paths, [
        { key: "protein" }, { key: "position" }, { key: "codon" }, { key: "mutant_codon" },
        { key: "ref_aa" }, { key: "alt_aa" }, { key: "n_nt_changes" },
        { key: "mutation_type_path" }, { key: "mu_path" }, { key: "mu_path_lo95" },
        { key: "mu_path_hi95" },
      ], prov));
  });
  slot.appendChild(btn);
}


function safeName(name) {
  return String(name).replace(/[^A-Za-z0-9._-]+/g, "_").slice(0, 60);
}

/** Provenance for a custom export: which model, which version, and what it rests on. */
function customProvenance(row, s, geom) {
  return exportProvenance(page.release, [
    `Custom atlas for submitted sequence: ${row.name}`,
    `Retrieved from the MRA portal; computed entirely in the browser.`,
    `Alignment to NL4-3 ${geom.refStart}-${geom.refEnd}, identity ${(s.identity * 100).toFixed(3)}%, ` +
      `${geom.orientation}.`,
    `Positions dropped for lying outside the modelled region ${region.first_position}-${region.last_position}: ` +
      `${s.droppedOutsideRegion}.`,
    `Spline fallback (genome mean) used at ${s.fallbackPositions} position(s).`,
    "THIS IS NOT A MEASUREMENT OF THIS VIRUS. It is the MRA model, fitted on NL4-3,",
    "applied to this sequence. Codon identity and local 3' context come from the submitted",
    "sequence; position effects were measured on NL4-3 and transfer only where it aligns.",
    "Indels are excluded: their predictor has no nucleotide context.",
  ]);
}

// Restore anything carried over from a previous run in this tab.
const carried = sessionStorage.getItem("mra.sequences");
if (carried && !fasta.value) fasta.value = carried;
