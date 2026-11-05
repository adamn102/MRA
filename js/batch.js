// SPDX-License-Identifier: MIT
//
// batch.js -- the batch-analysis tab.
//
// Compiles one atlas across a batch of sequences and compares them codon by codon. The batch
// is run through Analysis.batch(), which folds each sequence's atlas into two per-sample
// matrices and discards the atlas itself: retaining 250 full-genome atlases would be ~715 MB,
// the two matrices together are ~22 MB. Per-sequence tables are therefore not available here --
// that is what "Your sequence" is for.
//
// WHAT THIS TAB IS ACTUALLY MEASURING. The fitted model is identical for every submitted
// sequence, so a batch is not a set of different mutation rates. Two samples differ at a
// position only when they carry a different codon, or a different base immediately 3' of it;
// otherwise their rates are bit-for-bit equal. Everything presented here is therefore a
// statement about codon usage in the submitted set, which is why the views are keyed on the
// codon matrix (AA_STATS) rather than the nucleotide one.

import { bootstrapPage, privacyPanel, note, caveat, showError, showFatal, exportProvenance,
         shareButton, CUSTOM_CAVEAT } from "./page.js";
import { loadSmallAssets } from "./data.js";
import { Analysis, progressRow, MAX_BATCH_SEQUENCES } from "./analysis.js";
import { parseFasta } from "./align.js";
import { renderTable, formatRate, toCsv, downloadText } from "./table.js";
import { MraEngine } from "./mra_engine.js";
import { batchSiteStats, topDivergentSites, perSampleSummary,
         buildAaIndex, aaPositionStats, weightedDestinationSupply, codonSensitivity,
         samplesAtAaPosition, siteSupplyStats, perSampleRates, MUT_TYPES }
  from "./batch_core.js";
import { drawMutabilityTrack, drawMutabilityDivergenceTrack,
         drawPositionComposition, drawSupplyFlow } from "./batch_view.js";
import { geneSegments, targetFor as targetForIn } from "./explore_core.js";

const page = await bootstrapPage("batch.html");
document.getElementById("privacy-slot").appendChild(privacyPanel());
document.getElementById("caveat-slot").appendChild(caveat(CUSTOM_CAVEAT));
// The expected wait is stated up front: a full batch of genomes is most of a minute, and a
// user who does not know that will assume the page has hung.
document.getElementById("batch-note").textContent =
  `Up to ${MAX_BATCH_SEQUENCES} sequences per batch. A full batch takes roughly 10 seconds ` +
  `for pol amplicons and up to a minute for whole genomes; progress is shown as it runs.`;

let assets, engine, analysis, region;
try {
  assets = await loadSmallAssets(page.version);
  engine = new MraEngine(assets);
  region = assets.reference.modelled_region;
  analysis = new Analysis(assets);
  await analysis.start();
} catch (err) {
  showFatal("Could not load the atlas tables: " + err.message);
  throw err;
}

// Declared before anything below uses it: module code runs top to bottom, and a `const`
// referenced above its declaration throws on the temporal dead zone, which takes the
// whole module -- and so every control on the page -- down with it.
const el = (id) => document.getElementById(id);
const fasta = el("fasta");

const SEGMENTS = geneSegments(assets.reference.proteins);

// Which reading frames code at each nucleotide. buildFrameIndex() in explore_core filters by
// the codons one sequence actually covered; a batch has no single such set (the per-sequence
// atlases are discarded by design), so this indexes every codon in the reference instead.
// Protein / subunit targets for the coordinate box. Browse builds these from the amino-acid
// layer; this page does not load that layer, so the spans come from the reference itself --
// the same proteins, just sourced from what is already in memory.
const PROTEIN_SPANS = (() => {
  const m = new Map();
  for (const [name, rec] of Object.entries(assets.reference.proteins)) {
    if (!rec.codon_nt || !rec.codon_nt.length) continue;
    m.set(name, { min: rec.first_aa, max: rec.first_aa + rec.codon_nt.length - 1 });
  }
  return m;
})();
const TARGETS = engine.proteinTargets(PROTEIN_SPANS);
const TARGET_BY_KEY = new Map(TARGETS.map((t) => [t.key, t]));
{
  const sel = el("protein");
  for (const t of TARGETS) {
    const o = document.createElement("option");
    o.value = t.key; o.textContent = t.label;
    sel.appendChild(o);
  }
  sel.value = TARGETS.some((t) => t.key === "PR") ? "PR" : TARGETS[0].key;
}

/**
 * The nucleotide a protein-position lookup should land on.
 *
 * A codon spans three nucleotides and the track is per nucleotide, so a protein position alone
 * Returns null with a reason when the typed
 * position has no NL4-3 equivalent, which happens where HXB2 has residues NL4-3 lacks.
 */
function ntForProteinPosition(targetKey, ownPosition) {
  const target = TARGET_BY_KEY.get(targetKey);
  if (!target) return { nt: null, reason: "unknown protein" };
  const geneAa = engine.geneAaFromOwnPosition(target, ownPosition);
  if (geneAa === null) {
    return { nt: null, reason: `HXB2 ${targetKey} ${ownPosition} has no NL4-3 equivalent ` +
                               `(NL4-3 lacks that residue)` };
  }
  const rec = assets.reference.proteins[target.protein];
  const codon = rec && rec.codon_nt ? rec.codon_nt[geneAa - rec.first_aa] : null;
  if (!codon) {
    return { nt: null, reason: `${targetKey} ${ownPosition} is outside the modelled region` };
  }
  // The codon's first base. Which of the three is immaterial now that everything on this
  // page is keyed on the codon rather than the nucleotide; it only has to land inside it.
  return { nt: codon[0], geneAa, protein: target.protein, target };
}

const FRAME_AT = (() => {
  const m = new Map();
  for (const [protein, rec] of Object.entries(assets.reference.proteins)) {
    if (!rec.codon_nt) continue;
    rec.codon_nt.forEach((nts, i) => {
      for (const nt of nts) {
        let list = m.get(nt);
        if (!list) m.set(nt, (list = []));
        list.push({ protein, geneAa: rec.first_aa + i });
      }
    });
  }
  return m;
})();

// --- batch state, replaced wholesale by each run -----------------------------------
let BATCH = null;      // {matrix, names, summaries, geometry}
let STATS = [];         // batchSiteStats for the selected class (nucleotide; rate summary)
let AA_STATS = [];      // aaPositionStats -- what the tracks, scatter and panel read
let AA_INDEX_REF = null;
let markerCol = null;   // index into AA_STATS of the selected codon
let markerNt = 1600;
let markerSetByUser = false;   // once a position is chosen, a class change must not move it
let field = "nonsyn";
const CLASS_LABEL = { nonsyn: "all residue changes", ga: "G\u2192A only",
                      nonga: "non-G\u2192A only" };
// The nucleotide matrix carries only total and G->A, so it cannot express the non-G->A class;
// the per-sample coverage table, which is the only thing still reading it, falls back to total
// there. The codon views carry the exact split and are unaffected.
const ntFieldFor = (f) => (f === "ga" ? "ga" : "total");


// --- input helpers -------------------------------------------------------------------
el("load-file").addEventListener("click", () => el("file-input").click());
el("file-input").addEventListener("change", async (ev) => {
  const f = ev.target.files && ev.target.files[0];
  if (f) fasta.value = await f.text();
  ev.target.value = "";
});
el("clear").addEventListener("click", () => {
  fasta.value = "";
  el("results").hidden = true;
  el("run-status").innerHTML = "";
  el("progress-slot").innerHTML = "";
  BATCH = null; STATS = []; AA_STATS = []; markerCol = null; markerSetByUser = false;
  try { sessionStorage.removeItem("mra.batch"); } catch { /* nothing to clear */ }
});

// A worked example: the reference plus copies carrying a few codon changes, so the tab
// demonstrates its own point (codon identity moves supply) without needing real patient data.
// --- the worked example ---------------------------------------------------------------
// The LANL compendium rather than anything generated here: a real set is the only way to see
// what this page is for, because the interesting quantity -- how codon usage spreads
// mutational accessibility -- does not exist in sequences we invent. Fetched on demand rather
// than precached, since it is ~1.7 MB and most visits will not want it.
el("load-demo").addEventListener("click", async () => {
  const btn = el("load-demo");
  const prev = btn.textContent;
  btn.disabled = true;
  btn.textContent = "Loading\u2026";
  el("run-status").innerHTML = "";
  try {
    const resp = await fetch(`data/${page.version}/lanl_batch_demo.fasta`);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    const text = await resp.text();
    fasta.value = text;
    const n = (text.match(/^>/gm) || []).length;
    el("run-status").appendChild(note(
      `<strong>${n} HIV-1 group M genomes</strong> from the Los Alamos HIV Sequence Database ` +
      `2025 compendium \u2014 one reference strain per subtype and circulating recombinant ` +
      `form, named <code>subtype.country.year.name.accession</code>. Alignment gaps have been ` +
      `stripped, because each sequence is aligned to NL4-3 here independently. ` +
      `Press <em>Compile batch atlas</em>; whole genomes take roughly half a minute.`));
    el("run-status").appendChild(caveat(
      `Groups N, O and P and the SIVcpz/SIVgor outgroups in the published compendium are ` +
      `left out of this set. They sit at 74\u201378% identity to NL4-3 against 85\u201393% ` +
      `for group M, so their codons differ because the lineages diverged long before group M ` +
      `existed \u2014 two of them are not HIV-1 at all. Including them would fill the ` +
      `divergence track with between-lineage distance under a label that says codon usage. ` +
      `Even within group M, remember these are one reference strain per subtype, so the ` +
      `diversity here is between-subtype, not within-host.`));
  } catch (err) {
    showError(el("run-status"),
              "Could not load the compendium: " + err.message +
              ". It ships with the data release; check that data/" + page.version +
              "/lanl_batch_demo.fasta was deployed.");
  } finally {
    btn.disabled = false;
    btn.textContent = prev;
  }
});

// --- the run ---------------------------------------------------------------------------
el("run").addEventListener("click", async () => {
  const text = fasta.value.trim();
  el("run-status").innerHTML = "";
  if (!text) { showError(el("run-status"), "Paste a set of sequences first."); return; }
  const records = parseFasta(text).filter((r) => r.sequence.length >= 30);
  if (records.length < 2) {
    showError(el("run-status"), "A batch needs at least two sequences of 30 bases or more. " +
      "For a single sequence use the \u201cYour sequence\u201d tab, which gives its full tables.");
    return;
  }
  try { sessionStorage.setItem("mra.batch", text); } catch { /* quota: not essential */ }

  const bar = progressRow(el("progress-slot"));
  const t0 = performance.now();
  let res;
  try {
    res = await analysis.batch(records, (p) => bar.update(p));
  } catch (err) {
    showError(el("run-status"), "Batch analysis failed: " + err.message);
    return;
  }
  const ms = Math.round(performance.now() - t0);
  const analysed = res.summaries.filter((s) => s && !s.unalignable).length;
  bar.finish(`Done: ${analysed} of ${res.summaries.length} sequence(s) in ${ms} ms.`);
  BATCH = res;
  if (res.truncated) {
    el("run-status").appendChild(note(
      `Only the first ${MAX_BATCH_SEQUENCES} sequences were compiled. The cap is set from ` +
      `measured timing: a batch of ${MAX_BATCH_SEQUENCES} clinical pol amplicons takes about ` +
      `9 seconds and the same number of full-length genomes about 36 seconds.`));
  }
  const bad = res.summaries.map((s, i) => [s, i]).filter(([s]) => s && s.unalignable);
  if (bad.length) {
    el("run-status").appendChild(caveat(
      `<strong>${bad.length} sequence(s) could not be aligned</strong> and are excluded from ` +
      `every statistic below: ` +
      bad.map(([s]) => `${s.name} (${s.reason || "no anchor found"})`).join("; ") +
      `. They are not counted as zero.`));
  }
  if (!analysed) {
    showError(el("run-status"), "No sequence in this set could be aligned to the reference.");
    return;
  }
  el("results").hidden = false;
  recompute();
});

/** Rebuild the per-site statistics for the current class, then redraw everything. */
function recompute() {
  if (!BATCH) return;
  STATS = batchSiteStats(BATCH.matrix, { field: ntFieldFor(field) });
  // Mutability is a property of a codon, so everything the user looks at is keyed here.
  if (!AA_INDEX_REF) AA_INDEX_REF = buildAaIndex(assets.reference.proteins);
  AA_STATS = BATCH.aaMatrix ? aaPositionStats(BATCH.aaMatrix, AA_INDEX_REF, { field }) : [];
  const top = topDivergentSites(AA_STATS, { limit: 1 })[0];
  renderSummary();
  // Open on the most divergent site: with nothing selected yet, that is the position most
  // worth looking at, and it is what makes the batch differ from a single sequence.
  if (top && !markerSetByUser) markerNt = top.nt_start;
  if (!STATS.some((r) => r.ref_position === markerNt)) {
    markerNt = STATS.length ? STATS[Math.floor(STATS.length / 2)].ref_position : region.first_position;
  }
  el("track-slider").min = String(region.first_position);
  el("track-slider").max = String(region.last_position);
  el("track-slider").value = String(markerNt);
  el("position").value = String(markerNt);
  drawAll();
  renderBatchTable();
  renderSampleTable();
}

function renderSummary() {
  const analysed = BATCH.summaries.filter((s) => s && !s.unalignable).length;
  const diverging = AA_STATS.filter((r) => r.divergence !== null && r.divergence > 0).length;
  const top = topDivergentSites(AA_STATS, { limit: 1 })[0];
  const slot = el("batch-summary");
  const resVar = AA_STATS.filter((r) => r.n_aa > 1).length;
  const synVar = AA_STATS.filter((r) => r.n_aa === 1 && r.n_codons > 1).length;
  slot.innerHTML = "";
  slot.appendChild(note(
    `<strong>${analysed} sequences</strong> compiled over ` +
    `<strong>${AA_STATS.length.toLocaleString()}</strong> codons. ` +
    `Mutability differs between samples at <strong>${diverging.toLocaleString()}</strong> of ` +
    `them (${(100 * diverging / Math.max(1, AA_STATS.length)).toFixed(1)}%)` +
    (top ? `, most at <strong>${engine.hxb2Label(top.protein, top.position)}</strong>, where it ` +
           `spans ${Math.pow(2, top.divergence).toFixed(1)}-fold across ${top.n_covered} samples`
         : "") +
    `. <strong>${resVar.toLocaleString()}</strong> position${resVar === 1 ? "" : "s"} ` +
    `differ${resVar === 1 ? "s" : ""} in residue and a further ` +
    `<strong>${synVar.toLocaleString()}</strong> carry the same residue on different codons, ` +
    `where the protein is identical but the supply is not. The fitted model is the same for ` +
    `every sample \u2014 what differs is the codon each one carries.`));
}

/** Find the codon stat covering a nucleotide, so a track click lands on a position. */
function aaStatForNt(nt) {
  if (markerCol !== null && AA_STATS[markerCol] &&
      AA_STATS[markerCol].nt_start <= nt && nt < AA_STATS[markerCol].nt_start + 3) {
    return AA_STATS[markerCol];
  }
  let best = null, bestD = Infinity;
  for (let i = 0; i < AA_STATS.length; i++) {
    const d = Math.abs(AA_STATS[i].nt_start - nt);
    if (d < bestD) { bestD = d; best = i; }
  }
  markerCol = best;
  return best === null ? null : AA_STATS[best];
}

function drawAll() {
  const win = Number(el("track-window").value);
  const stat = aaStatForNt(markerNt);
  const label = stat ? engine.hxb2Label(stat.protein, stat.position) : `nt ${markerNt}`;
  const nSamples = BATCH.summaries.filter((s) => s && !s.unalignable).length;

  drawMutabilityTrack(el("batch-track"), {
    stats: AA_STATS, segments: SEGMENTS, first: region.first_position,
    last: region.last_position, win, marker: stat ? stat.nt_start : markerNt,
    markerLabel: label, nSamples, onSelect: (p) => selectPosition(p),
  });
  const med = AA_STATS.length
    ? AA_STATS.map((r) => r.median).sort((a, b) => a - b)[Math.floor(AA_STATS.length / 2)] : 0;
  const nVar = AA_STATS.filter((r) => r.n_aa > 1).length;
  const nSyn = AA_STATS.filter((r) => r.n_aa === 1 && r.n_codons > 1).length;
  el("track-caption").innerHTML =
    `One point per codon: the summed rate of one-step changes that alter the residue ` +
    `(<strong>${CLASS_LABEL[field] || field}</strong>), smoothed over ${win} codons. The band is the full min\u2013max across samples, not an ` +
    `interquartile range \u2014 at most positions every sample carries the same codon, so the ` +
    `quartiles coincide and only the extremes separate. Median across codons ` +
    `${formatRate(med)}. <strong>${nVar.toLocaleString()}</strong> positions vary in residue; ` +
    `a further <strong>${nSyn.toLocaleString()}</strong> carry the same residue on different ` +
    `codons, where supply differs but the protein does not. Click to select a position.`;

  const dwin = Number(el("divergence-window").value);
  drawMutabilityDivergenceTrack(el("divergence-track"), {
    stats: AA_STATS, first: region.first_position, last: region.last_position,
    win: dwin, marker: stat ? stat.nt_start : markerNt,
    onSelect: (p) => selectPosition(p),
  });
  el("divergence-caption").innerHTML =
    `Per-codon log<sub>2</sub>(max/min) of mutability across the samples covering it. Each ` +
    `stem is one codon where the samples would generate residue change at different rates ` +
    `\u2014 which happens whenever the codon or its 3\u2032 context differs, including when ` +
    `the residue is identical. ` +
    (dwin > 1
      ? `The line is the mean over ${dwin} codons, counting non-diverging codons as zero, so ` +
        `it shows how concentrated disagreement is rather than how large any single one is. `
      : `No smoothing is applied at a window of 1. `) +
    `Positions covered by one sample show nothing, because a fold range needs two.`;

  drawSite(label);
}

function drawSite(label) {
  const stat = aaStatForNt(markerNt);
  const head = el("site-head");
  head.innerHTML = "";
  if (!stat) {
    showError(head, "That position is not inside a modelled codon.");
    el("site-chart").innerHTML = "";
    el("codon-chart").innerHTML = "";
    el("site-sample-table").innerHTML = "";
    return;
  }
  const refCodon = refCodonFor(stat);
  const sens = codonSensitivity(BATCH.aaMatrix, stat);
  const sstats = siteSupplyStats(BATCH.aaMatrix, stat);
  head.appendChild(note(
    `<strong>${label}</strong> &middot; NL4-3 ${stat.protein} ${stat.position} ` +
    `&middot; codon starts at nt ${stat.nt_start.toLocaleString()}` +
    (refCodon ? ` &middot; NL4-3 codon <code>${refCodon}</code>` : "") +
    `. Batch median mutability ${formatRate(stat.median)}` +
    (sens.fold_range
      ? `; across the ${sens.n_codons} codons present it spans ` +
        `<strong>${sens.fold_range.toFixed(1)}-fold</strong>, with ` +
        `${(100 * sens.fraction_on_fastest).toFixed(0)}% of samples on the highest-supply codon.`
      : `; every covering sample carries the same codon, so there is no codon sensitivity ` +
        `to report here.`) +
    (sstats
      ? ` <strong>${(100 * sstats.syn_fraction).toFixed(0)}%</strong> of the supply here is ` +
        `synonymous; the residue-changing part is spread over the equivalent of ` +
        `<strong>${sstats.effective_destinations.toFixed(1)}</strong> of its ` +
        `${sstats.n_destinations} reachable residues.`
      : "")));

  // What the batch carries, as a single strip. Only the residue row: the codon composition is
  // the flow chart's own left column, and showing it twice is what made this section hard to
  // read in the first place.
  drawPositionComposition(el("site-chart"), el("site-caption"), stat,
                          { only: "aa", title: `What the batch carries at ${label}` });
  // Where the supply goes, as a flow. Band width is codon frequency times that codon's rate,
  // so the left column is each codon's share of SUPPLY rather than of sequences -- the gap
  // between those two numbers is the point of compiling a batch at all.
  drawSupplyFlow(el("codon-chart"), el("codon-caption"), stat, BATCH.aaMatrix,
                 { refCodon, weight: el("flow-weight").value, stats: sstats,
                   title: `Where mutation at ${label} goes` });
  renderSiteSampleTable(stat, label, refCodon);

  // Keep the coordinate box in step with a track click.
  const target = targetForIn(TARGETS, stat.protein, stat.position);
  if (target) {
    el("protein").value = target.key;
    if (el("coord-system").value === "protein") {
      el("position").value = String(engine.ownPositionFor(target, stat.position));
    }
  }
  if (el("coord-system").value === "nt") el("position").value = String(markerNt);
  el("coord-note").textContent =
    `${stat.n_covered} of ${BATCH.names.length} submitted sequences cover this codon.`;
}

/** The reference's own codon at a position, so the batch's codons can be set against it. */
function refCodonFor(stat) {
  const rec = assets.reference.proteins[stat.protein];
  if (!rec || !rec.codon_nt) return null;
  const nts = rec.codon_nt[stat.position - rec.first_aa];
  if (!nts) return null;
  return nts.map((n) => assets.reference.sequence[n - 1]).join("");
}

/**
 * Every sample at the selected codon. The composition bars show proportions and the supply
 * chart shows destinations; neither carries sample names, so they live here.
 */
function renderSiteSampleTable(stat, label, refCodon) {
  const rows = samplesAtAaPosition(BATCH.aaMatrix, stat, BATCH.names, { field: "nonsyn" });
  renderTable(el("site-sample-table"), {
    columns: [
      { key: "sample", label: "Sequence" },
      { key: "aa", label: "Residue",
        render: (r) => r.aa ? (r.aa === stat.majority_aa ? r.aa : r.aa + " \u25c0") : "\u2014" },
      { key: "codon", label: "Codon",
        render: (r) => !r.codon ? "\u2014"
          : `${r.codon} (${r.aa || "?"})` + (r.codon === refCodon ? " NL4-3" : "") },
      { key: "value", label: "Mutability", type: "rate" },
      { key: "ratio", label: "\u00d7 batch median", type: "num",
        render: (r) => r.ratio === null ? "\u2014" : r.ratio.toFixed(2) },
      { key: "status", label: "Status" },
    ],
    rows: rows.map((r) => ({
      ...r,
      status: !r.covered ? "does not cover this codon"
            : r.aa !== stat.majority_aa ? "minority residue"
            : r.codon !== stat.majority_codon ? "synonymous codon variant" : "",
    })),
    sort: { key: "value", dir: "desc" },
    filterPlaceholder: "filter by sequence\u2026",
    downloadName: `mra_${page.version}_batch_${String(label).replace(/[^A-Za-z0-9]+/g, "_")}` +
                  `_by_sample.csv`,
    provenance: batchProvenance(),
  });
}

function selectPosition(pos) {
  markerSetByUser = true;
  markerNt = Math.min(Math.max(Math.round(pos), region.first_position), region.last_position);
  el("track-slider").value = String(markerNt);
  el("position").value = String(markerNt);
  drawAll();
}

// --- controls --------------------------------------------------------------------------
el("track-window").addEventListener("input", () => {
  el("track-window-out").textContent = el("track-window").value;
  if (BATCH) drawAll();
});
function syncCoordControls() {
  const mode = el("coord-system").value;
  el("protein").disabled = mode !== "protein";
  el("position").disabled = mode === "divergent";
}
el("coord-system").addEventListener("change", syncCoordControls);
syncCoordControls();

// Re-draw only the position panel: the weighting changes how the flow is apportioned, not
// anything the tracks or tables show.
el("flow-weight").addEventListener("change", () => {
  if (!BATCH) return;
  const st = aaStatForNt(markerNt);
  drawSite(st ? engine.hxb2Label(st.protein, st.position) : `nt ${markerNt}`);
});
el("divergence-window").addEventListener("input", () => {
  el("divergence-window-out").textContent = el("divergence-window").value;
  if (BATCH) drawAll();
});
el("track-class").addEventListener("change", () => {
  field = el("track-class").value;
  recompute();
});
el("track-slider").addEventListener("input", () => selectPosition(Number(el("track-slider").value)));
el("go").addEventListener("click", () => {
  if (!BATCH) return;
  const mode = el("coord-system").value;
  if (mode === "divergent") {
    const top = topDivergentSites(AA_STATS, { limit: 1 })[0];
    if (top) selectPosition(top.nt_start);
    else showError(el("site-head"), "No site in this batch has any inter-sample divergence: " +
                   "every sample carries the same codon everywhere they overlap.");
    return;
  }
  const v = Number(el("position").value);
  if (!Number.isFinite(v) || v < 1) {
    showError(el("site-head"), "Enter a position.");
    return;
  }
  if (mode === "nt") { selectPosition(v); return; }
  const r = ntForProteinPosition(el("protein").value, v);
  if (r.nt === null) { showError(el("site-head"), r.reason); return; }
  selectPosition(r.nt);
});
el("download-png").addEventListener("click", async () => {
  const button = el("download-png");
  if (typeof html2canvas === "undefined") {
    showError(el("site-head"), "The PNG library did not load; reload the page and try again.");
    return;
  }
  const target = el("site-export");
  button.disabled = true;
  const label = button.textContent;
  button.textContent = "Generating\u2026";
  try {
    const canvas = await html2canvas(target, { backgroundColor: "#ffffff", scale: 2,
                                               logging: false, useCORS: true });
    const blob = await new Promise((res) => canvas.toBlob(res));
    if (!blob) throw new Error("the canvas could not be converted to an image");
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `mra_batch_nt${markerNt}.png`;
    a.click();
    URL.revokeObjectURL(url);
  } catch (err) {
    showError(el("site-head"), "Could not build the PNG: " + err.message);
  } finally {
    button.disabled = false;
    button.textContent = label;
  }
});

// --- tables ------------------------------------------------------------------------------
function batchProvenance() {
  const analysed = BATCH.summaries.filter((s) => s && !s.unalignable);
  return exportProvenance(page.release, [
    `Compiled batch atlas over ${analysed.length} sequence(s).`,
    `Sequences: ${analysed.map((s) => s.name).join(", ").slice(0, 600)}`,
    `Supply class: ${field === "ga" ? "G->A only" : "all substitution classes"}.`,
    `Computed entirely in the browser; no sequence left this machine.`,
    "THIS IS NOT A MEASUREMENT OF THESE VIRUSES. It is the MRA model, fitted on NL4-3,",
    "applied to each submitted sequence. Rates differ between samples only through codon",
    "identity and local 3' context; position effects were measured on NL4-3.",
    "Sites a sample does not cover are excluded from that site's statistics, not counted as 0.",
    "Divergence is log2(max/min) across covering samples; null where a site's minimum is 0",
    "or fewer than two samples cover it.",
  ]);
}

/** Per-site rows, shared by the on-page table and the combined export. */
/**
 * The compiled atlas, one row per codon. Shared by the on-page table and the CSV export so
 * the two cannot drift. `mutability` is the one-step rate of residue change; synonymous and
 * stop supply are carried separately because they are different quantities, not components
 * of the same one.
 */
function siteRows() {
  return AA_STATS.map((r) => {
    const sens = codonSensitivity(BATCH.aaMatrix, r);
    const ss = siteSupplyStats(BATCH.aaMatrix, r);
    return {
      syn_fraction: ss ? ss.syn_fraction : null,
      stop_fraction: ss ? ss.stop_fraction : null,
      n_destinations: ss ? ss.n_destinations : null,
      effective_destinations: ss ? ss.effective_destinations : null,
      evenness: ss && ss.evenness !== null ? ss.evenness : null,
      top_destination: ss ? ss.top_destination : "",
      top_share: ss ? ss.top_share : null,
      protein: r.protein,
      hxb2_label: engine.hxb2Label(r.protein, r.position),
      nl43_position: r.position,
      nt_start: r.nt_start,
      n_covered: r.n_covered,
      majority_aa: r.majority_aa,
      majority_aa_freq: r.majority_aa_freq,
      n_aa: r.n_aa,
      entropy: r.entropy,
      majority_codon: r.majority_codon,
      majority_codon_freq: r.majority_codon_freq,
      n_codons: r.n_codons,
      mutability_median: r.median,
      mutability_min: r.min,
      mutability_max: r.max,
      codon_fold_range: sens.fold_range,
      fraction_on_fastest_codon: sens.fraction_on_fastest,
      divergence: r.divergence,
    };
  });
}

function renderBatchTable() {
  renderTable(el("batch-table"), {
    columns: [
      { key: "hxb2_label", label: "Position (HXB2)" },
      { key: "n_covered", label: "Samples", type: "num" },
      { key: "majority_aa", label: "Residue" },
      { key: "majority_aa_freq", label: "Residue freq", type: "num",
        render: (r) => (100 * r.majority_aa_freq).toFixed(0) + "%" },
      { key: "n_aa", label: "Residues", type: "num" },
      { key: "entropy", label: "Entropy (bits)", type: "num",
        render: (r) => r.entropy.toFixed(2) },
      { key: "majority_codon", label: "Codon" },
      { key: "n_codons", label: "Codons", type: "num" },
      { key: "mutability_median", label: "Mutability", type: "rate" },
      { key: "codon_fold_range", label: "Codon fold range", type: "num",
        render: (r) => r.codon_fold_range === null ? "\u2014"
                                                   : r.codon_fold_range.toFixed(2) },
      { key: "syn_fraction", label: "syn share", type: "num",
        title: "Fraction of this site's supply that does not change the residue",
        render: (r) => r.syn_fraction === null ? "\u2014"
                     : (100 * r.syn_fraction).toFixed(0) + "%" },
      { key: "effective_destinations", label: "eff. residues", type: "num",
        title: "exp(entropy) of the residue-changing supply \u2014 the effective number of " +
               "destination residues. Near 1 when one change takes almost all of it; near " +
               "the reachable count when the supply is spread evenly.",
        render: (r) => r.effective_destinations === null ? "\u2014"
                     : r.effective_destinations.toFixed(1) },
      { key: "divergence", label: "log2(max/min)", type: "num",
        render: (r) => r.divergence === null ? "\u2014" : r.divergence.toFixed(3) },
    ],
    rows: siteRows(),
    sort: { key: "divergence", dir: "desc" },
    filterPlaceholder: "filter by position, residue or codon\u2026",
    downloadName: `mra_${page.version}_batch_codon_atlas.csv`,
    provenance: batchProvenance(),
  });
}

/** Per-sample rows, shared by the on-page table and the combined export. */
function sampleRows() {
  const summ = perSampleSummary(BATCH.matrix, STATS, BATCH.names,
                                { field: ntFieldFor(field) });
  return summ.map((r, i) => {
    const s = BATCH.summaries[i] || {};
    const g = BATCH.geometry[i] || {};
    return {
      sample: r.sample,
      excluded: s.unalignable ? "yes" : "",
      sites_covered: r.sites_covered,
      fraction_of_region: r.fraction_of_region,
      identity: s.identity != null ? s.identity : null,
      orientation: g.orientation || "",
      ratio_to_batch: r.ratio_to_batch,
    };
  });
}

function renderSampleTable() {
  const rows = sampleRows();
  renderTable(el("sample-table"), {
    columns: [
      { key: "sample", label: "Sample" },
      { key: "sites_covered", label: "Sites covered", type: "num" },
      { key: "fraction_of_region", label: "% of region", type: "num",
        render: (r) => (100 * r.fraction_of_region).toFixed(1) },
      { key: "identity", label: "Identity to NL4-3", type: "num",
        render: (r) => r.identity === null ? "\u2014" : (100 * r.identity).toFixed(2) + "%" },
      { key: "orientation", label: "Orientation" },
      { key: "ratio_to_batch", label: "Rate vs batch median", type: "num",
        render: (r) => r.ratio_to_batch === null ? "\u2014" : r.ratio_to_batch.toFixed(3),
        title: "Median across sites of this sample's rate divided by the batch median" },
      { key: "excluded", label: "Excluded" },
    ],
    rows,
    sort: { key: "sites_covered", dir: "desc" },
    filterPlaceholder: "filter by sample\u2026",
    downloadName: `mra_${page.version}_batch_samples.csv`,
    provenance: batchProvenance(),
  });
}

/** Per-sequence rate summary: bulk, each substitution class, and G->A called out. */
/**
 * Per-sequence rate summary. No longer shown on the page -- it is a per-nucleotide-class
 * quantity, which sits oddly next to views keyed on codons -- but it is still written to the
 * combined CSV, where a reader wanting bulk and per-type rates can find it.
 */
function ratesRows() {
  return perSampleRates(BATCH.matrix, BATCH.names).map((r, i) => {
    const sm = BATCH.summaries[i] || {};
    const row = { sample: r.sample, sites_covered: r.sites_covered,
                  excluded: sm.unalignable ? "could not align" : "",
                  bulk_rate: r.bulk_rate, ga_rate: r.ga_rate, ga_share: r.ga_share,
                  transition_rate: r.transition_rate, transversion_rate: r.transversion_rate,
                  expected_per_cycle: r.bulk_rate === null ? null
                                    : r.bulk_rate * r.sites_covered };
    for (const t of MUT_TYPES) row[t] = r[t];
    return row;
  });
}


// --- the whole analysis as one CSV ------------------------------------------------------
// Three tables of different shapes cannot share a header row, so they are written as labelled
// blocks in one file rather than joined. Each block repeats the provenance so a block pasted
// out on its own still says what it is and what it is not.
/**
 * Every destination residue at every covered position, batch-weighted.
 *
 * Capped by default: at 2,984 positions and ~6 one-step destinations each this is ~18,000
 * rows, which is fine for a CSV but would be slow to build on every redraw, so it is only
 * built when the export runs.
 */
function destinationRows() {
  const out = [];
  for (const st of AA_STATS) {
    const label = engine.hxb2Label(st.protein, st.position);
    for (const d of weightedDestinationSupply(BATCH.aaMatrix, st)) {
      if (!(d.weighted_rate > 0)) continue;
      out.push({
        protein: st.protein, hxb2_label: label, nl43_position: st.position,
        ref_aa: st.majority_aa, alt_aa: d.alt_aa, n_covered: st.n_covered,
        weighted_rate: d.weighted_rate, rate_min: d.rate_min, rate_max: d.rate_max,
        fold_range_across_codons: d.fold_range,
        fraction_reaching: d.fraction_reaching,
        is_stop: d.is_stop, is_synonymous: d.is_synonymous,
      });
    }
  }
  return out;
}

el("download-all").addEventListener("click", () => {
  const prov = batchProvenance();
  const parts = [];
  const block = (title, rows, columns) => {
    parts.push(`# ===== ${title} =====`);
    parts.push(toCsv(rows, columns, prov));
    parts.push("");
  };
  // toCsv() writes each column's `key` as the header, so these are key-only specs: a `label`
  // here would look like a rename and do nothing.
  const cols = (...keys) => keys.map((k) => ({ key: k }));
  block("PER-SEQUENCE RATE SUMMARY", ratesRows(),
        cols("sample", "sites_covered", "bulk_rate", "expected_per_cycle", "ga_rate", "ga_share",
             "transition_rate", "transversion_rate", ...MUT_TYPES, "excluded"));
  block("PER-CODON COMPILED ATLAS", siteRows(),
        cols("protein", "hxb2_label", "nl43_position", "nt_start", "n_covered",
             "majority_aa", "majority_aa_freq", "n_aa", "entropy",
             "majority_codon", "majority_codon_freq", "n_codons",
             "mutability_median", "mutability_min", "mutability_max",
             "codon_fold_range", "fraction_on_fastest_codon", "divergence",
             "syn_fraction", "stop_fraction", "n_destinations", "effective_destinations",
             "evenness", "top_destination", "top_share"));
  // Where the supply goes, one row per (position, destination residue), weighted by the
  // batch's own codon frequencies. This is the block a reader needs to reproduce a
  // population-weighted rate for a specific amino-acid change.
  block("PER-DESTINATION WEIGHTED SUPPLY", destinationRows(),
        cols("protein", "hxb2_label", "nl43_position", "ref_aa", "alt_aa", "n_covered",
             "weighted_rate", "rate_min", "rate_max", "fold_range_across_codons",
             "fraction_reaching", "is_stop", "is_synonymous"));
  block("PER-SAMPLE COVERAGE", sampleRows(),
        cols("sample", "sites_covered", "fraction_of_region", "identity", "orientation",
             "ratio_to_batch", "excluded"));
  downloadText(`mra_${page.version}_batch_analysis.csv`, parts.join("\n"));
});

el("share-slot").appendChild(shareButton({ version: page.version,
                                           extra: { nt: () => String(markerNt) } }));

// A position named in the URL wins over the default once a batch exists.
const q = new URLSearchParams(window.location.search);
if (q.get("nt")) {
  const v = Number(q.get("nt"));
  if (Number.isFinite(v)) { markerNt = v; markerSetByUser = true; }
}
try {
  const saved = sessionStorage.getItem("mra.batch");
  if (saved) fasta.value = saved;
} catch { /* sessionStorage unavailable: nothing to restore */ }
