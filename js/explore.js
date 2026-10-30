// SPDX-License-Identifier: MIT
//
// explore.js -- the interactive genome view.
//
// One smoothed supply track across the whole modelled region, and a pie of where a single
// substitution's supply goes at whichever codon you select. The track and the pie read the
// same published atlas tables; nothing here is recomputed from the model.
//
// COORDINATES. The track's x-axis is the NL4-3 reference nucleotide, because that is the
// coordinate the atlas is indexed on and it is continuous across the genome. Results are
// reported in HXB2 amino-acid numbering, which is what a clinical or publication reader
// expects. The coordinate box accepts either; the reading-frame selector resolves the
// overlapping genes (gag/pol, tat/rev/vpu/env) that a single nucleotide can belong to.

import { bootstrapPage, note, showError, shareButton } from "./page.js";
import { loadSmallAssets, loadNucleotideLayer, loadAminoAcidLayer } from "./data.js";
import { renderTable, formatRate } from "./table.js";
import { genomeTrack, pieChart } from "./chart.js";
import { MraEngine } from "./mra_engine.js";
import { perSiteSupply, smoothSeries, buildFrameIndex, geneSegments, oneStepSlices } from "./explore_core.js";

const page = await bootstrapPage("explore.html");
document.getElementById("scope-slot").appendChild(note(
  "The track shows mutational <em>supply</em> &mdash; how often a change is generated &mdash; " +
  "not how often it is observed in patients. Supply and selection are different quantities."));

const assets = await loadSmallAssets(page.version);
const ntLayer = await loadNucleotideLayer(page.version);
const aaLayer = await loadAminoAcidLayer(page.version);
const engine = new MraEngine({
  bundle: assets.bundle, reference: assets.reference,
  intervals: assets.intervals, hxb2: assets.hxb2,
});

const region = assets.reference.modelled_region;
const FIRST = region.first_position, LAST = region.last_position;

// ---------------------------------------------------------------------------------
// Per-site supply, summed over substitution classes. Indels are excluded: they are a
// different process with a different denominator, and mixing them into a "substitution
// supply" track would misreport both.
const perSite = perSiteSupply(ntLayer.byPosition);
const SITES = [...perSite.keys()].sort((a, b) => a - b);

const smoothed = (field, win) => smoothSeries(SITES, perSite, field, win);

// ---------------------------------------------------------------------------------
// Reverse index: nucleotide -> every reading frame that codes at that nucleotide.
const frameAt = buildFrameIndex(assets.reference.proteins, aaLayer.byPosition);

const GENE_SPANS = geneSegments(assets.reference.proteins);

// ---------------------------------------------------------------------------------
const elWin   = document.getElementById("explore-window");
const elWinOut= document.getElementById("explore-window-out");
const elClass = document.getElementById("explore-class");
const elSlider= document.getElementById("explore-slider");
const elPos   = document.getElementById("explore-position");
const elSys   = document.getElementById("explore-coord-system");
const elFrame = document.getElementById("explore-protein");

elSlider.min = String(FIRST); elSlider.max = String(LAST);

const CLASS_LABEL = { total: "total substitution supply", nonga: "non-G\u2192A supply", ga: "G\u2192A supply" };
const CLASS_COLOUR = { total: "#2f4b7c", nonga: "#3f8f7a", ga: "#c2561a" };

let selectedNt = 6804;          // Env 460 in HXB2 numbering -- a documented bNAb contact
let selectedFrame = null;

function drawTrack() {
  const win = Number(elWin.value);
  const cls = elClass.value;
  elWinOut.textContent = String(win);
  const pts = smoothed(cls, win);
  const hx = selectedFrame ? engine.hxb2Label(selectedFrame.protein, selectedFrame.geneAa) : null;
  genomeTrack(document.getElementById("explore-track"), [
    { label: CLASS_LABEL[cls], colour: CLASS_COLOUR[cls], points: pts, width: 1.1 },
  ], {
    title: `Smoothed ${CLASS_LABEL[cls]} along the genome`,
    xMin: FIRST, xMax: LAST, geneBand: GENE_SPANS,
    width: 1020, height: 300,
    yLabel: "supply per site per replication",
    marker: selectedNt, markerLabel: hx || String(selectedNt),
    onSelect: (pos) => selectPosition(pos, "track"),
  });

  const vals = pts.map((p) => p[1]);
  const med = vals.slice().sort((a, b) => a - b)[Math.floor(vals.length / 2)];
  const hereIdx = pts.findIndex((p) => p[0] >= selectedNt);
  const here = hereIdx >= 0 ? pts[hereIdx][1] : null;
  document.getElementById("explore-track-caption").innerHTML =
    `Centred running mean of the per-site rate over a ${win} nt window, summed across ` +
    `substitution classes${cls === "total" ? "" : " and restricted to " + CLASS_LABEL[cls]}. ` +
    `Insertions and deletions are excluded. The band above the plot shows coding segments; ` +
    `overlapping reading frames sit on separate rows, and spliced genes (tat, rev) appear as ` +
    `one box per exon. The gene holding the marker is filled solid. Genome median ${formatRate(med)}` +
    (here ? `, value at the marker ${formatRate(here)}` : "") + `. ` +
    `<strong>G\u2192A carries most of the total supply</strong>, so the total track partly ` +
    `traces G density &mdash; use the class selector to separate the two.`;
}

/** Resolve the reading frames at a nucleotide and populate the frame selector. */
function refreshFrames(nt) {
  const frames = frameAt.get(nt) || [];
  elFrame.innerHTML = "";
  for (const f of frames) {
    const o = document.createElement("option");
    o.value = f.protein + "|" + f.geneAa;
    o.textContent = engine.hxb2Label(f.protein, f.geneAa) + "  (" + f.protein + " " + f.geneAa + ")";
    elFrame.appendChild(o);
  }
  elFrame.disabled = frames.length < 2;
  return frames;
}

function selectPosition(nt, source) {
  nt = Math.min(Math.max(Math.round(nt), FIRST), LAST);
  selectedNt = nt;
  const frames = refreshFrames(nt);
  selectedFrame = frames.length ? frames[0] : null;
  if (source !== "slider") elSlider.value = String(nt);
  if (source !== "box" && elSys.value === "nt") elPos.value = String(nt);
  if (elSys.value === "hxb2" && selectedFrame) {
    const h = engine.hxb2Position(selectedFrame.protein, selectedFrame.geneAa);
    if (h !== null) elPos.value = String(h);
  }
  render();
}

function render() {
  drawTrack();
  showCodon();
}

// ---------------------------------------------------------------------------------
function showCodon() {
  const head = document.getElementById("explore-head");
  const pie = document.getElementById("explore-pie");
  const cap = document.getElementById("explore-pie-caption");
  const tableSlot = document.getElementById("explore-table");
  head.innerHTML = ""; pie.innerHTML = ""; cap.textContent = ""; tableSlot.innerHTML = "";

  if (!selectedFrame) {
    showError(head, `NL4-3 nucleotide ${selectedNt.toLocaleString()} is not inside a modelled ` +
                    `codon. The atlas covers ${FIRST.toLocaleString()}\u2013${LAST.toLocaleString()}, ` +
                    `and nef is modelled only to position 8,213.`);
    return;
  }
  const { protein, geneAa, codonIndex } = selectedFrame;
  const rows = aaLayer.byPosition.get(protein + "|" + geneAa);
  if (!rows) { showError(head, `No atlas rows at ${protein} ${geneAa}.`); return; }

  const protRec = assets.reference.proteins[protein];
  const codon = protRec.codon_seq[codonIndex];
  const codonNt = protRec.codon_nt[codonIndex];
  const hxb2 = engine.hxb2Label(protein, geneAa);
  const site = perSite.get(selectedNt);
  const which = codonNt.indexOf(selectedNt) + 1;

  head.innerHTML =
    `<p><strong>${hxb2}</strong> &mdash; reference amino acid <strong>${rows[0].ref_aa}</strong>, ` +
    `codon <code>${codon}</code>. ${protein} position ${geneAa}; NL4-3 nucleotides ` +
    `${codonNt[0].toLocaleString()}\u2013${codonNt[2].toLocaleString()}. ` +
    `Selected nucleotide ${selectedNt.toLocaleString()} is codon base ${which}` +
    (site ? `, per-site supply ${formatRate(site.total)} ` +
            `(G\u2192A ${(100 * site.ga / site.total).toFixed(0)}% of it)` : "") + `.</p>`;

  const oneStep = oneStepSlices(rows);
  if (oneStep.length) {
    pieChart(pie, oneStep.map((r) => ({
      label: `${r.alt_aa}${r.synonymous ? "*" : ""}`, value: r.mu_aa, highlight: false,
    })), { title: `Single-nucleotide amino-acid changes at ${hxb2}`, width: 680 });
    cap.innerHTML =
      "Each slice is an amino-acid change reachable from this codon by a <strong>single</strong> " +
      "nucleotide substitution, sized by its share of the one-step supply at this codon. " +
      "<code>*</code> marks a synonymous change. The slice shares are a property of the whole " +
      "codon, not of the single nucleotide selected on the track. Two- and three-step changes " +
      "are in the table below. Values are read from the atlas, not recomputed.";
  } else {
    cap.innerHTML = "No amino-acid change is reachable from this codon by a single substitution; " +
                    "the table below lists the 2- and 3-step changes.";
  }

  renderTable(tableSlot, {
    rows: rows.slice().sort((a, b) => b.mu_aa - a.mu_aa).map((r) => ({
      hxb2_change: engine.hxb2Mutation(protein, geneAa).replace(/\s\d+$/, "") + " " +
                   r.ref_aa + (engine.hxb2Position(protein, geneAa) ?? geneAa) + r.alt_aa,
      aa_change: r.aa_change, alt_aa: r.alt_aa, synonymous: r.synonymous ? "yes" : "",
      mu_aa: r.mu_aa, interval: `${formatRate(r.mu_aa_lo95)} \u2013 ${formatRate(r.mu_aa_hi95)}`,
      min_nt_changes: r.min_nt_changes, n_paths: r.n_paths,
    })),
    columns: [
      { key: "hxb2_change", label: "Mutation (HXB2)" },
      { key: "aa_change", label: "Change (gene numbering)" },
      { key: "synonymous", label: "Synonymous" },
      { key: "mu_aa", label: "Supply", type: "rate" },
      { key: "interval", label: "95% interval" },
      { key: "min_nt_changes", label: "Min nt changes", type: "num" },
      { key: "n_paths", label: "Paths", type: "num" },
    ],
    sort: { key: "mu_aa", dir: "desc" },
  });
}

// ---------------------------------------------------------------------------------
// PNG export of the title-through-pie block.
document.getElementById("explore-download").addEventListener("click", async () => {
  const btn = document.getElementById("explore-download");
  const target = document.getElementById("explore-export");
  if (typeof html2canvas === "undefined") {
    showError(document.getElementById("explore-head"),
      "The image library did not load, so PNG export is unavailable. Reload the page to retry.");
    return;
  }
  btn.disabled = true; btn.textContent = "Generating\u2026";
  try {
    const canvas = await html2canvas(target, { backgroundColor: "#ffffff", scale: 2, logging: false });
    const blob = await new Promise((res) => canvas.toBlob(res, "image/png"));
    if (!blob) throw new Error("the browser returned an empty image");
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    const hx = selectedFrame ? engine.hxb2Label(selectedFrame.protein, selectedFrame.geneAa)
                             : ("nt" + selectedNt);
    a.href = url;
    a.download = `mra_explore_${hx.replace(/\s+/g, "_")}_${page.version}.png`;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    URL.revokeObjectURL(url);
  } catch (err) {
    showError(document.getElementById("explore-head"), "PNG export failed: " + err.message);
  } finally {
    btn.disabled = false; btn.textContent = "Download as PNG";
  }
});

// ---------------------------------------------------------------------------------
elWin.addEventListener("input", () => { elWinOut.textContent = elWin.value; });
elWin.addEventListener("change", drawTrack);
elClass.addEventListener("change", drawTrack);
elSlider.addEventListener("input", () => selectPosition(Number(elSlider.value), "slider"));
elFrame.addEventListener("change", () => {
  const [protein, geneAa] = elFrame.value.split("|");
  selectedFrame = { protein, geneAa: Number(geneAa),
                    codonIndex: Number(geneAa) - assets.reference.proteins[protein].first_aa };
  render();
});

function fromBox() {
  const raw = Number(elPos.value);
  if (!Number.isFinite(raw) || raw < 1) return;
  if (elSys.value === "nt") { selectPosition(raw, "box"); return; }
  // HXB2 amino-acid input: resolve against whichever gene the current frame sits in.
  const protein = selectedFrame ? selectedFrame.protein : "Env";
  const nl43Aa = engine.nl43PositionFromHxb2(protein, raw);
  if (nl43Aa === null) {
    showError(document.getElementById("explore-head"),
      `HXB2 ${protein} position ${raw} has no NL4-3 equivalent.`);
    return;
  }
  const rec = assets.reference.proteins[protein];
  const idx = nl43Aa - rec.first_aa;
  if (idx < 0 || idx >= rec.codon_nt.length) {
    showError(document.getElementById("explore-head"),
      `HXB2 ${protein} position ${raw} is outside the modelled region.`);
    return;
  }
  selectPosition(rec.codon_nt[idx][1], "box");       // centre base of that codon
}
elPos.addEventListener("change", fromBox);
elPos.addEventListener("keydown", (e) => { if (e.key === "Enter") fromBox(); });
elSys.addEventListener("change", () => {
  if (elSys.value === "nt") elPos.value = String(selectedNt);
  else if (selectedFrame) {
    const h = engine.hxb2Position(selectedFrame.protein, selectedFrame.geneAa);
    elPos.value = h === null ? "" : String(h);
  }
});

document.getElementById("share-slot").appendChild(
  shareButton({ version: page.version, extra: { nt: () => String(selectedNt) } }));

document.getElementById("explore-status").textContent =
  `${SITES.length.toLocaleString()} modelled nucleotide positions loaded.`;
selectPosition(selectedNt, "init");
