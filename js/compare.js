// SPDX-License-Identifier: MIT
//
// compare.js -- the view that makes the codon-identity argument concrete.
//
// It leads with CODON DIFFERENCES, not with rate differences, because the codon is the
// thing the submitted sequence actually determines. Two viruses carrying the same amino
// acid at the same position but using different synonymous codons have different sets of
// reachable amino acids and different mutabilities for each. That is the scientific point
// of the custom-atlas feature, so the first table on this page shows, for every position
// where the codon differs, which amino acids became reachable in one step and which
// stopped being reachable -- a statement that does not depend on any rate at all.
//
// The rate comparison comes second, and the genome track third.

import { bootstrapPage, privacyPanel, note, caveat, showError, exportProvenance,
         CUSTOM_CAVEAT } from "./page.js";
import { loadSmallAssets, loadAminoAcidLayer, loadNucleotideLayer } from "./data.js";
import { Analysis, progressRow } from "./analysis.js";
import { parseFasta } from "./align.js";
import { GENETIC_CODE } from "./mra_engine.js";
import { renderTable, formatRate } from "./table.js";
import { dotPlot, genomeTrack } from "./chart.js";

const page = await bootstrapPage("compare.html");
const region = page.release.modelled_region;
document.getElementById("link-sequence").href = "sequence.html?v=" + page.version.replace(/^v/, "");
document.getElementById("privacy-slot").appendChild(privacyPanel());
document.getElementById("caveat-slot").appendChild(note(CUSTOM_CAVEAT));

const statusSlot = document.getElementById("run-status");
let assets, analysis, aaLayer;
try {
  assets = await loadSmallAssets(page.version);
  analysis = new Analysis(assets);
  await analysis.start();
  aaLayer = await loadAminoAcidLayer(page.version);
  statusSlot.appendChild(caveat("Ready. The published NL4-3 amino-acid layer is loaded for comparison."));
} catch (err) {
  showError(statusSlot, "Could not load the atlas: " + err.message);
  throw err;
}

const fasta = document.getElementById("fasta");
const carried = sessionStorage.getItem("mra.sequences");
if (carried) {
  fasta.value = carried;
  document.getElementById("which").textContent = "Carried over from the sequence page in this tab.";
}
document.getElementById("load-file").addEventListener("click",
  () => document.getElementById("file-input").click());
document.getElementById("file-input").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (file) fasta.value = await file.text();
});
document.getElementById("clear").addEventListener("click", () => {
  fasta.value = "";
  sessionStorage.removeItem("mra.sequences");
  for (const id of ["codon-table", "access-chart", "access-table", "track-chart"]) {
    document.getElementById(id).innerHTML = "";
  }
  for (const id of ["codon-heading", "codon-intro", "access-heading", "track-heading"]) {
    document.getElementById(id).hidden = true;
  }
});

// --- one-step reachability of a codon, which needs no rates at all ----------------
const BASES = ["A", "C", "G", "T"];
function oneStepAminoAcids(codon) {
  const out = new Set();
  for (let slot = 0; slot < 3; slot++) {
    for (const b of BASES) {
      if (b === codon[slot]) continue;
      const alt = codon.slice(0, slot) + b + codon.slice(slot + 1);
      out.add(GENETIC_CODE[alt]);
    }
  }
  out.delete(GENETIC_CODE[codon]);     // the codon's own amino acid is not a "change"
  return out;
}

let current = null;    // {records, out, detail}

document.getElementById("run").addEventListener("click", async () => {
  const text = fasta.value.trim();
  if (!text) { showError(document.getElementById("codon-table"), "Paste a sequence first."); return; }
  const records = parseFasta(text).filter((r) => r.sequence.length >= 30);
  if (!records.length) {
    showError(document.getElementById("codon-table"), "No sequence of at least 30 bases found.");
    return;
  }
  // The comparison is a one-sequence view; a batch would make the tables meaningless.
  const first = [records[0]];
  if (records.length > 1) {
    document.getElementById("which").textContent =
      `Comparing the first of ${records.length} sequences: ${records[0].name}.`;
  }
  const bar = progressRow(document.getElementById("progress-slot"));
  let out, detail;
  try {
    out = await analysis.run(first, (p) => bar.update(p));
    if (out.summaries[0].unalignable) {
      bar.finish("Could not place the sequence.");
      showError(document.getElementById("codon-table"), out.summaries[0].reason);
      return;
    }
    detail = await analysis.detail(0);
  } catch (err) {
    showError(document.getElementById("codon-table"), "Analysis failed: " + err.message);
    return;
  }
  bar.finish(`Analysed ${records[0].name}.`);
  current = { record: records[0], out, detail };
  renderCodonDifferences();
  renderAccessibility();
  await setupTrack();
});

// ---------------------------------------------------------------------------------
function renderCodonDifferences() {
  const { detail } = current;
  const seen = new Map();
  for (const r of detail.aa) {
    const key = r.protein + "|" + r.position;
    if (seen.has(key)) continue;
    const protRec = assets.reference.proteins[r.protein];
    const refCodon = protRec.codon_seq[r.position - protRec.first_aa];
    if (refCodon === r.codon) continue;                 // identical codon: nothing to say
    const refAa = GENETIC_CODE[refCodon], yourAa = GENETIC_CODE[r.codon];
    const refSet = oneStepAminoAcids(refCodon), yourSet = oneStepAminoAcids(r.codon);
    const gained = [...yourSet].filter((a) => !refSet.has(a)).sort();
    const lost = [...refSet].filter((a) => !yourSet.has(a)).sort();
    seen.set(key, {
      hxb2_label: r.hxb2_label, protein: r.protein, position: r.position,
      ref_codon: refCodon, your_codon: r.codon,
      ref_aa: refAa, your_aa: yourAa,
      kind: refAa === yourAa ? "synonymous codon change" : "amino-acid difference",
      gained: gained.join(" ") || "\u2014",
      lost: lost.join(" ") || "\u2014",
      n_gained: gained.length, n_lost: lost.length,
      query_nt_start: r.query_nt_start,
    });
  }
  const rows = [...seen.values()];
  document.getElementById("codon-heading").hidden = false;
  const intro = document.getElementById("codon-intro");
  intro.hidden = false;
  const syn = rows.filter((r) => r.kind === "synonymous codon change").length;
  const reshaped = rows.filter((r) => r.kind === "synonymous codon change" && (r.n_gained || r.n_lost)).length;
  intro.innerHTML =
    `<strong>${rows.length.toLocaleString()} analysed codons differ from NL4-3</strong>, of which ` +
    `${syn.toLocaleString()} encode the <em>same amino acid</em> by a different codon. ` +
    `<strong>${reshaped.toLocaleString()} of those synonymous differences change which amino acids are ` +
    `reachable in a single nucleotide substitution.</strong> Those rows are the codon-identity effect ` +
    `in its purest form: the protein sequence is unchanged, and the set of accessible next steps is not. ` +
    `Gained and lost are computed from the genetic code alone and do not depend on any rate estimate.`;

  renderTable(document.getElementById("codon-table"), {
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
    downloadName: `mra_${page.version}_codon_differences.csv`,
    provenance: exportProvenance(page.release,
      [`Codon differences between submitted sequence "${current.record.name}" and NL4-3.`,
       "Gained/lost reachability is derived from the genetic code, not from rate estimates."]),
  });
}

// ---------------------------------------------------------------------------------
function renderAccessibility() {
  const { detail } = current;
  const rows = [];
  for (const r of detail.aa) {
    const published = (aaLayer.byPosition.get(r.protein + "|" + r.position) || [])
      .find((p) => p.alt_aa === r.alt_aa);
    const protRec = assets.reference.proteins[r.protein];
    const refCodon = protRec.codon_seq[r.position - protRec.first_aa];
    const codonDiffers = refCodon !== r.codon;
    if (!published) {
      // Reachable from your codon but absent from the published table for this position:
      // only possible when your codon differs, and worth showing as newly available.
      rows.push({ hxb2_label: r.hxb2_label, aa_change: r.aa_change, alt_aa: r.alt_aa,
                  your_codon: r.codon, ref_codon: refCodon, codon_differs: codonDiffers,
                  mu_custom: r.mu_aa, mu_nl43: null, log2fc: null,
                  min_nt_custom: r.min_nt_changes, min_nt_nl43: null,
                  change: "not in the NL4-3 table at this position" });
      continue;
    }
    const log2fc = Math.log2(r.mu_aa / published.mu_aa);
    rows.push({
      hxb2_label: r.hxb2_label, aa_change: r.aa_change, alt_aa: r.alt_aa,
      your_codon: r.codon, ref_codon: refCodon, codon_differs: codonDiffers,
      mu_custom: r.mu_aa, mu_nl43: published.mu_aa, log2fc,
      min_nt_custom: r.min_nt_changes, min_nt_nl43: published.min_nt_changes,
      change: r.min_nt_changes < published.min_nt_changes ? "fewer substitutions needed"
            : r.min_nt_changes > published.min_nt_changes ? "more substitutions needed"
            : "",
    });
  }

  document.getElementById("access-heading").hidden = false;
  const withFc = rows.filter((r) => r.log2fc !== null && Number.isFinite(r.log2fc));
  const top = withFc.slice().sort((a, b) => Math.abs(b.log2fc) - Math.abs(a.log2fc)).slice(0, 25);
  dotPlot(document.getElementById("access-chart"),
    top.map((r) => ({ label: `${r.hxb2_label} ${r.aa_change}`, value: r.log2fc,
                      note: `${r.ref_codon}\u2192${r.your_codon}` })),
    { title: "The 25 amino-acid changes whose accessibility shifts most", width: 740 });
  const shifted = withFc.filter((r) => Math.abs(r.log2fc) > 1).length;
  const stepChanged = rows.filter((r) => r.change && r.change.includes("substitutions")).length;
  document.getElementById("access-caption").innerHTML =
    `${withFc.length.toLocaleString()} amino-acid changes are present in both atlases. ` +
    `${shifted.toLocaleString()} shift by more than two-fold, and ${stepChanged.toLocaleString()} ` +
    `require a different <em>number</em> of nucleotide substitutions than in NL4-3 \u2014 which is a ` +
    `categorical change in accessibility, not a quantitative one. Codon shown at the right of each row.`;

  renderTable(document.getElementById("access-table"), {
    columns: [
      { key: "hxb2_label", label: "HXB2" },
      { key: "aa_change", label: "Change" },
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
    downloadName: `mra_${page.version}_accessibility_vs_nl43.csv`,
    provenance: exportProvenance(page.release,
      [`Accessibility of amino-acid changes in "${current.record.name}" relative to NL4-3.`,
       "log2 FC = log2(mu_aa custom / mu_aa published NL4-3).",
       "A difference in min_nt_changes is a categorical change in accessibility."]),
  });
  document.getElementById("access-table").appendChild(caveat(
    "A fold change here reflects your codon and your local context. It is not evidence that " +
    "your virus mutates at a different rate: the underlying per-site model is the same one, " +
    "fitted on NL4-3."));
}

// ---------------------------------------------------------------------------------
let ntLayer = null;
async function setupTrack() {
  const select = document.getElementById("track-type");
  if (!select.options.length) {
    const types = ["A->C", "A->G", "A->T", "C->A", "C->G", "C->T",
                   "G->A", "G->C", "G->T", "T->A", "T->C", "T->G"];
    for (const t of types) {
      const opt = document.createElement("option");
      opt.value = t; opt.textContent = t;
      select.appendChild(opt);
    }
    select.value = "G->A";
    select.addEventListener("change", drawTrack);
    document.getElementById("track-draw").addEventListener("click", drawTrack);
  }
  document.getElementById("track-heading").hidden = false;
  if (!ntLayer) {
    document.getElementById("track-caption").textContent = "Loading the published nucleotide layer\u2026";
    ntLayer = await loadNucleotideLayer(page.version);    // ~0.8 MB gzipped, loaded on demand
  }
  drawTrack();
}

/** Running mean in log space over a window of positions, for readability at genome scale. */
function smooth(points, window) {
  if (window <= 1 || points.length < 3) return points;
  const out = [];
  for (let i = 0; i < points.length; i++) {
    let sum = 0, n = 0;
    for (let j = i; j >= 0 && points[i][0] - points[j][0] <= window / 2; j--) { sum += Math.log(points[j][1]); n++; }
    for (let j = i + 1; j < points.length && points[j][0] - points[i][0] <= window / 2; j++) { sum += Math.log(points[j][1]); n++; }
    out.push([points[i][0], Math.exp(sum / n)]);
  }
  return out;
}

function drawTrack() {
  if (!current || !ntLayer) return;
  const type = document.getElementById("track-type").value;
  const window = Math.max(1, Number(document.getElementById("track-smooth").value) || 1);

  const nl43 = [];
  for (const r of ntLayer.records) {
    if (r.mut_type === type && r.bio_rate_median > 0) nl43.push([r.ref_position, r.bio_rate_median]);
  }
  nl43.sort((a, b) => a[0] - b[0]);

  const yours = current.detail.nt
    .filter((r) => r.mut_type === type && r.rate > 0)
    .map((r) => [r.ref_position, r.rate])
    .sort((a, b) => a[0] - b[0]);

  genomeTrack(document.getElementById("track-chart"), [
    { label: "NL4-3 (published)", colour: "#9aa5ad", points: smooth(nl43, window), width: 1 },
    { label: "your sequence", colour: "#c2561a", points: smooth(yours, window), width: 1.2 },
  ], {
    title: `${type} rate along the genome`,
    xMin: region.first_position, xMax: region.last_position,
    spans: geneSpans(),
    width: 1020, height: 280,
  });

  const onlyYours = yours.length - nl43.filter((p) => yours.some((q) => q[0] === p[0])).length;
  document.getElementById("track-caption").innerHTML =
    `Running geometric mean over a ${window} nt window. NL4-3 has ${nl43.length.toLocaleString()} ` +
    `sites where ${type} is possible; your sequence has ${yours.length.toLocaleString()} within the ` +
    `region it covers. <strong>The two series are defined at different sets of positions</strong>, ` +
    `because a substitution only exists where the sequence actually carries that base \u2014 which is ` +
    `itself part of what differs between viruses. Grey shading marks genes; the ` +
    `track is restricted to the modelled region ${region.first_position.toLocaleString()}\u2013` +
    `${region.last_position.toLocaleString()}.`;
}

/** Gene spans for shading, from the per-position annotation shipped with the reference. */
function geneSpans() {
  const levels = assets.reference.annotation_levels;
  const codes = assets.reference.annotation_by_position;
  const wanted = ["gag", "pol", "vif", "vpr", "tat", "rev", "env", "nef"];
  const spans = [];
  for (const gene of wanted) {
    const idx = levels.indexOf(gene);
    if (idx < 0) continue;
    let start = null;
    for (let p = region.first_position; p <= region.last_position; p++) {
      const here = codes[p - 1] === idx;
      if (here && start === null) start = p;
      if (!here && start !== null) { spans.push({ start, end: p - 1, label: "" }); start = null; }
    }
    if (start !== null) spans.push({ start, end: region.last_position, label: "" });
    // Label the widest span of each gene only, so the strip stays readable.
    const mine = spans.filter((s) => s.label === "");
    if (mine.length) {
      mine.sort((a, b) => (b.end - b.start) - (a.end - a.start));
      mine[0].label = gene;
      for (const s of mine.slice(1)) s.label = " ";
    }
  }
  return spans;
}
