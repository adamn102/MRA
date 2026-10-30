// SPDX-License-Identifier: MIT
//
// sequence.js -- submit a sequence, get a custom atlas.
//
// The scientific framing is load-bearing here, not decoration. This view must never
// let a reader mistake a modelled rate for a measured one, so:
//   - the standing caveat is above the input, not buried under the results;
//   - every result row carries flags naming which components it rests on;
//   - the codon-identity component is named first, because it is the one the submitted
//     sequence genuinely determines;
//   - a sequence that aligns entirely outside the modelled region gets a refusal that
//     names the covered region, not an empty table.

import { bootstrapPage, privacyPanel, note, caveat, showError, exportProvenance,
         CUSTOM_CAVEAT, NEF_NOTE } from "./page.js";
import { loadSmallAssets } from "./data.js";
import { Analysis, progressRow, MAX_BATCH_SEQUENCES } from "./analysis.js";
import { parseFasta } from "./align.js";
import { renderTable, formatRate, toCsv, downloadText } from "./table.js";

const page = await bootstrapPage("sequence.html");
const region = page.release.modelled_region;

document.getElementById("privacy-slot").appendChild(privacyPanel());
document.getElementById("caveat-slot").appendChild(note(CUSTOM_CAVEAT));
document.getElementById("batch-note").textContent =
  `Up to ${MAX_BATCH_SEQUENCES} sequences per batch.`;

const statusSlot = document.getElementById("run-status");
let assets, analysis;
try {
  // The sequence tool runs on the small assets alone -- about 360 KB gzipped. The
  // multi-megabyte browse tables are never fetched on this page.
  assets = await loadSmallAssets(page.version);
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
  document.getElementById("summary-table").innerHTML = "";
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
  const summarySlot = document.getElementById("summary-table");
  summarySlot.innerHTML = ""; detailSlot.innerHTML = "";
  if (!text) { showError(summarySlot, "Paste a sequence first."); return; }

  const records = parseFasta(text).filter((r) => r.sequence.length >= 30);
  if (!records.length) {
    showError(summarySlot, "No sequence of at least 30 bases was found in that input.");
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
    showError(summarySlot, "Analysis failed: " + err.message);
    return;
  }
  const ms = Math.round(performance.now() - t0);
  bar.finish(`Done: ${out.summaries.length} sequence(s) in ${ms} ms.`);
  document.getElementById("results-heading").hidden = false;
  if (out.truncated) {
    summarySlot.appendChild(note(
      `Only the first ${MAX_BATCH_SEQUENCES} sequences were analysed. The cap is a named ` +
      `constant (<code>MAX_BATCH_SEQUENCES</code> in <code>js/align.js</code>) and is set from ` +
      `measured timing, not guessed.`));
  }

  renderSummaries(out, records, summarySlot, detailSlot);
});

function renderSummaries(out, records, summarySlot, detailSlot) {
  const rows = out.summaries.map((s, i) => ({
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

  renderTable(summarySlot, {
    columns: [
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
      { key: "index", label: "",
        render: (r) => {
          const btn = document.createElement("button");
          btn.className = "secondary";
          btn.textContent = "Show atlas";
          btn.disabled = r.status !== "analysed";
          btn.addEventListener("click", () => showDetail(r.index, r, out, detailSlot));
          return btn;
        } },
    ],
    rows,
    filterPlaceholder: rows.length > 10 ? "filter sequences" : false,
  });

  // Anything that could not be analysed gets an explanation, never an empty table.
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
    summarySlot.appendChild(div);
  }

  if (rows.length === 1 && rows[0].status === "analysed") showDetail(0, rows[0], out, detailSlot);
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

async function showDetail(index, row, out, slot) {
  slot.innerHTML = "";
  slot.appendChild(Object.assign(document.createElement("h3"),
    { textContent: `Custom atlas for ${row.name}` }));
  const loading = caveat("Fetching the full tables\u2026");
  slot.appendChild(loading);

  let detail;
  try { detail = await analysis.detail(index); }
  catch (err) { showError(slot, "Could not retrieve results: " + err.message); return; }
  loading.remove();

  const geom = out.geometry[index];
  const s = row.summary;
  const proteinBits = Object.entries(s.proteins || {})
    .filter(([, v]) => v.codonsAnalysed > 0)
    .map(([k, v]) => `${k} ${v.firstAnalysed}\u2013${v.lastAnalysed} (${v.codonsAnalysed} of ${v.codonsInReference} modelled codons)`)
    .join("; ");
  slot.appendChild(note(
    `Placed on NL4-3 ${geom.refStart.toLocaleString()}&ndash;${geom.refEnd.toLocaleString()} at ` +
    `${(s.identity * 100).toFixed(2)}% identity, ${(geom.coverageOfQuery * 100).toFixed(1)}% of the ` +
    `submitted bases placed. Proteins covered: ${proteinBits || "none"}. ` +
    (s.droppedOutsideRegion
      ? `<strong>${s.droppedOutsideRegion.toLocaleString()} aligned positions were dropped</strong> for ` +
        `lying outside the modelled region ${region.first_position.toLocaleString()}&ndash;` +
        `${region.last_position.toLocaleString()}. `
      : "") +
    (geom.unalignedSpans.length
      ? `${geom.unalignedSpans.length} span(s) of your sequence could not be aligned: ` +
        geom.unalignedSpans.slice(0, 6).map((x) => `${x.queryStart}\u2013${x.queryEnd}`).join(", ") +
        (geom.unalignedSpans.length > 6 ? ", \u2026" : "") + ". "
      : "") +
    ((s.proteins && s.proteins.Nef && s.proteins.Nef.codonsAnalysed > 0) ? NEF_NOTE : "")));

  if (s.contextUnavailable || s.gaContextUnavailable) {
    slot.appendChild(caveat(
      `${s.contextUnavailable} site(s) had no readable 3&prime; base and fell back to the ` +
      `type-only rate; ${s.gaContextUnavailable} G&rarr;A site(s) were omitted because their ` +
      `3&prime; 2-mer could not be read \u2014 a G&rarr;A rate is not defined by this model without it. ` +
      `This happens at the very end of a submitted sequence.`));
  }

  // --- amino-acid layer ----------------------------------------------------------
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

  // --- nucleotide layer ----------------------------------------------------------
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

  // --- download everything -------------------------------------------------------
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
