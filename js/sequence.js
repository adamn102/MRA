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
import { MraEngine } from "./mra_engine.js";
import { perSiteSupply, groupByPosition, buildFrameIndex, geneSegments, targetFor }
  from "./explore_core.js";
import { drawSupplyTrack, drawPositionPie } from "./track_view.js";

const page = await bootstrapPage("sequence.html");
const region = page.release.modelled_region;

document.getElementById("privacy-slot").appendChild(privacyPanel());
document.getElementById("caveat-slot").appendChild(note(CUSTOM_CAVEAT));
document.getElementById("batch-note").textContent =
  `Up to ${MAX_BATCH_SEQUENCES} sequences per batch.`;

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


/**
 * The interactive genome view for one submitted sequence: the same track and pie as the
 * published atlas, but every number comes from THIS sequence's custom atlas rather than from
 * NL4-3. Ticks along the baseline mark where the submission differs from the reference, so
 * the genome can be scanned for variants and each one opened.
 *
 * Elements are built as local references rather than by id, because a batch renders one of
 * these per sequence and ids would collide.
 */
function renderInteractiveView(slot, detail, row) {
  const heading = document.createElement("h4");
  heading.textContent = "Scan the genome";
  slot.appendChild(heading);

  // --- per-site supply, from this sequence's own rates -----------------------------
  const perSite = perSiteSupply(groupByPosition(detail.nt), "rate");
  const SITES = [...perSite.keys()].sort((a, b) => a - b);
  if (!SITES.length) {
    slot.appendChild(caveat("No modelled nucleotide position was covered, so there is " +
                            "nothing to plot for this sequence."));
    return;
  }

  // --- where this sequence differs from NL4-3 --------------------------------------
  // The custom atlas stores the SUBMITTED base at each aligned reference position, so a
  // mismatch against the reference string is exactly a substitution in the submission.
  // Indels are not included: they have no single reference position to mark.
  const refSeq = assets.reference.sequence;
  const baseAt = new Map();
  for (const r of detail.nt) if (!baseAt.has(r.ref_position)) baseAt.set(r.ref_position, r.ref_nucleotide);
  const variants = [...baseAt.entries()]
    .filter(([pos, base]) => base !== refSeq[pos - 1]).map(([pos]) => pos).sort((a, b) => a - b);

  // --- codon lookup for this sequence ----------------------------------------------
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
  const frameAt = buildFrameIndex(assets.reference.proteins, aaByPosition);
  const SEGMENTS = geneSegments(assets.reference.proteins);
  const targets = engine.proteinTargets(proteinSpans);
  const targetByKey = new Map(targets.map((t) => [t.key, t]));

  // --- controls ---------------------------------------------------------------------
  const ctrl = document.createElement("div");
  ctrl.className = "controls";
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
    `<div><button class="sq-go">Go</button></div>` +
    `<div><label>Reading frame here</label><select class="sq-frame"></select></div>`;
  slot.appendChild(ctrl);

  const elWin = ctrl.querySelector(".sq-win"), elWinOut = ctrl.querySelector(".sq-win-out");
  const elClass = ctrl.querySelector(".sq-class"), elSys = ctrl.querySelector(".sq-sys");
  const elTarget = ctrl.querySelector(".sq-target"), elPos = ctrl.querySelector(".sq-pos");
  const elFrame = ctrl.querySelector(".sq-frame");
  for (const t of targets) {
    const o = document.createElement("option");
    o.value = t.key; o.textContent = t.label;
    elTarget.appendChild(o);
  }

  // --- figure ------------------------------------------------------------------------
  const fig = document.createElement("figure");
  fig.className = "chart";
  const trackSlot = document.createElement("div");
  const sliderWrap = document.createElement("div");
  sliderWrap.className = "track-slider";
  const elSlider = document.createElement("input");
  elSlider.type = "range"; elSlider.step = "1";
  elSlider.min = String(SITES[0]); elSlider.max = String(SITES[SITES.length - 1]);
  elSlider.setAttribute("aria-label", "Genome position (NL4-3 nucleotide)");
  sliderWrap.appendChild(elSlider);
  const figCap = document.createElement("figcaption");
  fig.append(trackSlot, sliderWrap, figCap);
  slot.appendChild(fig);

  const head = document.createElement("div");
  const pieFig = document.createElement("figure");
  pieFig.className = "chart";
  const pieSlot = document.createElement("div");
  const pieCap = document.createElement("figcaption");
  pieFig.append(pieSlot, pieCap);
  slot.append(head, pieFig);

  // --- state ---------------------------------------------------------------------------
  let selected = null, markerNt = variants.length ? variants[0] : SITES[0];
  let variantIdx = 0;

  function draw() {
    const win = Number(elWin.value);
    elWinOut.textContent = String(win);
    figCap.innerHTML = drawSupplyTrack(trackSlot, {
      perSite, sites: SITES, win, cls: elClass.value, segments: SEGMENTS,
      first: region.first_position, last: region.last_position,
      marker: markerNt,
      markerLabel: selected ? engine.hxb2Label(selected.protein, selected.geneAa) : String(markerNt),
      variants,
      title: `Smoothed supply along ${row.name}`,
      source: "computed from the codons and 3\u2032 context in your sequence, not from NL4-3",
      onSelect: (pos) => selectNt(pos),
    });
  }

  function selectNt(nt) {
    const frames = frameAt.get(Math.round(nt)) || [];
    markerNt = Math.round(nt);
    elSlider.value = String(markerNt);
    if (elSys.value === "nt") elPos.value = String(markerNt);
    if (!frames.length) {
      selected = null;
      head.innerHTML = "";
      showError(head, `NL4-3 nucleotide ${markerNt.toLocaleString()} is covered by your ` +
                      `sequence but is not inside a codon this atlas models.`);
      pieSlot.innerHTML = ""; pieCap.textContent = "";
      draw();
      return;
    }
    const keep = selected && frames.find((f) => f.protein === selected.protein);
    apply((keep || frames[0]).protein, (keep || frames[0]).geneAa, markerNt);
  }

  function apply(protein, geneAa, preferNt) {
    const rec = assets.reference.proteins[protein];
    const codonNt = rec.codon_nt[geneAa - rec.first_aa];
    selected = { protein, geneAa };
    markerNt = (preferNt && codonNt.includes(preferNt)) ? preferNt : codonNt[1];
    elSlider.value = String(markerNt);

    elFrame.innerHTML = "";
    for (const f of (frameAt.get(markerNt) || [])) {
      const o = document.createElement("option");
      o.value = f.protein + "|" + f.geneAa;
      o.textContent = `${engine.hxb2Label(f.protein, f.geneAa)}  (${f.protein} ${f.geneAa})`;
      elFrame.appendChild(o);
    }
    elFrame.value = protein + "|" + geneAa;
    elFrame.disabled = (frameAt.get(markerNt) || []).length < 2;

    const t = targetFor(targets, protein, geneAa);
    if (t && elSys.value === "protein") {
      elTarget.value = t.key;
      elPos.value = String(engine.ownPositionFor(t, geneAa));
    }

    const rows = aaByPosition.get(protein + "|" + geneAa) || [];
    const label = engine.hxb2Label(protein, geneAa);
    const yourCodon = rows.length ? rows[0].codon : "?";
    const refCodon = rec.codon_seq[geneAa - rec.first_aa];
    const differs = yourCodon !== refCodon;
    head.innerHTML =
      `<p><strong>${label}</strong> &mdash; your codon <code>${yourCodon}</code>` +
      (differs ? ` <strong>(NL4-3 has <code>${refCodon}</code>)</strong>` : " (same as NL4-3)") +
      `, amino acid <strong>${rows.length ? rows[0].ref_aa : "?"}</strong>. ` +
      `${protein} position ${geneAa}; NL4-3 nucleotides ` +
      `${codonNt[0].toLocaleString()}\u2013${codonNt[2].toLocaleString()}; ` +
      `selected nucleotide ${markerNt.toLocaleString()} is codon base ` +
      `${codonNt.indexOf(markerNt) + 1}.</p>`;
    drawPositionPie(pieSlot, pieCap, rows, label, { custom: true });
    draw();
  }

  function jump() {
    head.innerHTML = "";
    if (elSys.value === "variant") {
      if (!variants.length) { showError(head, "This sequence matches NL4-3 at every modelled position."); return; }
      const want = Number(elPos.value);
      variantIdx = Number.isInteger(want) && want >= 1 && want <= variants.length
        ? want - 1 : (variantIdx + 1) % variants.length;
      elPos.value = String(variantIdx + 1);
      selectNt(variants[variantIdx]);
      return;
    }
    const raw = Number(elPos.value);
    if (!Number.isInteger(raw) || raw < 1) { showError(head, "Enter a positive whole number."); return; }
    if (elSys.value === "nt") {
      if (!perSite.has(raw)) {
        showError(head, `NL4-3 nucleotide ${raw.toLocaleString()} is not covered by this ` +
                        `sequence. It covers ${SITES[0].toLocaleString()}\u2013` +
                        `${SITES[SITES.length - 1].toLocaleString()} of the reference.`);
        return;
      }
      selectNt(raw);
      return;
    }
    const t = targetByKey.get(elTarget.value);
    if (!t) { showError(head, "Pick a protein or subunit."); return; }
    const geneAa = engine.geneAaFromOwnPosition(t, raw);
    if (geneAa === null || !aaByPosition.has(t.protein + "|" + geneAa)) {
      const span = proteinSpans.get(t.protein);
      showError(head, `${t.key} ${raw} is not covered by this sequence` +
        (span ? ` \u2014 ${t.protein} ${span.min}\u2013${span.max} was analysed.` : "."));
      return;
    }
    apply(t.protein, geneAa, null);
  }

  elWin.addEventListener("input", () => { elWinOut.textContent = elWin.value; });
  elWin.addEventListener("change", draw);
  elClass.addEventListener("change", draw);
  elSlider.addEventListener("input", () => selectNt(Number(elSlider.value)));
  ctrl.querySelector(".sq-go").addEventListener("click", jump);
  elPos.addEventListener("keydown", (e) => { if (e.key === "Enter") jump(); });
  elTarget.addEventListener("change", jump);
  elFrame.addEventListener("change", (e) => {
    const [p, g] = e.target.value.split("|");
    apply(p, Number(g), markerNt);
  });
  elSys.addEventListener("change", () => {
    const mode = elSys.value;
    elTarget.disabled = mode !== "protein";
    elPos.value = mode === "nt" ? String(markerNt)
      : mode === "variant" ? String(variantIdx + 1)
      : (selected ? String(engine.ownPositionFor(targetFor(targets, selected.protein, selected.geneAa), selected.geneAa)) : "1");
  });

  slot.appendChild(caveat(
    variants.length
      ? `<strong>${variants.length.toLocaleString()} position(s)</strong> in the modelled region ` +
        `differ from NL4-3; each is a red tick on the baseline. Set &ldquo;Jump to&rdquo; to ` +
        `&ldquo;Difference from NL4-3&rdquo; and press Go to step through them. A tick marks a ` +
        `substitution, not an indel \u2014 an indel has no single reference position to mark. ` +
        `These are differences in <em>your</em> sequence, not predictions: the track shows how ` +
        `readily each site mutates, which is a separate question from whether it already has.`
      : "This sequence matches NL4-3 at every modelled position, so there are no difference " +
        "ticks. The track still shows how readily each site mutates."));

  elTarget.disabled = false;
  selectNt(markerNt);
}

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

  // --- interactive genome view ------------------------------------------------------
  renderInteractiveView(slot, detail, row);

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
