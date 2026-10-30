// SPDX-License-Identifier: MIT
//
// browse.js -- the genome track, the position dashboard, and the paste-a-list analysis.
//
// One page, one selected position. The track is the spatial view and the navigator: click it,
// drag the slider, or type a coordinate, and all three converge on the same codon, whose
// dashboard appears below. (The track and the dashboard used to be separate tabs.)
//
// Everything shown here is READ from a shipped data file. Nothing published is recomputed in
// JavaScript: mu_aa, its interval, min_nt_changes and n_paths all come straight out of
// data/<version>/aa_level_supply.csv, and the individual path rates out of the per-protein
// shard. The only arithmetic is the track's running mean and turning a rate into a slice angle.
//
// COORDINATES. The track's x-axis is the NL4-3 nucleotide, because that is what the atlas is
// indexed on and it is continuous across the genome. Everything reported back is in HXB2
// amino-acid numbering, which is what a clinical or publication reader expects, with the NL4-3
// figures kept alongside for traceability.

import { bootstrapPage, note, caveat, showError, exportProvenance, shareButton, NEF_NOTE } from "./page.js";
import { loadSmallAssets, loadAminoAcidLayer, loadNucleotideLayer, loadPathShard } from "./data.js";
import { renderTable, formatRate } from "./table.js";
import { barChart, pieChart, genomeTrack } from "./chart.js";
import { MraEngine } from "./mra_engine.js";
import { resolveMutation } from "./mutations.js";
import { perSiteSupply, smoothSeries, buildFrameIndex, geneSegments, oneStepSlices,
         targetFor as targetForIn } from "./explore_core.js";


/**
 * Download the current dashboard as a PNG: query info + pie/bar chart.
 */
async function downloadDashboardPNG() {
  const button = document.getElementById("download-png");
  if (!button) return;
  
  button.disabled = true;
  button.textContent = "Generating...";
  
  try {
    // Check if html2canvas loaded
    if (typeof html2canvas === "undefined") {
      throw new Error("html2canvas library failed to load. Please refresh the page and try again.");
    }
    
    // Collect only the Position Dashboard section (title through pie chart)
    const dashboard = document.getElementById("dashboard-export");
    if (!dashboard) {
      throw new Error("Dashboard export section not found on page. Has the position been loaded?");
    }
    
    // Use html2canvas to render the dashboard as an image
    const canvas = await html2canvas(dashboard, {
      backgroundColor: "#ffffff",
      scale: 2,                    // 2x for crisp output
      logging: false,
      useCORS: true,
      allowTaint: true,
      proxy: null,                 // Disable proxy to avoid CORS issues
      timeout: 30000,              // 30 second timeout
    });
    
    // Convert canvas to PNG blob and download (properly awaited)
    return new Promise((resolve, reject) => {
      canvas.toBlob((blob) => {
        try {
          if (!blob) {
            throw new Error("Failed to create blob from canvas.");
          }
          
          const url = URL.createObjectURL(blob);
          const a = document.createElement("a");
          const date = new Date().toISOString().split("T")[0];  // YYYY-MM-DD
          a.href = url;
          a.download = `mra_atlas_${date}.png`;
          document.body.appendChild(a);
          a.click();
          document.body.removeChild(a);
          URL.revokeObjectURL(url);
          
          button.disabled = false;
          button.textContent = "Download as PNG";
          resolve();
        } catch (err) {
          reject(err);
        }
      }, "image/png");
    });
  } catch (err) {
    console.error("PNG download failed:", err);
    button.disabled = false;
    button.textContent = "Download failed";
    alert("PNG download failed: " + err.message);
    setTimeout(() => { button.textContent = "Download as PNG"; }, 3000);
  }
}

const page = await bootstrapPage("browse.html");
const region = page.release.modelled_region;
const FIRST = region.first_position, LAST = region.last_position;
document.getElementById("scope-slot").appendChild(note(
  `Positions ${FIRST.toLocaleString()}&ndash;${LAST.toLocaleString()} of the NL4-3 reference. ` +
  NEF_NOTE + " The track shows mutational <em>supply</em> &mdash; how often a change is " +
  "generated &mdash; not how often it is observed in patients."));

const status = document.getElementById("dashboard-status");
let assets, aaLayer, engine;
try {
  assets = await loadSmallAssets(page.version);
  engine = new MraEngine({ bundle: assets.bundle, reference: assets.reference,
                           intervals: assets.intervals, hxb2: assets.hxb2 });
  aaLayer = await loadAminoAcidLayer(page.version);      // ~1.2 MB gzipped, loaded once
  status.textContent = `${aaLayer.records.length.toLocaleString()} amino-acid changes loaded.`;
} catch (err) {
  showError(status, "Could not load the atlas tables: " + err.message);
  throw err;
}

// --- the protein / subunit selector ---------------------------------------------
// Gene products come from the atlas itself; the subunits (PR, RT, RNase, IN and the Gag
// cleavage products) come from the workbook's alignment sheet, in amino-acid units.
const select = document.getElementById("protein");
const targets = engine.proteinTargets(aaLayer.proteins);
for (const t of targets) {
  const opt = document.createElement("option");
  opt.value = t.key; opt.textContent = t.label;
  select.appendChild(opt);
}
const targetByKey = new Map(targets.map((t) => [t.key, t]));

// Frames are indexed off the amino-acid layer, so the track can only ever point at codons
// the atlas can answer for.
const frameAt = buildFrameIndex(assets.reference.proteins, aaLayer.byPosition);
const GENE_SEGMENTS = geneSegments(assets.reference.proteins);

const targetFor = (protein, geneAa) => targetForIn(targets, protein, geneAa);
const ownPositionFor = (target, geneAa) => engine.ownPositionFor(target, geneAa);

// --- selection state --------------------------------------------------------------
let selected = null;          // {target, geneAa, codonIndex, codonNt}
let markerNt = 1600;          // where the track's marker sits

// --- the genome track ---------------------------------------------------------------
// The nucleotide layer is several megabytes, so it loads AFTER the dashboard has painted.
// Until it arrives the track slot shows a note and its controls stay disabled.
let perSite = null, SITES = null;
const elWin = document.getElementById("track-window");
const elWinOut = document.getElementById("track-window-out");
const elClass = document.getElementById("track-class");
const elSlider = document.getElementById("track-slider");
elSlider.min = String(FIRST); elSlider.max = String(LAST);
for (const el of [elWin, elClass, elSlider]) el.disabled = true;

const CLASS_LABEL = { total: "total substitution supply", nonga: "non-G\u2192A supply", ga: "G\u2192A supply" };
const CLASS_COLOUR = { total: "#2f4b7c", nonga: "#3f8f7a", ga: "#c2561a" };

function drawTrack() {
  if (!perSite) return;
  const win = Number(elWin.value), cls = elClass.value;
  elWinOut.textContent = String(win);
  const pts = smoothSeries(SITES, perSite, cls, win);
  genomeTrack(document.getElementById("genome-track"), [
    { label: CLASS_LABEL[cls], colour: CLASS_COLOUR[cls], points: pts, width: 1.1 },
  ], {
    title: `Smoothed ${CLASS_LABEL[cls]} along the genome`,
    xMin: FIRST, xMax: LAST, geneBand: GENE_SEGMENTS,
    width: 1020, height: 300, yLabel: "supply per site per replication",
    marker: markerNt,
    markerLabel: selected ? engine.hxb2Label(selected.target.protein, selected.geneAa) : String(markerNt),
    onSelect: (pos) => selectNucleotide(pos),
  });

  const vals = pts.map((p) => p[1]).sort((a, b) => a - b);
  const med = vals[Math.floor(vals.length / 2)];
  const at = pts.find((p) => p[0] >= markerNt);
  const site = perSite.get(markerNt);
  document.getElementById("track-caption").innerHTML =
    `Centred running mean of the per-site rate over a ${win} nt window, summed across ` +
    `substitution classes${cls === "total" ? "" : " and restricted to " + CLASS_LABEL[cls]}. ` +
    `Insertions and deletions are excluded. The band above the plot shows coding segments; ` +
    `overlapping reading frames sit on separate rows, and spliced genes (tat, rev) appear as ` +
    `one box per exon, with the gene holding the marker filled solid. Genome median ` +
    `${formatRate(med)}` + (at ? `, value at the marker ${formatRate(at[1])}` : "") +
    (site ? `, of which G\u2192A is ${(100 * site.ga / site.total).toFixed(0)}% at this exact site` : "") +
    `. <strong>G\u2192A carries most of the total supply</strong>, so the total track partly ` +
    `traces G density &mdash; use &ldquo;Supply shown&rdquo; to separate the two.`;
}

async function loadTrack() {
  const slot = document.getElementById("genome-track");
  slot.innerHTML = "<p class=\"caveat\">Loading the nucleotide layer for the track\u2026</p>";
  try {
    const ntLayer = await loadNucleotideLayer(page.version);
    perSite = perSiteSupply(ntLayer.byPosition);
    SITES = [...perSite.keys()].sort((a, b) => a - b);
    for (const el of [elWin, elClass, elSlider]) el.disabled = false;
    drawTrack();
  } catch (err) {
    showError(slot, "Could not load the nucleotide layer, so the track is unavailable: " +
                    err.message + " The dashboard below still works.");
  }
}

// --- choosing a position -------------------------------------------------------------
/** From a nucleotide: resolve the reading frames there and show the first (or the one we are
 *  already in, if it still codes there -- so dragging along gag/pol does not flip frames). */
function selectNucleotide(nt) {
  nt = Math.min(Math.max(Math.round(nt), FIRST), LAST);
  const frames = frameAt.get(nt) || [];
  if (!frames.length) {
    markerNt = nt;
    elSlider.value = String(nt);
    if (document.getElementById("coord-system").value === "nt") {
      document.getElementById("position").value = String(nt);
    }
    drawTrack();
    showError(document.getElementById("dashboard-head"),
      `NL4-3 nucleotide ${nt.toLocaleString()} is not inside a modelled codon. ` +
      `The atlas covers ${FIRST.toLocaleString()}\u2013${LAST.toLocaleString()}, and nef is ` +
      `modelled only to 8,213.`);
    return;
  }
  const keep = selected && frames.find((f) => f.protein === selected.target.protein);
  const frame = keep || frames[0];
  applySelection(frame.protein, frame.geneAa, nt);
}

/** The single place that moves the selection. Everything else routes through here. */
function applySelection(protein, geneAa, preferNt) {
  const target = targetFor(protein, geneAa);
  if (!target) return;
  const rec = assets.reference.proteins[protein];
  const codonIndex = geneAa - rec.first_aa;
  const codonNt = rec.codon_nt[codonIndex];
  selected = { target, geneAa, codonIndex, codonNt };
  markerNt = (preferNt && codonNt.includes(preferNt)) ? preferNt : codonNt[1];

  elSlider.value = String(markerNt);
  refreshFrameSelector();
  syncCoordinateBox();
  drawTrack();
  showPosition(target, geneAa);
  const url = new URL(window.location.href);
  url.searchParams.set("protein", target.key);
  url.searchParams.set("pos", String(ownPositionFor(target, geneAa)));
  url.searchParams.set("nt", String(markerNt));
  url.searchParams.set("v", page.version.replace(/^v/, ""));
  window.history.replaceState({}, "", url);
  refreshShare();
}

function refreshFrameSelector() {
  const el = document.getElementById("reading-frame");
  const frames = frameAt.get(markerNt) || [];
  el.innerHTML = "";
  for (const f of frames) {
    const o = document.createElement("option");
    o.value = f.protein + "|" + f.geneAa;
    o.textContent = `${engine.hxb2Label(f.protein, f.geneAa)}  (${f.protein} ${f.geneAa})`;
    el.appendChild(o);
  }
  if (selected) el.value = selected.target.protein + "|" + selected.geneAa;
  el.disabled = frames.length < 2;
}

function syncCoordinateBox() {
  if (!selected) return;
  const protMode = document.getElementById("coord-system").value === "protein";
  select.disabled = !protMode;
  if (protMode) {
    select.value = selected.target.key;
    document.getElementById("position").value =
      String(ownPositionFor(selected.target, selected.geneAa));
  } else {
    document.getElementById("position").value = String(markerNt);
  }
}

// --- the dashboard ---------------------------------------------------------------
const shardCache = new Map();
async function shardFor(protein) {
  if (!shardCache.has(protein)) shardCache.set(protein, loadPathShard(page.version, protein));
  return shardCache.get(protein);
}

/** Render the dashboard for a gene position. `target` only affects how it is labelled. */
async function showPosition(target, genePosition) {
  const key = target.key;
  const ownPosition = ownPositionFor(target, genePosition);
  const head = document.getElementById("dashboard-head");
  const chart = document.getElementById("dashboard-chart");
  const caption = document.getElementById("dashboard-caption");
  const tableSlot = document.getElementById("dashboard-table");
  const pathSlot = document.getElementById("dashboard-paths");
  head.innerHTML = ""; chart.innerHTML = ""; caption.textContent = "";
  tableSlot.innerHTML = ""; pathSlot.innerHTML = "";
  document.getElementById("paths-heading").hidden = true;

  const rows = aaLayer.byPosition.get(target.protein + "|" + genePosition);
  if (!rows) {
    const span = aaLayer.proteins.get(target.protein);
    showError(head,
      `${key} position ${ownPosition} is outside the modelled region. ` +
      (target.protein === "Nef"
        ? "The model covers 72 of roughly 206 nef codons, because it ends at reference position 8,213."
        : `${target.protein} is modelled over positions ${span.min}\u2013${span.max} ` +
          `(${key} 1\u2013${target.maxOwn || span.max}).`));
    return;
  }

  const protRec = assets.reference.proteins[target.protein];
  const codon = protRec.codon_seq[genePosition - protRec.first_aa];
  const refAa = rows[0].ref_aa;
  const hxb2 = engine.hxb2Label(target.protein, genePosition);
  const refNt = protRec.codon_nt[genePosition - protRec.first_aa];
  const site = perSite ? perSite.get(markerNt) : null;

  head.innerHTML =
    `<p><strong>${hxb2}</strong> &mdash; reference amino acid <strong>${refAa}</strong>, ` +
    `codon <code>${codon}</code>. ` +
    `${target.protein} position ${genePosition}` +
    (target.offset ? ` (${key} ${ownPosition})` : "") +
    `; NL4-3 reference nucleotides ${refNt[0].toLocaleString()}\u2013${refNt[2].toLocaleString()}` +
    (site ? `. Selected nucleotide ${markerNt.toLocaleString()} is codon base ` +
            `${refNt.indexOf(markerNt) + 1}, per-site supply ${formatRate(site.total)}` : "") +
    `.</p>`;

  const sorted = rows.slice().sort((a, b) => b.mu_aa - a.mu_aa);
  // Pie of the ONE-nucleotide-accessible amino-acid changes: each slice is a reachable
  // target amino acid, sized by its share of the single-step mutational supply at this
  // position. One-step changes are the ones actually generated at appreciable rate, so
  // the pie shows where supply goes when a single substitution occurs here.
  const oneStep = oneStepSlices(rows);
  if (oneStep.length) {
    pieChart(chart, oneStep.map((r) => ({
      label: `${r.alt_aa}${r.synonymous ? "*" : ""}`,
      value: r.mu_aa, highlight: false,
    })), { title: `Single-nucleotide amino-acid changes at ${hxb2}`, width: 680 });
    caption.innerHTML =
      "Each slice is an amino-acid change reachable from the reference codon by a " +
      "<strong>single</strong> nucleotide substitution, sized by its share of the one-step " +
      "mutational supply at this position. <code>*</code> marks a synonymous change. The " +
      "shares are a property of the whole codon, not of the single nucleotide marked on the " +
      "track. Two- and three-step changes are in the table below but not in the pie, because " +
      "they are generated far more rarely. Values are read from the atlas, not recomputed.";
  } else {
    chart.innerHTML = "";
    caption.innerHTML =
      "No amino-acid change is reachable from this codon by a single nucleotide substitution; " +
      "see the table below for the 2- and 3-step changes.";
  }

  renderTable(tableSlot, {
    columns: [
      { key: "hxb2_change", label: "Mutation (HXB2)", title: "Clinical / publication numbering" },
      { key: "alt_aa", label: "To" },
      { key: "mu_aa", label: "mu_aa", type: "rate", title: "Summed rate over all paths to this amino acid" },
      { key: "interval", label: "95% interval" },
      { key: "min_nt_changes", label: "min nt", type: "num" },
      { key: "n_paths", label: "paths", type: "num" },
      { key: "synonymous", label: "syn" },
      { key: "aa_change", label: "NL4-3 aa change", title: "Same change in NL4-3 gene numbering" },
      { key: "ref_position", label: "NL4-3 nt", type: "num", title: "Internal coordinate, for traceability" },
    ],
    rows: sorted.map((r) => ({
      ...r, hxb2,
      interval: `${formatRate(r.mu_aa_lo95)} \u2013 ${formatRate(r.mu_aa_hi95)}`,
      ref_position: refNt[0],
    })),
    sort: { key: "mu_aa", dir: "desc" },
    downloadName: `mra_${page.version}_${key}_${ownPosition}.csv`,
    provenance: exportProvenance(page.release, [
      `Position dashboard: ${hxb2} (${target.protein} ${genePosition}), reference codon ${codon}.`,
      "Values read from aa_level_supply.csv; not recomputed.",
    ]),
  });

  // Individual paths, from this protein's shard -- loaded only now.
  document.getElementById("paths-heading").hidden = false;
  pathSlot.textContent = "Loading codon paths\u2026";
  try {
    const shard = await shardFor(target.protein);
    const paths = (shard.byPosition.get(genePosition) || []).slice()
      .sort((a, b) => b.mu_path - a.mu_path);
    renderTable(pathSlot, {
      columns: [
        { key: "mutant_codon", label: "Mutant codon" },
        { key: "alt_aa", label: "To" },
        { key: "n_nt_changes", label: "nt changes", type: "num" },
        { key: "mutation_type_path", label: "Substitutions" },
        { key: "mu_path", label: "mu_path", type: "rate" },
        { key: "interval", label: "95% interval" },
      ],
      rows: paths.map((p) => ({ ...p,
        interval: `${formatRate(p.mu_path_lo95)} \u2013 ${formatRate(p.mu_path_hi95)}` })),
      sort: { key: "mu_path", dir: "desc" },
      downloadName: `mra_${page.version}_${key}_${ownPosition}_paths.csv`,
      provenance: exportProvenance(page.release, [
        `Codon paths at ${hxb2} from codon ${codon}. All 63 alternative codons are listed.`,
      ]),
    });
    pathSlot.appendChild(caveat(
      "A path rate is the product of the per-site rates at the changed positions; its interval " +
      "is the same product over the marginal bounds. <code>mu_aa</code> above is the sum over " +
      "every path to that amino acid, including 2- and 3-substitution routes."));
  } catch (err) {
    showError(pathSlot, "Could not load codon paths: " + err.message);
  }
}

// --- coordinate entry ---------------------------------------------------------------
/** Resolve whatever is in the boxes and move the selection there. */
function jumpFromBoxes() {
  const head = document.getElementById("dashboard-head");
  const raw = Number(document.getElementById("position").value);
  if (!Number.isInteger(raw) || raw < 1) {
    showError(head, "Enter a positive whole number.");
    return;
  }
  if (document.getElementById("coord-system").value === "nt") {
    if (raw < FIRST || raw > LAST) {
      showError(head, `Nucleotide ${raw.toLocaleString()} is outside the modelled region ` +
                      `${FIRST.toLocaleString()}\u2013${LAST.toLocaleString()}.`);
      return;
    }
    selectNucleotide(raw);
    return;
  }
  const target = targetByKey.get(select.value);
  if (!target) { showError(head, "Pick a protein or subunit."); return; }
  if (target.maxOwn && raw > target.maxOwn) {
    showError(head, `${target.key} is ${target.maxOwn} residues long, so position ${raw} ` +
                    `does not exist.`);
    return;
  }
  const geneAa = engine.geneAaFromOwnPosition(target, raw);
  if (geneAa === null) {
    showError(head, `${target.key} ${raw} is an HXB2 position with no NL4-3 equivalent ` +
                    `(an insertion relative to HXB2), so the atlas has no row for it.`);
    return;
  }
  const rec = assets.reference.proteins[target.protein];
  const idx = geneAa - rec.first_aa;
  if (idx < 0 || idx >= rec.codon_nt.length ||
      !aaLayer.byPosition.has(target.protein + "|" + geneAa)) {
    showPosition(target, geneAa);         // renders the specific out-of-region message
    return;
  }
  applySelection(target.protein, geneAa, null);
}

document.getElementById("go").addEventListener("click", jumpFromBoxes);
document.getElementById("position").addEventListener("keydown",
  (e) => { if (e.key === "Enter") jumpFromBoxes(); });
select.addEventListener("change", jumpFromBoxes);
document.getElementById("coord-system").addEventListener("change", syncCoordinateBox);
document.getElementById("reading-frame").addEventListener("change", (e) => {
  const [protein, geneAa] = e.target.value.split("|");
  applySelection(protein, Number(geneAa), markerNt);
});
elSlider.addEventListener("input", () => selectNucleotide(Number(elSlider.value)));
elWin.addEventListener("input", () => { elWinOut.textContent = elWin.value; });
elWin.addEventListener("change", drawTrack);
elClass.addEventListener("change", drawTrack);
document.getElementById("download-png")?.addEventListener("click", downloadDashboardPNG);

// A shareable link reproduces the exact view, including the atlas version.
function refreshShare() {
  const slot = document.getElementById("share-slot");
  slot.innerHTML = "";
  slot.appendChild(shareButton({ version: page.version,
    extra: selected
      ? { protein: selected.target.key,
          pos: String(ownPositionFor(selected.target, selected.geneAa)),
          nt: String(markerNt) }
      : {} }));
}

// --- table analysis ---------------------------------------------------------------
// The parser lives in js/mutations.js so it can be tested without a DOM; see
// test/mutation_cases.mjs.
const mutationContext = {
  engine,
  proteins: [...aaLayer.proteins.keys()],
  lookup: (protein, genePosition) => aaLayer.byPosition.get(protein + "|" + genePosition),
};

function analyse() {
  const slot = document.getElementById("analysis-result");
  const lines = document.getElementById("mutation-list").value
    .split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
  if (!lines.length) { showError(slot, "Paste at least one mutation."); return; }
  const results = lines.map((l) => resolveMutation(l, mutationContext));
  const good = results.filter((r) => !r.error);
  const bad = results.filter((r) => r.error);

  slot.innerHTML = "";
  renderTable(slot, {
    columns: [
      { key: "line", label: "Input" },
      { key: "hxb2_change", label: "Mutation (HXB2)" },
      { key: "min_nt_changes", label: "min nt", type: "num",
        title: "Minimum number of nucleotide substitutions required" },
      { key: "mu_aa", label: "mu_aa", type: "rate" },
      { key: "interval", label: "95% interval" },
      { key: "n_paths", label: "paths", type: "num" },
      { key: "note", label: "Note" },
    ],
    rows: good.map((r) => ({ ...r,
      interval: `${formatRate(r.mu_aa_lo95)} \u2013 ${formatRate(r.mu_aa_hi95)}` })),
    sort: { key: "mu_aa", dir: "desc" },
    filterPlaceholder: false,
    downloadName: `mra_${page.version}_mutation_list.csv`,
    provenance: exportProvenance(page.release, ["Table analysis of a user-supplied mutation list."]),
  });
  if (bad.length) {
    const div = document.createElement("div");
    div.className = "error-box";
    div.innerHTML = `<strong>${bad.length} line(s) could not be looked up.</strong><ul>` +
      bad.map((r) => `<li><code>${r.line.replace(/</g, "&lt;")}</code> &mdash; ${r.error}</li>`).join("") +
      "</ul>";
    slot.appendChild(div);
  }
}

document.getElementById("analyse").addEventListener("click", analyse);
document.getElementById("load-file").addEventListener("click",
  () => document.getElementById("file-input").click());
document.getElementById("file-input").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (!file) return;
  // Read locally with FileReader; the file is never uploaded.
  document.getElementById("mutation-list").value = await file.text();
  analyse();
});

// --- open on whatever the URL asked for -------------------------------------------
// `nt` is the unambiguous coordinate and wins when present. `protein`+`pos` is kept for
// links made before the track existed; `pos` is read in the subunit's own numbering.
const qs = new URLSearchParams(window.location.search);
const wantNt = Number(qs.get("nt"));
const wantProtein = qs.get("protein");
const wantPos = Number(qs.get("pos"));
if (Number.isInteger(wantNt) && wantNt >= FIRST && wantNt <= LAST && frameAt.has(wantNt)) {
  selectNucleotide(wantNt);
} else if (wantProtein && targetByKey.has(wantProtein) && wantPos > 0) {
  select.value = wantProtein;
  document.getElementById("position").value = String(wantPos);
  jumpFromBoxes();
} else {
  select.value = "PR";
  document.getElementById("position").value = "46";
  jumpFromBoxes();
}

// The track is heavy, so it comes last -- the dashboard is already usable by now.
loadTrack();
