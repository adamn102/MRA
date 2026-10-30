// SPDX-License-Identifier: MIT
//
// browse.js -- the position dashboard and the paste-a-list table analysis.
//
// Everything shown here is READ from a shipped data file. Nothing published is
// recomputed in JavaScript: mu_aa, its interval, min_nt_changes and n_paths all come
// straight out of data/<version>/aa_level_supply.csv, and the individual path rates out
// of the per-protein shard. The only arithmetic is turning a rate into a bar length.

import { bootstrapPage, note, caveat, showError, exportProvenance, shareButton, NEF_NOTE } from "./page.js";
import { loadSmallAssets, loadAminoAcidLayer, loadPathShard } from "./data.js";
import { renderTable, formatRate } from "./table.js";
import { barChart, pieChart } from "./chart.js";
import { MraEngine } from "./mra_engine.js";
import { resolveMutation } from "./mutations.js";


/**
 * Download the current dashboard as a PNG: query info + pie/bar chart.
 */
async function downloadDashboardPNG() {
  const button = document.getElementById("download-png");
  if (!button) return;
  
  button.disabled = true;
  button.textContent = "Generating...";
  
  try {
    // Collect the entire dashboard section: query info + chart
    const dashboard = document.querySelector("main") || document.body;
    
    // Use html2canvas to render the dashboard as an image
    const canvas = await html2canvas(dashboard, {
      backgroundColor: "#ffffff",
      scale: 2,                    // 2x for crisp output
      logging: false,
      useCORS: true,
      allowTaint: true,
    });
    
    // Convert canvas to PNG blob and download
    canvas.toBlob((blob) => {
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
    });
  } catch (err) {
    console.error("PNG download failed:", err);
    button.disabled = false;
    button.textContent = "Download failed";
    setTimeout(() => { button.textContent = "Download as PNG"; }, 3000);
  }
}
const page = await bootstrapPage("browse.html");
const region = page.release.modelled_region;
document.getElementById("scope-slot").appendChild(note(
  `Positions ${region.first_position.toLocaleString()}&ndash;${region.last_position.toLocaleString()} ` +
  `of the NL4-3 reference. ` + NEF_NOTE));

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
const targets = [];
for (const protein of [...aaLayer.proteins.keys()].sort()) {
  const span = aaLayer.proteins.get(protein);
  targets.push({ key: protein, label: `${protein}  (1\u2013${span.max})`, protein, offset: 0 });
  for (const sp of engine.subProteins[protein] || []) {
    targets.push({ key: sp.abbrev, protein,
                   label: `\u2003${sp.abbrev} \u2014 ${sp.name} (1\u2013${sp.lastGeneAa - sp.offset})`,
                   offset: sp.offset, maxOwn: sp.lastGeneAa - sp.offset });
  }
}
for (const t of targets) {
  const opt = document.createElement("option");
  opt.value = t.key; opt.textContent = t.label;
  select.appendChild(opt);
}
const targetByKey = new Map(targets.map((t) => [t.key, t]));

// --- the dashboard ---------------------------------------------------------------
const shardCache = new Map();
async function shardFor(protein) {
  if (!shardCache.has(protein)) shardCache.set(protein, loadPathShard(page.version, protein));
  return shardCache.get(protein);
}

async function showPosition(key, ownPosition) {
  const target = targetByKey.get(key);
  const head = document.getElementById("dashboard-head");
  const chart = document.getElementById("dashboard-chart");
  const caption = document.getElementById("dashboard-caption");
  const tableSlot = document.getElementById("dashboard-table");
  const pathSlot = document.getElementById("dashboard-paths");
  head.innerHTML = ""; chart.innerHTML = ""; caption.textContent = "";
  tableSlot.innerHTML = ""; pathSlot.innerHTML = "";
  document.getElementById("paths-heading").hidden = true;

  const genePosition = ownPosition + target.offset;
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

  head.innerHTML =
    `<p><strong>${hxb2}</strong> &mdash; reference amino acid <strong>${refAa}</strong>, ` +
    `codon <code>${codon}</code>. ` +
    `${target.protein} position ${genePosition}` +
    (target.offset ? ` (${key} ${ownPosition})` : "") +
    `; NL4-3 reference nucleotides ${refNt[0].toLocaleString()}\u2013${refNt[2].toLocaleString()}.</p>`;

  const sorted = rows.slice().sort((a, b) => b.mu_aa - a.mu_aa);
  // Pie of the ONE-nucleotide-accessible amino-acid changes: each slice is a reachable
  // target amino acid, sized by its share of the single-step mutational supply at this
  // position. One-step changes are the ones actually generated at appreciable rate, so
  // the pie shows where supply goes when a single substitution occurs here.
  const oneStep = sorted.filter((r) => r.min_nt_changes === 1);
  if (oneStep.length) {
    pieChart(chart, oneStep.map((r) => ({
      label: `${r.alt_aa}${r.synonymous ? "*" : ""}`,
      value: r.mu_aa, highlight: false,
    })), { title: `Single-nucleotide amino-acid changes at ${hxb2}`, width: 680 });
    caption.innerHTML =
      "Each slice is an amino-acid change reachable from the reference codon by a " +
      "<strong>single</strong> nucleotide substitution, sized by its share of the one-step " +
      "mutational supply at this position. <code>*</code> marks a synonymous change. " +
      "Two- and three-step changes are in the table below but not in the pie, because they " +
      "are generated far more rarely. Values are read from the atlas, not recomputed.";
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

document.getElementById("go").addEventListener("click", () => {
  const key = select.value;
  const posInput = Number(document.getElementById("position").value);
  const coordSystem = document.getElementById("coord-system")?.value || "NL4-3";
  
  if (!Number.isInteger(posInput) || posInput < 1) {
    showError(document.getElementById("dashboard-head"), "Enter a positive integer position.");
    return;
  }
  
  // Convert HXB2 to NL4-3 if needed
  let pos = posInput;
  if (coordSystem === "HXB2") {
    // Reverse lookup: find NL4-3 position that maps to this HXB2 position
    const nl43Pos = engine.nl43PositionFromHxb2(key, posInput);
    if (nl43Pos === null) {
      showError(document.getElementById("dashboard-head"), `HXB2 position ${posInput} not found in ${key}. Try NL4-3 numbering.`);
      return;
    }
    pos = nl43Pos;
  }
  const url = new URL(window.location.href);
  url.searchParams.set("protein", key);
  url.searchParams.set("pos", String(pos));
  url.searchParams.set("v", page.version.replace(/^v/, ""));
  window.history.replaceState({}, "", url);
  showPosition(key, pos);
});
document.getElementById("download-png")?.addEventListener("click", downloadDashboardPNG);
document.getElementById("position").addEventListener("keydown",
  (e) => { if (e.key === "Enter") document.getElementById("go").click(); });

// A shareable link reproduces the exact view, including the atlas version.
const refreshShare = () => {
  const slot = document.getElementById("share-slot");
  slot.innerHTML = "";
  slot.appendChild(shareButton({ version: page.version,
    extra: { protein: select.value, pos: document.getElementById("position").value } }));
};
select.addEventListener("change", refreshShare);
document.getElementById("position").addEventListener("input", refreshShare);
refreshShare();

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
const qs = new URLSearchParams(window.location.search);
const wantProtein = qs.get("protein");
const wantPos = Number(qs.get("pos"));
if (wantProtein && targetByKey.has(wantProtein) && wantPos > 0) {
  select.value = wantProtein;
  document.getElementById("position").value = String(wantPos);
  showPosition(wantProtein, wantPos);
} else {
  select.value = "PR";
  showPosition("PR", 46);
}
refreshShare();
