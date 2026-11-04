// SPDX-License-Identifier: MIT
//
// tables.js -- the table-analysis tab.
//
// Takes a pasted list of amino-acid mutations and returns the atlas rate for each. This was a
// section of "Browse the atlas" until it earned its own tab: browsing is a one-position-at-a-
// time activity, while this is a batch lookup against a list a user already has, and the two
// were competing for the same page.
//
// The parser itself lives in js/mutations.js so it can be tested without a DOM; see
// test/mutation_cases.mjs. This module is the page around it.

import { bootstrapPage, note, showError, exportProvenance } from "./page.js";
import { loadSmallAssets, loadAminoAcidLayer } from "./data.js";
import { renderTable, formatRate } from "./table.js";
import { MraEngine } from "./mra_engine.js";
import { resolveMutation } from "./mutations.js";

const page = await bootstrapPage("tables.html");
const region = page.release.modelled_region;
document.getElementById("scope-slot").appendChild(note(
  `Positions ${region.first_position.toLocaleString()}\u2013` +
  `${region.last_position.toLocaleString()} of the NL4-3 reference. Rates are mutational ` +
  `<em>supply</em> \u2014 how often a change is generated \u2014 not how often it is ` +
  `observed in patients.`));

const status = document.getElementById("tables-status");
let assets, aaLayer, engine;
try {
  assets = await loadSmallAssets(page.version);
  engine = new MraEngine({ bundle: assets.bundle, reference: assets.reference,
                           intervals: assets.intervals, hxb2: assets.hxb2 });
  aaLayer = await loadAminoAcidLayer(page.version);
  status.textContent = "";
} catch (err) {
  showError(status, "Could not load the atlas tables: " + err.message);
  throw err;
}

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
