// SPDX-License-Identifier: MIT
//
// tables_custom.js -- the custom-query section of the MRA Table page.
//
// Takes a pasted list of amino-acid mutations and returns the atlas rate for each. It sits
// last on the page, after the curated tables: a reader who wants PI or NNRTI finds it already
// compiled above, and this is for the list they brought with them.
//
// The parser lives in js/mutations.js so it can be tested without a DOM; see
// test/mutation_cases.mjs. This module is the section around it, and it is a function rather
// than top-level code so that one page can host it alongside the curated half without the two
// competing to bootstrap the page.

import { showError, exportProvenance } from "./page.js";
import { renderTable, formatRate } from "./table.js";
import { resolveMutation } from "./mutations.js";

export function initCustomQuery({ page, engine, aaLayer }) {
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
      provenance: exportProvenance(page.release,
        ["Custom mutation list analysed against the atlas."]),
    });
    if (bad.length) {
      const div = document.createElement("div");
      div.className = "error-box";
      div.innerHTML = `<strong>${bad.length} line(s) could not be looked up.</strong><ul>` +
        bad.map((r) => `<li><code>${r.line.replace(/</g, "&lt;")}</code> &mdash; ${r.error}</li>`)
          .join("") + "</ul>";
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
}
