// SPDX-License-Identifier: MIT
//
// reference.js -- the "Reference and Download" tab.
//
// The Overview and Download tabs were two halves of the same thing: one said what the numbers
// mean, the other handed over the files, and a reader needed both before trusting either. They
// are now one page, and it is where the portal's fine print lives -- the privacy statement,
// the supply-versus-frequency distinction, the coordinate caveat and the notes on reading the
// curated tables -- so no other page has to repeat it.
//
// This module is index.js and download.js merged. The download half is unchanged: a file must
// be self-identifying, so every CSV goes out with a header naming the atlas version, the model
// and the retrieval date, and is otherwise byte-identical to the shipped asset.

import { bootstrapPage, privacyPanel, note, caveat, exportProvenance, NEF_NOTE } from "./page.js";
import { assetUrl, verifyAssetChecksum, citationText } from "./version.js";
import { renderTable, downloadText } from "./table.js";

const page = await bootstrapPage("reference.html");
document.getElementById("privacy-slot").appendChild(privacyPanel());

const region = page.release.modelled_region;
document.getElementById("scope-slot").appendChild(note(
  `The portal covers <strong>NL4-3 reference positions ${region.first_position.toLocaleString()}` +
  `&ndash;${region.last_position.toLocaleString()}</strong>, where position 1 is the gag start ` +
  `codon. That is the extent of the fitted model. ` + NEF_NOTE
));

// Every link carries its version explicitly, so a copied URL keeps resolving to this release.
const v = page.version.replace(/^v/, "");
for (const [id, href] of [["link-browse", "index.html"], ["link-tables", "tables.html"],
                          ["link-sequence", "sequence.html"], ["link-batch", "batch.html"]]) {
  document.getElementById(id).href = href + "?v=" + v;
}

const release = page.release;

document.getElementById("version-info").appendChild(note(
  `Atlas <strong>${release.version}</strong>, released ${release.release_date || "\u2014"}. ` +
  `Model <code>${release.model}</code>. ` +
  (release.doi ? `DOI ${release.doi}.` : "DOI not yet assigned.") +
  ` Covers NL4-3 reference positions ${release.modelled_region.first_position.toLocaleString()}` +
  `&ndash;${release.modelled_region.last_position.toLocaleString()}.`));

const DESCRIPTIONS = {
  "nt_level_supply.csv":
    "Published nucleotide layer: one row per reference position and mutation type, including insertions and deletions.",
  "aa_level_supply.csv":
    "Published amino-acid mutability: one row per protein, position and target amino acid.",
  "mra_web_model.json":
    "The fitted model. Coefficients, context deltas, G->A intercepts and six position splines.",
  "interval_ratios.json":
    "95% interval ratios, so a rate computed for a submitted sequence can carry an interval.",
  "lanl_flow.json":
    "Precomputed batch over the 182 LANL group M compendium genomes: per amino-acid position, the codons circulating virus carries and the rate each supplies to each destination residue. Backs the supply-flow chart on the landing page.",
  "nl43_reference.json":
    "Reference sequence, per-position protein assignment and codon frames (including the tat/rev splice junctions).",
  "hxb2_align.json":
    "HXB2 <-> NL4-3 amino-acid mapping for Env and Vpu, plus protein offsets within genes.",
  "citation.json": "Machine-readable citation block.",
  "manifest.json": "SHA-256 and gzipped size of every asset in this release.",
};

async function fetchAndDownload(filename) {
  const resp = await fetch(assetUrl(page.version, filename));
  if (!resp.ok) throw new Error("HTTP " + resp.status);
  const body = await resp.text();
  const stamped = filename.endsWith(".csv")
    ? exportProvenance(release, [`Source file: ${filename} (unmodified below this header).`]) + body
    : body;
  downloadText(`mra_${page.version}_${filename.replace(/\//g, "_")}`, stamped,
               filename.endsWith(".csv") ? "text/csv" : "application/json");
}

const assets = Object.entries(release.assets || {});
const rows = assets.map(([filename, meta]) => ({
  filename,
  description: DESCRIPTIONS[filename] ||
    (filename.startsWith("aa_paths/")
      ? `Codon-to-codon path table for ${filename.replace("aa_paths/", "").replace(".csv", "")}: every one of the 63 alternative codons per codon.`
      : ""),
  gzip_kb: Math.round(meta.gzip_bytes / 1024),
  sha256: meta.sha256,
}));

renderTable(document.getElementById("asset-table"), {
  columns: [
    { key: "filename", label: "File" },
    { key: "description", label: "What it is" },
    { key: "gzip_kb", label: "KB (gz)", type: "num" },
    { key: "sha256", label: "SHA-256",
      render: (r) => `<span class="mono" title="${r.sha256}">${r.sha256.slice(0, 12)}\u2026</span>` },
    { key: "filename", label: "",
      render: (r) => {
        const btn = document.createElement("button");
        btn.className = "secondary";
        btn.textContent = "Download";
        btn.addEventListener("click", async () => {
          btn.disabled = true; btn.textContent = "checking\u2026";
          try {
            const check = await verifyAssetChecksum(page.version, r.filename, release);
            if (check.status === "MISMATCH") {
              btn.textContent = "checksum mismatch";
              return;
            }
            await fetchAndDownload(r.filename);
            btn.textContent = check.status === "ok" ? "verified \u2713" : "downloaded";
          } catch (err) {
            btn.textContent = "failed: " + err.message;
          } finally {
            setTimeout(() => { btn.disabled = false; btn.textContent = "Download"; }, 2500);
          }
        });
        return btn;
      } },
  ],
  rows: rows.filter((r) => r.filename !== "mra_web_model.json"),
  sort: { key: "filename", dir: "asc" },
  filterPlaceholder: "filter files",
});

const bundleRow = rows.find((r) => r.filename === "mra_web_model.json");
if (bundleRow) {
  const p = document.createElement("p");
  const btn = document.createElement("button");
  btn.textContent = `Download mra_web_model.json (${bundleRow.gzip_kb} KB gzipped)`;
  btn.addEventListener("click", () => fetchAndDownload("mra_web_model.json"));
  p.appendChild(btn);
  document.getElementById("bundle-slot").appendChild(p);
  document.getElementById("bundle-slot").appendChild(caveat(
    `SHA-256 <code>${bundleRow.sha256}</code>. The bundle carries its own ` +
    `<code>version</code> field; the portal refuses to load one whose embedded version ` +
    `disagrees with the directory it came from.`));
}

document.getElementById("dictionary").appendChild(note(
  'The full column-by-column dictionary is <a href="docs/data_dictionary.md">' +
  "docs/data_dictionary.md</a> in the repository. The essentials: rates are " +
  "<strong>per site per replication</strong>, control-corrected; <code>bio_rate_median</code> " +
  "is the posterior median and is the column the amino-acid layer is built from; " +
  "<code>bio_rate_mean</code> is the posterior mean and is <em>not</em> interchangeable with it, " +
  "because the mean of a rate is not the exponential of the mean of the log rate; " +
  "<code>mu_aa</code> is summed over every path to that amino acid, including 2- and " +
  "3-substitution routes; rates in the shipped CSVs are rounded to 6 significant figures."));

const pre = document.createElement("pre");
pre.className = "citation";
pre.textContent = citationText(release);
document.getElementById("citation-slot").appendChild(pre);
const copy = document.createElement("button");
copy.className = "secondary";
copy.textContent = "Copy citation";
copy.addEventListener("click", async () => {
  try { await navigator.clipboard.writeText(citationText(release)); copy.textContent = "Copied"; }
  catch { copy.textContent = "Select the text above to copy"; }
  setTimeout(() => { copy.textContent = "Copy citation"; }, 2500);
});
document.getElementById("citation-slot").appendChild(copy);
document.getElementById("citation-slot").appendChild(caveat(
  "Atlas data are licensed CC-BY-4.0, so attribution is a licence condition, not a courtesy. " +
  "Portal code is MIT."));
