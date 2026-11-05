// SPDX-License-Identifier: MIT
//
// drm.js -- the resistance-and-escape page.
//
// Two catalogs, one engine. The DRM half is a database join: every number is read from
// drm_catalog.csv, which was rebuilt from the Stanford SDRD exports and joined to the
// published atlas. The bNAb half is curated from primary literature, every row citing a
// study that measured an effect on neutralisation. Clicking any
// row runs the codon-sensitivity engine live, in the browser, on the small model bundle.
//
// The page is careful about three distinctions that are easy to blur:
//   - mutational SUPPLY (the rate the variant is generated) vs clinical FREQUENCY (how
//     often it is seen, which is supply filtered through selection). Both columns shown.
//   - codon OPTIMIZATION (rate varies among synonymous codons at equal accessibility) vs
//     a change in the NUMBER of nucleotide steps (categorical). The drill-down separates
//     them, as the engine does.
//   - a specific escape change vs glycan SEQUON LOSS (summed over three codons). Glycan
//     bNAb rows carry both.

import { bootstrapPage, note, caveat, showError, exportProvenance, shareButton } from "./page.js";
import { loadSmallAssets, loadDrmCatalog, loadBnabCatalog, loadLanlCodonFrequencies } from "./data.js";
import { renderTable, formatRate } from "./table.js";
import { barChart } from "./chart.js";
import { MraEngine } from "./mra_engine.js";
import { codonSensitivity, codonSensitivitySummary, formatFold, sequonLossRate }
  from "./codon_sensitivity.js";

const page = await bootstrapPage("drm.html");
const region = page.release.modelled_region;
document.getElementById("scope-slot").appendChild(note(
  "Rates are the control-corrected biological mutation rate per site per replication, from " +
  "the NL4-3 atlas (positions " + region.first_position.toLocaleString() + "\u2013" +
  region.last_position.toLocaleString() + "). Resistance positions are reported in standard " +
  "protein numbering (PR, RT, IN, CA). Clinical frequencies are Stanford SDR-DB sequence " +
  "counts; <code>freq = (mutants + 0.5) / sequences</code>."));

const status = document.getElementById("drm-status");
let assets, engine, drm, bnab, lanl;
try {
  assets = await loadSmallAssets(page.version);
  engine = new MraEngine({ bundle: assets.bundle, reference: assets.reference,
                           intervals: assets.intervals, hxb2: assets.hxb2 });
  [drm, bnab] = await Promise.all([loadDrmCatalog(page.version), loadBnabCatalog(page.version)]);
  lanl = await loadLanlCodonFrequencies(page.version);   // null when not shipped
  status.textContent = `${new Set(drm.records.map((r) => r.mutation)).size} DRMs across ` +
    `${drm.taxonomy.classes.length} drug classes loaded.`;
} catch (err) {
  showError(status, "Could not load the catalog: " + err.message);
  throw err;
}

// --- drug class / stratum selectors ----------------------------------------------
const classSel = document.getElementById("drug-class");
const stratSel = document.getElementById("drug-stratum");
for (const c of drm.taxonomy.classes) {
  const opt = document.createElement("option");
  opt.value = c.class;
  opt.textContent = `${c.class} \u2014 ${c.label} (${c.n_mutations})`;
  classSel.appendChild(opt);
}

function populateStrata(cls) {
  stratSel.innerHTML = "";
  const c = drm.taxonomy.classes.find((x) => x.class === cls);
  for (const s of c.strata) {
    const opt = document.createElement("option");
    opt.value = s.rx;
    opt.textContent = s.rx === "naive" ? "treatment-naive (baseline)"
                    : s.kind === "class" ? `${s.rx} (class, pooled)` : s.rx;
    stratSel.appendChild(opt);
  }
  // default to the first on-treatment stratum, falling back to naive
  const treated = c.strata.find((s) => s.rx !== "naive");
  stratSel.value = treated ? treated.rx : "naive";
}

function renderDrm() {
  const cls = classSel.value;
  const stratum = stratSel.value;
  const accessAll = document.getElementById("access-filter").value === "all";

  // one row per unique mutation; clinical frequency from the chosen stratum and from naive
  const inClass = drm.records.filter((r) => r.drug_class === cls);
  const byMut = new Map();
  for (const r of inClass) {
    if (!byMut.has(r.mutation)) byMut.set(r.mutation, {});
    byMut.get(r.mutation)[r.rx] = r;
  }
  const rows = [];
  for (const [mutation, strata] of byMut) {
    const any = strata[stratum] || strata.naive || Object.values(strata)[0];
    if (any.mechanism !== "substitution") {
      // indel-coded DRM: real, but no mu_aa. Show it, labelled, never with a fake rate.
      rows.push({
        mutation, hxb2: any.protein + " " + any.ref_aa + any.protein_position + any.alt_aa,
        alt_aa: any.alt_aa, mechanism: any.mechanism, mu_aa: null, interval: "\u2014",
        min_nt_changes: null, n_paths: null,
        freq_treated: strata[stratum] ? strata[stratum].freq : null,
        freq_naive: strata.naive ? strata.naive.freq : null,
        nl43_codon: "\u2014", note: any.rate_note || "in-frame codon indel; not a substitution",
        _row: any, _indel: true,
      });
      continue;
    }
    if (!accessAll && any.min_nt_changes !== 1) continue;
    rows.push({
      mutation,
      hxb2: any.protein + " " + any.ref_aa + any.protein_position + any.alt_aa,
      alt_aa: any.alt_aa, mechanism: "substitution",
      mu_aa: any.mu_aa, interval: `${formatRate(any.mu_aa_lo95)} \u2013 ${formatRate(any.mu_aa_hi95)}`,
      min_nt_changes: any.min_nt_changes, n_paths: any.n_paths,
      freq_treated: strata[stratum] ? strata[stratum].freq : null,
      freq_naive: strata.naive ? strata.naive.freq : null,
      nl43_codon: any.nl43_codon,
      note: any.nl43_matches_consensus ? "" : `NL4-3 carries ${any.atlas_ref_aa}, not consensus ${any.ref_aa}`,
      ...codonFields(codonSensitivitySummary(engine, {
        gene: any.gene, genePosition: Number(any.gene_position),
        refAa: any.ref_aa, altAa: any.alt_aa })),
      _row: any,
    });
  }

  const label = (s) => s === "naive" ? "naive" : s;
  renderTable(document.getElementById("drm-table"), {
    columns: [
      { key: "hxb2", label: "Mutation" },
      { key: "mu_aa", label: "mutation rate", type: "rate",
        title: "mu_aa: summed mutability to this amino acid, per replication" },
      { key: "interval", label: "95% interval" },
      { key: "min_nt_changes", label: "min nt", type: "num" },
      { key: "n_paths", label: "paths", type: "num" },
      { key: "freq_treated", label: `freq (${label(stratum)})`, type: "rate",
        title: "clinical sequence frequency in this treatment stratum" },
      { key: "freq_naive", label: "freq (naive)", type: "rate" },
      ...CODON_COLUMNS,
      { key: "note", label: "Note" },
    ],
    rows: rows.sort((a, b) => (b.mu_aa || 0) - (a.mu_aa || 0)),
    sort: { key: "mu_aa", dir: "desc" },
    filterPlaceholder: "filter mutations",
    downloadName: `mra_${page.version}_DRM_${cls}_${stratum}.csv`,
    provenance: exportProvenance(page.release, [
      `Drug-resistance mutations, class ${cls}, clinical frequency stratum ${stratum}.`,
      "Rates from the NL4-3 atlas; frequencies from Stanford SDR-DB.",
      "Indel-coded DRMs are listed without a rate: they are not substitutions.",
    ]),
  });

  const indels = rows.filter((r) => r._indel);
  if (indels.length) {
    document.getElementById("drm-table").appendChild(caveat(
      `${indels.length} listed ${cls} DRM(s) are in-frame codon insertions or deletions ` +
      `(Stanford's ~ and # codes). They are real resistance mechanisms but are not amino-acid ` +
      `substitutions, so the model estimates no mutability for them and none is shown.`));
  }
}

classSel.addEventListener("change", () => { populateStrata(classSel.value); renderDrm(); refreshShare(); });
stratSel.addEventListener("change", () => { renderDrm(); refreshShare(); });
document.getElementById("access-filter").addEventListener("change", renderDrm);

function refreshShare() {
  const slot = document.getElementById("drm-share");
  slot.innerHTML = "";
  slot.appendChild(shareButton({ version: page.version,
    extra: { class: classSel.value, rx: stratSel.value } }));
}

// --- bNAb selector ----------------------------------------------------------------
document.getElementById("bnab-note").appendChild(note(
  "<strong>Curated from primary literature.</strong> Each row is a substitution with a cited " +
  "study showing a neutralisation effect \u2014 escape selected in vivo, a rebound virus, or a " +
  "reported resistance assignment. Structural contacts and glycans without a measured effect are " +
  "excluded. <strong>Evidence</strong> is <em>primary</em> (reduces neutralisation on its own), " +
  "<em>secondary</em> (only with a partner mutation) or <em>clinical</em> (emerged in a treated " +
  "participant). Escape is published in HXB2 numbering against varied backgrounds; where NL4-3 " +
  "lacks the published wild-type residue the row is kept with its citation but no rate, since the " +
  "atlas rates only changes from the NL4-3 codon. One mutation per row: a combined genotype rate " +
  "would assume sequential accumulation, which recombination can bypass."));

const bnabSel = document.getElementById("bnab-class");
const bnabClasses = [...new Set(bnab.records.map((r) => r.epitope_class))].sort();
for (const c of bnabClasses) {
  const opt = document.createElement("option");
  opt.value = c; opt.textContent = c;
  bnabSel.appendChild(opt);
}

function renderBnab() {
  const cls = bnabSel.value;
  const inClass = bnab.records.filter((r) => r.epitope_class === cls);
  const sequonByHxb2 = new Map(bnab.sequons.map((s) => [s.hxb2_position, s]));
  const rows = [];
  for (const r of inClass) {
    // The rate is recomputed live from the engine rather than read from the CSV column, so
    // the table and the codon-sensitivity drill-down can never disagree. The CSV value is
    // the same number; test/bnab_cases.mjs asserts the two agree.
    let mu = null, lo = null, hi = null, minnt = null;
    if (r.applies_to_nl43 && r.nl43_env_position && r.escape_aa) {
      const s2 = codonSensitivity(engine, { gene: "Env", genePosition: r.nl43_env_position,
        refAa: r.wt_aa, altAa: r.escape_aa });
      if (s2.applicable && s2.nl43_rate) { mu = s2.nl43_rate; minnt = s2.min_nt_changes; }
    }
    if (mu === null && r.applies_to_nl43) { mu = r.mu_aa; minnt = r.min_nt_changes; }
    if (r.applies_to_nl43) { lo = r.mu_aa_lo95; hi = r.mu_aa_hi95; }
    rows.push({
      antibody: r.antibody,
      change: `${r.wt_aa}${r.hxb2_position}${r.escape_aa}`,
      evidence: r.evidence_strength,
      effect: r.effect,
      mu,
      interval: (lo && hi) ? `${formatRate(lo)} \u2013 ${formatRate(hi)}` : "\u2014",
      steps: minnt,
      n_star: mu ? Math.round(1 / mu) : null,
      nl43: r.applies_to_nl43 ? `Env ${r.nl43_env_position}` : "\u2014",
      citation: r.citation,
      doi: r.citation_doi,
      note: r.note || "",
      ...codonFields(r.applies_to_nl43 && r.nl43_env_position
        ? codonSensitivitySummary(engine, { gene: "Env",
            genePosition: Number(r.nl43_env_position), refAa: r.wt_aa, altAa: r.escape_aa })
        : null),
      // Several escape substitutions work by destroying an N-linked glycan sequon. The total
      // supply for losing that glycan by ANY route is larger than this one substitution --
      // three codons can destroy it -- so it belongs beside the single-substitution rate
      // rather than one click away, where a reader comparing the two would never see it.
      ...sequonFields(sequonByHxb2.get(r.hxb2_position) || null, r),
      _row: r,
      _sequon: sequonByHxb2.get(r.hxb2_position) || null,
    });
  }

  renderTable(document.getElementById("bnab-table"), {
    columns: [
      { key: "antibody", label: "Antibody" },
      { key: "change", label: "Escape (HXB2)" },
      { key: "evidence", label: "Evidence" },
      { key: "mu", label: "escape rate", type: "rate",
        title: "mu_aa for this substitution from the NL4-3 codon; blank where NL4-3 " +
               "does not carry the published wild-type residue" },
      { key: "interval", label: "95% interval" },
      { key: "steps", label: "nt steps" },
      { key: "n_star", label: "N* for one copy",
        title: "replication events for one expected de novo copy, 1/rate" },
      { key: "nl43", label: "NL4-3" },
      { key: "citation", label: "Source",
        render: (r) => {
          const a = document.createElement("a");
          a.href = "https://doi.org/" + r.doi;
          a.textContent = r.citation;
          a.target = "_blank"; a.rel = "noopener noreferrer";
          return a;
        } },
      ...CODON_COLUMNS,
      { key: "sequon_loss", label: "glycan loss (any route)", type: "rate",
        title: "Where this position carries an N-X-S/T sequon, the summed rate of destroying " +
               "it by any single substitution in any of its three codons. Larger than the " +
               "single escape rate beside it, because several changes destroy the motif." },
      { key: "effect", label: "Reported effect" },
      { key: "note", label: "Note" },
    ],
    rows: rows.sort((a, b) => (b.mu || 0) - (a.mu || 0)),
    sort: { key: "mu", dir: "desc" },
    filterPlaceholder: false,
    downloadName: `mra_${page.version}_bNAb_escape_${cls}.csv`,
    provenance: exportProvenance(page.release, [
      `bNAb escape, epitope class ${cls}. Curated from primary literature; every row cites a`,
      "study that measured an effect on neutralisation. See docs/BNAB_PROVENANCE.md.",
      "Rows without a rate: NL4-3 does not carry the published wild-type residue.",
      "One mutation per row; no combined multi-mutation rates (recombination not modelled).",
    ]),
  });
}
bnabSel.addEventListener("change", renderBnab);

// --- codon sensitivity, inline ----------------------------------------------------
// Every row carries its own codon decomposition rather than hiding it behind a control. The
// whole catalogue costs ~14 ms to compute, so there is nothing to defer, and the comparison a
// reader wants -- is this mutation rate a property of the virus or of how this virus happens
// to spell the wild-type residue? -- is only visible when it sits next to the rate itself.

/** Flatten a summary into row fields, so the table and its CSV export carry the same values. */
function codonFields(sum) {
  if (!sum) return { ref_codons: "", codon_fold: null, codon_fold_text: "\u2014",
                     n_ref_codons: null, nl43_codon_rank: "" };
  return {
    ref_codons: sum.text,
    codon_fold: sum.n_codons > 1 ? sum.fold_range : null,
    codon_fold_text: formatFold(sum.fold_range, sum.n_codons),
    n_ref_codons: sum.n_codons,
    nl43_codon_rank: sum.n_codons > 1 && sum.nl43_rank
      ? `${sum.nl43_rank} of ${sum.n_codons}` : "",
  };
}

/** Glycan sequon-loss total for a row, where its position carries a sequon. */
function sequonFields(sequon, row) {
  if (!sequon || !row.applies_to_nl43) return { sequon_loss: null, sequon_motif: "" };
  const loss = sequonLossRate(engine, sequon);
  if (!loss || !loss.applicable) return { sequon_loss: null, sequon_motif: "" };
  return { sequon_loss: loss.mu_sequon_loss, sequon_motif: loss.motif };
}

const CODON_COLUMNS = [
  { key: "ref_codons", label: "Wild-type codon options (rate)",
    title: "Every synonymous codon the wild-type residue could be written as, and the rate " +
           "this change would have from each. * marks the codon NL4-3 carries; (2nt)/(3nt) " +
           "marks a codon needing more than one substitution to reach the target.",
    render: (r) => {
      const span = document.createElement("span");
      span.className = "codon-options";
      span.textContent = r.ref_codons || "\u2014";
      return span;
    } },
  { key: "codon_fold", label: "codon sensitivity", type: "num",
    title: "Fold range (max/min) of this change's rate across every synonymous codon for the " +
           "wild-type residue. 1x or blank means the residue has a single codon, so codon " +
           "usage cannot move this rate. Large values are usually driven by a codon that " +
           "needs an extra nucleotide, which the codon list shows.",
    render: (r) => r.codon_fold_text },
];


// --- open on URL params -----------------------------------------------------------
const qs = new URLSearchParams(window.location.search);
populateStrata(classSel.value);
if (qs.get("class") && drm.taxonomy.classes.some((c) => c.class === qs.get("class"))) {
  classSel.value = qs.get("class");
  populateStrata(classSel.value);
  if (qs.get("rx")) stratSel.value = qs.get("rx");
}
renderDrm(); renderBnab(); refreshShare();
