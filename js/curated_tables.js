// SPDX-License-Identifier: MIT
//
// curated_tables.js -- the resistance and antibody-escape half of the MRA Table page.
//
// Two catalogs, one engine. The DRM half is a database join: every number is read from
// drm_catalog.csv, which was rebuilt from the Stanford SDRD exports and joined to the
// published atlas. The bNAb half is curated from primary literature, every row citing a study
// that measured an effect on neutralisation. Clicking any row runs the codon-sensitivity
// engine live, in the browser, on the small model bundle.
//
// WHY EVERY TABLE IS RENDERED AT ONCE. This was one table behind a drug-class selector, which
// made the catalog look like one table that changes rather than the eight it is, and left no
// address for "the NNRTI table" to be linked to or cited. Each class now gets its own
// anchored section, built from the TABLE_SECTIONS list in page.js so that the nav dropdown,
// the index at the top of the page and the anchors here cannot drift apart. The whole
// catalogue costs ~14 ms of codon-sensitivity computation, so there is nothing to defer.
//
// Treatment stratum stays a per-section control rather than a page-wide one, because the
// strata are class-specific: the PI table's strata are PI regimens and mean nothing under
// NNRTI.
//
// The page is careful about three distinctions that are easy to blur:
//   - mutational SUPPLY (the rate the variant is generated) vs clinical FREQUENCY (how often
//     it is seen, which is supply filtered through selection). Both columns shown.
//   - codon OPTIMIZATION (rate varies among synonymous codons at equal accessibility) vs a
//     change in the NUMBER of nucleotide steps (categorical). The drill-down separates them.
//   - a specific escape change vs glycan SEQUON LOSS (summed over three codons). Glycan bNAb
//     rows carry both.

import { note, caveat, exportProvenance, shareButton, TABLE_SECTIONS } from "./page.js";
import { renderTable, formatRate } from "./table.js";
import { codonSensitivity, codonSensitivitySummary, formatFold, sequonLossRate }
  from "./codon_sensitivity.js";

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

/** Build the <section> shell for one entry in TABLE_SECTIONS. */
function sectionShell(host, section, subtitle) {
  const el = document.createElement("section");
  el.className = "table-section";
  el.id = section.id;
  const h = document.createElement("h2");
  h.textContent = section.label;
  el.appendChild(h);
  if (subtitle) {
    const p = document.createElement("p");
    p.className = "caveat";
    p.innerHTML = subtitle;
    el.appendChild(p);
  }
  host.appendChild(el);
  return el;
}

/** A labelled <select>, returned with its element so the caller can read and watch it. */
function controlSelect(row, labelText, options, initial) {
  const wrap = document.createElement("div");
  const label = document.createElement("label");
  label.textContent = labelText;
  const sel = document.createElement("select");
  for (const o of options) {
    const opt = document.createElement("option");
    opt.value = o.value; opt.textContent = o.label;
    sel.appendChild(opt);
  }
  if (initial !== undefined) sel.value = initial;
  label.htmlFor = sel.id = "sel-" + Math.random().toString(36).slice(2, 10);
  wrap.appendChild(label); wrap.appendChild(sel);
  row.appendChild(wrap);
  return sel;
}

/**
 * Render the eight curated sections into `host`, in TABLE_SECTIONS order.
 *
 * Returns the list of section ids actually rendered, so the caller can build an index that
 * cannot name a section the page does not have.
 */
export function initCuratedTables({ page, engine, drm, bnab, host }) {
  const rendered = [];
  const sequonByHxb2 = new Map(bnab.sequons.map((s) => [s.hxb2_position, s]));

  const bnabNote =
    "<strong>Curated from primary literature.</strong> Each row is a substitution with a " +
    "cited study showing a neutralisation effect \u2014 escape selected in vivo, a rebound " +
    "virus, or a reported resistance assignment. Structural contacts and glycans without a " +
    "measured effect are excluded. <strong>Evidence</strong> is <em>primary</em> (reduces " +
    "neutralisation on its own), <em>secondary</em> (only with a partner mutation) or " +
    "<em>clinical</em> (emerged in a treated participant). Escape is published in HXB2 " +
    "numbering against varied backgrounds; where NL4-3 lacks the published wild-type residue " +
    "the row is kept with its citation but no rate. One mutation per row: a combined genotype " +
    "rate would assume sequential accumulation, which recombination can bypass.";

  for (const section of TABLE_SECTIONS) {
    if (section.kind === "drm") {
      const cls = drm.taxonomy.classes.find((c) => c.class === section.key);
      if (!cls) continue;                       // a class the shipped catalog does not carry
      renderDrmSection({ page, engine, drm, host, section, cls });
      rendered.push(section.id);
    } else if (section.kind === "bnab") {
      const records = bnab.records.filter((r) => r.epitope_class === section.key);
      if (!records.length) continue;
      renderBnabSection({ page, engine, bnab, host, section, records, sequonByHxb2,
                          intro: rendered.some((id) => id.startsWith("bnab-")) ? null : bnabNote });
      rendered.push(section.id);
    }
  }
  return rendered;
}

// --- one drug class ---------------------------------------------------------------------
function renderDrmSection({ page, engine, drm, host, section, cls }) {
  const el = sectionShell(host, section,
    `${cls.label}. ${cls.n_mutations} documented mutations in the Stanford catalog.`);

  const row = document.createElement("div");
  row.className = "controls";
  el.appendChild(row);

  const strata = cls.strata || [];
  const stratSel = controlSelect(row, "Treatment stratum (for clinical frequency)",
    strata.map((s) => ({ value: s.rx,
      label: s.rx === "naive" ? "treatment-naive (baseline)"
           : s.kind === "class" ? `${s.rx} (class, pooled)` : s.rx })));
  const treated = strata.find((s) => s.rx !== "naive");
  stratSel.value = treated ? treated.rx : "naive";

  const accessSel = controlSelect(row, "Accessibility", [
    { value: "1", label: "1-nucleotide changes only" },
    { value: "all", label: "include 2- and 3-nucleotide changes" },
  ], "1");

  const shareSlot = document.createElement("div");
  row.appendChild(shareSlot);

  const tableSlot = document.createElement("div");
  el.appendChild(tableSlot);

  function draw() {
    const stratum = stratSel.value;
    const accessAll = accessSel.value === "all";

    // one row per unique mutation; clinical frequency from the chosen stratum and from naive
    const inClass = drm.records.filter((r) => r.drug_class === section.key);
    const byMut = new Map();
    for (const r of inClass) {
      if (!byMut.has(r.mutation)) byMut.set(r.mutation, {});
      byMut.get(r.mutation)[r.rx] = r;
    }
    const rows = [];
    for (const [mutation, byRx] of byMut) {
      const any = byRx[stratum] || byRx.naive || Object.values(byRx)[0];
      if (any.mechanism !== "substitution") {
        // indel-coded DRM: real, but no mu_aa. Show it, labelled, never with a fake rate.
        rows.push({
          mutation, hxb2: any.protein + " " + any.ref_aa + any.protein_position + any.alt_aa,
          alt_aa: any.alt_aa, mechanism: any.mechanism, mu_aa: null, interval: "\u2014",
          min_nt_changes: null, n_paths: null,
          freq_treated: byRx[stratum] ? byRx[stratum].freq : null,
          freq_naive: byRx.naive ? byRx.naive.freq : null,
          nl43_codon: "\u2014",
          note: any.rate_note || "in-frame codon indel; not a substitution",
          ...codonFields(null),
          _row: any, _indel: true,
        });
        continue;
      }
      if (!accessAll && any.min_nt_changes !== 1) continue;
      rows.push({
        mutation,
        hxb2: any.protein + " " + any.ref_aa + any.protein_position + any.alt_aa,
        alt_aa: any.alt_aa, mechanism: "substitution",
        mu_aa: any.mu_aa,
        interval: `${formatRate(any.mu_aa_lo95)} \u2013 ${formatRate(any.mu_aa_hi95)}`,
        min_nt_changes: any.min_nt_changes, n_paths: any.n_paths,
        freq_treated: byRx[stratum] ? byRx[stratum].freq : null,
        freq_naive: byRx.naive ? byRx.naive.freq : null,
        nl43_codon: any.nl43_codon,
        note: any.nl43_matches_consensus
          ? "" : `NL4-3 carries ${any.atlas_ref_aa}, not consensus ${any.ref_aa}`,
        ...codonFields(codonSensitivitySummary(engine, {
          gene: any.gene, genePosition: Number(any.gene_position),
          refAa: any.ref_aa, altAa: any.alt_aa })),
        _row: any,
      });
    }

    const label = (s) => s === "naive" ? "naive" : s;
    renderTable(tableSlot, {
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
      downloadName: `mra_${page.version}_DRM_${section.key}_${stratum}.csv`,
      provenance: exportProvenance(page.release, [
        `Drug-resistance mutations, class ${section.key}, clinical frequency stratum ${stratum}.`,
        "Rates from the NL4-3 atlas; frequencies from Stanford SDR-DB.",
        "Indel-coded DRMs are listed without a rate: they are not substitutions.",
      ]),
    });

    const indels = rows.filter((r) => r._indel);
    if (indels.length) {
      tableSlot.appendChild(caveat(
        `${indels.length} listed ${section.key} DRM(s) are in-frame codon insertions or ` +
        `deletions (Stanford's ~ and # codes). They are real resistance mechanisms but are ` +
        `not amino-acid substitutions, so the model estimates no mutability for them and ` +
        `none is shown.`));
    }

    shareSlot.innerHTML = "";
    shareSlot.appendChild(shareButton({ version: page.version,
      extra: { class: section.key, rx: stratum } }));
  }

  stratSel.addEventListener("change", draw);
  accessSel.addEventListener("change", draw);
  draw();
}

// --- one bNAb epitope class -------------------------------------------------------------
function renderBnabSection({ page, engine, host, section, records, sequonByHxb2, intro }) {
  const el = sectionShell(host, section, null);
  if (intro) el.appendChild(note(intro));
  const tableSlot = document.createElement("div");
  el.appendChild(tableSlot);

  const rows = [];
  for (const r of records) {
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
    const sequon = sequonByHxb2.get(r.hxb2_position) || null;
    let sequonLoss = { sequon_loss: null, sequon_motif: "" };
    if (sequon && r.applies_to_nl43) {
      const loss = sequonLossRate(engine, sequon);
      if (loss && loss.applicable) {
        sequonLoss = { sequon_loss: loss.mu_sequon_loss, sequon_motif: loss.motif };
      }
    }
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
      ...sequonLoss,
      _row: r,
      _sequon: sequon,
    });
  }

  renderTable(tableSlot, {
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
    downloadName: `mra_${page.version}_bNAb_escape_${section.key}.csv`,
    provenance: exportProvenance(page.release, [
      `bNAb escape, epitope class ${section.key}. Curated from primary literature; every row`,
      "cites a study that measured an effect on neutralisation. See docs/BNAB_PROVENANCE.md.",
      "Rows without a rate: NL4-3 does not carry the published wild-type residue.",
      "One mutation per row; no combined multi-mutation rates (recombination not modelled).",
    ]),
  });
}
