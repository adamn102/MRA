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
import { codonSensitivity } from "./codon_sensitivity.js";

// The five columns every table on this page carries. Defined once so the drug-resistance and
// antibody-escape tables cannot drift apart, and so a reader comparing a DRM to a bNAb escape
// is comparing like with like. `mutationKey` is the only thing that differs: the DRM tables
// name a mutation in protein numbering, the bNAb tables an escape change in HXB2 numbering.
const RATE_COLUMNS = (mutationKey) => [
  { key: mutationKey, label: "Mutation" },
  { key: "mu_aa", label: "Amino Acid Mutability (Mutation Rate)", type: "rate",
    title: "mu_aa: rate of reaching this amino acid from the NL4-3 codon, summed over every " +
           "nucleotide path, per site per replication" },
  { key: "interval", label: "95% Confidence Interval" },
  { key: "min_nt_changes", label: "Minimum Nucleotides", type: "num",
    title: "Fewest nucleotide substitutions needed to reach this amino acid" },
  { key: "n_paths", label: "Paths", type: "num",
    title: "Number of distinct nucleotide paths from this codon to this amino acid" },
];

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

  const bnabNote =
    "<strong>Curated from primary literature.</strong> Every change listed was reported to " +
    "affect neutralisation \u2014 escape selected in vivo, a rebound virus, or a reported " +
    "resistance assignment. Structural contacts and glycans without a measured effect are " +
    "excluded. The antibody, the evidence strength (primary, secondary or clinical) and the " +
    "citation for each change are in <code>data/" + page.version + "/bnab_catalog.csv</code> " +
    "and <code>docs/BNAB_PROVENANCE.md</code>. Escape is published in HXB2 numbering against " +
    "varied backgrounds; changes whose published wild-type residue NL4-3 does not carry have " +
    "no atlas rate and are counted under each table rather than listed. One change per row: " +
    "a combined genotype rate would assume sequential accumulation, which recombination can " +
    "bypass.";

  for (const section of TABLE_SECTIONS) {
    if (section.kind === "drm") {
      const cls = drm.taxonomy.classes.find((c) => c.class === section.key);
      if (!cls) continue;                       // a class the shipped catalog does not carry
      renderDrmSection({ page, engine, drm, host, section, cls });
      rendered.push(section.id);
    } else if (section.kind === "bnab") {
      const records = bnab.records.filter((r) => r.epitope_class === section.key);
      if (!records.length) continue;
      renderBnabSection({ page, engine, host, section, records,
                          intro: rendered.some((id) => id.startsWith("bnab-")) ? null : bnabNote });
      rendered.push(section.id);
    }
  }
  return rendered;
}

// --- one drug class ---------------------------------------------------------------------
function renderDrmSection({ page, engine, drm, host, section, cls }) {
  // Every column here is a property of the mutation and the atlas, so the table is fixed:
  // there is nothing left for a control to vary. The treatment-stratum and accessibility
  // selectors went with the clinical-frequency columns they fed.
  const inClass = drm.records.filter((r) => r.drug_class === section.key);
  const byMut = new Map();
  for (const r of inClass) {
    if (!byMut.has(r.mutation)) byMut.set(r.mutation, {});
    byMut.get(r.mutation)[r.rx] = r;
  }

  const rows = [];
  let indelsExcluded = 0;
  for (const [mutation, byRx] of byMut) {
    // One record per mutation is enough now: the treatment strata differed only in the
    // clinical-frequency columns, and the rate is the same in every one of them.
    const any = byRx.naive || Object.values(byRx)[0];
    // Indel-coded DRMs (Stanford's ~ and # codes) are real resistance mechanisms but are not
    // amino-acid substitutions, so the model estimates no mutability for them. They are
    // excluded rather than listed rate-less; the count is reported under the table so the
    // omission is visible instead of silent.
    if (any.mechanism !== "substitution") { indelsExcluded++; continue; }
    rows.push({
      mutation,
      hxb2: any.protein + " " + any.ref_aa + any.protein_position + any.alt_aa,
      mu_aa: any.mu_aa,
      interval: `${formatRate(any.mu_aa_lo95)} \u2013 ${formatRate(any.mu_aa_hi95)}`,
      min_nt_changes: any.min_nt_changes,
      n_paths: any.n_paths,
    });
  }

  const el = sectionShell(host, section,
    `${cls.label}. ${rows.length} substitution mutations, including 1-, 2- and ` +
    `3-nucleotide changes.`);

  const row = document.createElement("div");
  row.className = "controls";
  el.appendChild(row);
  const shareSlot = document.createElement("div");
  row.appendChild(shareSlot);
  shareSlot.appendChild(shareButton({ version: page.version, extra: { class: section.key } }));

  const tableSlot = document.createElement("div");
  el.appendChild(tableSlot);

  renderTable(tableSlot, {
    columns: RATE_COLUMNS("hxb2"),
    rows: rows.sort((a, b) => (b.mu_aa || 0) - (a.mu_aa || 0)),
    sort: { key: "mu_aa", dir: "desc" },
    filterPlaceholder: "filter mutations",
    downloadName: `mra_${page.version}_DRM_${section.key}.csv`,
    provenance: exportProvenance(page.release, [
      `Drug-resistance mutations, class ${section.key}.`,
      "Rates are mutational supply from the NL4-3 atlas, per site per replication.",
      "1-, 2- and 3-nucleotide changes are all included.",
      "Indel-coded DRMs are excluded: they are not amino-acid substitutions.",
    ]),
  });

  if (indelsExcluded) {
    tableSlot.appendChild(caveat(
      `${indelsExcluded} further ${section.key} DRM(s) in the Stanford catalog are in-frame ` +
      `codon insertions or deletions. They are real resistance mechanisms, but they are not ` +
      `amino-acid substitutions, so the model estimates no mutability for them and they are ` +
      `not listed here.`));
  }
}

// --- one bNAb epitope class -------------------------------------------------------------
function renderBnabSection({ page, engine, host, section, records, intro }) {
  // Collapsed to one row per escape change. The catalog carries one row per antibody, so a
  // position targeted by two antibodies appeared twice; with the antibody column gone those
  // would read as a duplicated row. The rate is a property of the codon, not of the antibody,
  // so the duplicates carry the same number and collapsing them loses no rate information.
  const byChange = new Map();
  for (const r of records) {
    const change = `${r.wt_aa}${r.hxb2_position}${r.escape_aa}`;
    if (!byChange.has(change)) byChange.set(change, r);
  }

  const rows = [];
  let notInNl43 = 0;
  for (const [change, r] of byChange) {
    // Where NL4-3 does not carry the published wild-type residue the atlas has no rate for
    // the change, because its rates are all from the NL4-3 codon. Such rows are counted and
    // reported under the table rather than listed with an empty rate.
    if (!r.applies_to_nl43 || !r.nl43_env_position) { notInNl43++; continue; }
    // Recomputed live from the engine rather than read from the CSV column, so this table and
    // the codon-sensitivity engine can never disagree; test/bnab_cases.mjs asserts they match.
    let mu = null, minnt = null;
    const s2 = codonSensitivity(engine, { gene: "Env", genePosition: r.nl43_env_position,
      refAa: r.wt_aa, altAa: r.escape_aa });
    if (s2.applicable && s2.nl43_rate) { mu = s2.nl43_rate; minnt = s2.min_nt_changes; }
    if (mu === null) { mu = r.mu_aa; minnt = r.min_nt_changes; }
    if (mu === null || mu === undefined) { notInNl43++; continue; }
    rows.push({
      change,
      mu_aa: mu,
      interval: (r.mu_aa_lo95 && r.mu_aa_hi95)
        ? `${formatRate(r.mu_aa_lo95)} \u2013 ${formatRate(r.mu_aa_hi95)}` : "\u2014",
      min_nt_changes: minnt,
      n_paths: r.n_paths !== undefined && r.n_paths !== null && r.n_paths !== ""
        ? Number(r.n_paths) : null,
    });
  }

  const el = sectionShell(host, section, null);
  if (intro) el.appendChild(note(intro));
  const tableSlot = document.createElement("div");
  el.appendChild(tableSlot);

  renderTable(tableSlot, {
    columns: RATE_COLUMNS("change"),
    rows: rows.sort((a, b) => (b.mu_aa || 0) - (a.mu_aa || 0)),
    sort: { key: "mu_aa", dir: "desc" },
    filterPlaceholder: false,
    downloadName: `mra_${page.version}_bNAb_escape_${section.key}.csv`,
    provenance: exportProvenance(page.release, [
      `bNAb escape, epitope class ${section.key}. Curated from primary literature; every`,
      "change listed was reported to affect neutralisation. Per-row antibody, evidence",
      "strength and citation are in data/<version>/bnab_catalog.csv and",
      "docs/BNAB_PROVENANCE.md.",
      "One change per row; no combined multi-mutation rates (recombination not modelled).",
    ]),
  });

  if (notInNl43) {
    tableSlot.appendChild(caveat(
      `${notInNl43} further curated ${section.key} escape change(s) are not listed: NL4-3 ` +
      `does not carry the published wild-type residue at those positions, so the atlas has ` +
      `no rate for the change. They remain in bnab_catalog.csv with their citations.`));
  }
}
