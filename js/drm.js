// SPDX-License-Identifier: MIT
//
// drm.js -- the resistance-and-escape page.
//
// Two catalogs, one engine. The DRM half is a database join: every number is read from
// drm_catalog.csv, which was rebuilt from the Stanford SDRD exports and joined to the
// published atlas. The bNAb half is PROVISIONAL and marked so on every row. Clicking any
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
import { codonSensitivity, sequonLossRate } from "./codon_sensitivity.js";

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
      { key: "nl43_codon", label: "NL4-3 codon" },
      { key: "note", label: "Note" },
      { key: "mutation", label: "",
        render: (r) => {
          if (r._indel) return "\u2014";
          const b = document.createElement("button");
          b.className = "secondary"; b.textContent = "codon sensitivity";
          b.addEventListener("click", () => drillDrm(r._row));
          return b;
        } },
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
document.getElementById("bnab-provisional").appendChild(note(
  "<strong>Provisional.</strong> These antibody-escape annotations are a literature draft, " +
  "not a database export, and are not yet signed off. Positions and glycan dependencies are " +
  "the well-known ones, but every row is marked provisional and the citations are not yet " +
  "filled in. See <code>docs/BNAB_PROVENANCE.md</code>. The mutation rates themselves are " +
  "computed by the same validated engine as everything else."));

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
    const base = {
      antibody: r.antibody, mechanism: r.mechanism,
      position: "Env " + r.hxb2_position + (r.wt_aa ? " (" + r.wt_aa + ")" : ""),
      nl43_ok: r.nl43_matches_expected,
      note: r.nl43_matches_expected ? "" : `NL4-3 has ${r.nl43_wt_aa}, not ${r.wt_aa}`,
      _row: r,
    };
    if (r.mechanism === "glycan_sequon") {
      const sq = sequonByHxb2.get(r.hxb2_position);
      const loss = sq ? sequonLossRate(engine, sq) : null;
      rows.push({ ...base,
        change: `lose N${r.hxb2_position} glycan`,
        mu: loss ? loss.mu_sequon_loss : null,
        interval: loss ? `${formatRate(loss.lo95)} \u2013 ${formatRate(loss.hi95)}` : "\u2014",
        kind: "sequon loss (3 codons)",
        _sequon: sq, _loss: loss });
    } else {
      // specific contact change: compute its mu_aa directly if NL4-3 carries the WT
      const nl = r.nl43_env_position;
      let mu = null, interval = "\u2014", minnt = null;
      if (nl && r.escape_aa) {
        const s = codonSensitivity(engine, { gene: "Env", genePosition: nl,
          refAa: r.wt_aa, altAa: r.escape_aa });
        if (s.applicable && s.nl43_rate) { mu = s.nl43_rate; minnt = s.min_nt_changes; }
      }
      rows.push({ ...base,
        change: r.wt_aa + r.hxb2_position + (r.escape_aa || "?"),
        mu, interval, kind: "specific change",
        min_nt_changes: minnt });
    }
  }

  renderTable(document.getElementById("bnab-table"), {
    columns: [
      { key: "antibody", label: "Antibody" },
      { key: "change", label: "Escape" },
      { key: "kind", label: "Kind" },
      { key: "mu", label: "escape rate", type: "rate",
        title: "specific change: mu_aa; sequon loss: summed over all three codons" },
      { key: "interval", label: "95% interval" },
      { key: "position", label: "Env (HXB2)" },
      { key: "note", label: "Note" },
      { key: "antibody", label: "",
        render: (r) => {
          const b = document.createElement("button");
          b.className = "secondary"; b.textContent = "codon sensitivity";
          b.addEventListener("click", () => drillBnab(r._row, r._sequon, r._loss));
          return b;
        } },
    ],
    rows: rows.sort((a, b) => (b.mu || 0) - (a.mu || 0)),
    sort: { key: "mu", dir: "desc" },
    filterPlaceholder: false,
    downloadName: `mra_${page.version}_bNAb_${cls}_PROVISIONAL.csv`,
    provenance: exportProvenance(page.release, [
      `PROVISIONAL bNAb escape, epitope class ${cls}. Not signed off; see BNAB_PROVENANCE.md.`,
      "Glycan rows report summed sequon-loss rate over three codons.",
    ]),
  });
}
bnabSel.addEventListener("change", renderBnab);

// --- codon-sensitivity drill-down -------------------------------------------------
function drillDrm(row) {
  const s = codonSensitivity(engine, { gene: row.gene, genePosition: row.gene_position,
    refAa: row.ref_aa, altAa: row.alt_aa });
  showDrill(`${row.protein} ${row.ref_aa}${row.protein_position}${row.alt_aa}`, s, null);
}
function drillBnab(row, sequon, loss) {
  if (row.mechanism === "glycan_sequon" && loss) {
    showSequonDrill(row, loss);
  } else if (row.nl43_env_position && row.escape_aa) {
    const s = codonSensitivity(engine, { gene: "Env", genePosition: row.nl43_env_position,
      refAa: row.wt_aa, altAa: row.escape_aa });
    showDrill(`${row.antibody}: Env ${row.wt_aa}${row.hxb2_position}${row.escape_aa}`, s, null);
  }
}

function showDrill(title, s, _) {
  document.getElementById("drill-heading").hidden = false;
  const intro = document.getElementById("drill-intro");
  intro.hidden = false;
  const chart = document.getElementById("drill-chart");
  const detail = document.getElementById("drill-detail");
  chart.innerHTML = ""; detail.innerHTML = "";
  if (!s.applicable) { intro.textContent = s.reason; return; }

  intro.innerHTML =
    `<strong>${title}.</strong> The rate of this change if the wild-type amino acid were ` +
    `encoded by each of its synonymous codons, holding the neighbouring amino acid fixed. ` +
    `NL4-3 uses <code>${s.nl43Codon}</code>` +
    (s.nl43_rank_among_peers ? `, which ranks ${s.nl43_rank_among_peers} of ${s.n_peers_at_nl43_step} ` +
      `(1 = lowest supply) among codons reaching the target in ${s.nl43_step_count} step(s).` : ".") +
    ` Codon-optimization fold range: <strong>${s.site_fold_range.toFixed(1)}\u00d7</strong>` +
    (s.step_count_varies
      ? `. Some synonymous codons need a different number of nucleotide changes (a categorical ` +
        `accessibility difference); including those the full range is ${s.full_fold_range.toExponential(1)}\u00d7.`
      : ".") +
    (s.flank_varied
      ? ` Varying the 3\u2032 flanking codon synonymously changes the rate by up to ` +
        `${s.flank_fold_range.toFixed(2)}\u00d7.` : "");

  // Optional LANL weighting: how common each synonymous codon is in circulating virus.
  const lanlFor = (codon) => {
    if (!lanl || !lanl.present) return null;
    const hit = lanl.byKey.get(`${s.gene}|${s.genePosition}|${codon}`);
    return hit ? hit.fraction : 0;
  };
  const anyLanl = lanl && lanl.present && s.by_site_codon.some((r) => lanlFor(r.codon) !== null);
  barChart(chart, s.by_site_codon.map((r) => {
    const f = lanlFor(r.codon);
    const steps = r.min_nt_changes > s.min_nt_changes ? "more steps" : "";
    const freq = anyLanl ? ` ${((f || 0) * 100).toFixed(0)}% LANL` : "";
    return {
      label: r.codon + (r.is_nl43 ? " (NL4-3)" : "") + (r.min_nt_changes > s.min_nt_changes ? ` [${r.min_nt_changes}nt]` : ""),
      value: r.mu_aa, lo: r.lo95, hi: r.hi95, highlight: r.is_nl43,
      note: [steps, freq].filter(Boolean).join("  "),
    };
  }), { title: "mutation rate per synonymous wild-type codon", width: 720,
        xLabel: "mu_aa (per site per replication)" });
  if (anyLanl) {
    // LANL-weighted mean rate and the fraction of circulating virus on a high-supply codon
    const median = s.by_site_codon.map((r) => r.mu_aa).sort((a, b) => a - b)[Math.floor(s.by_site_codon.length / 2)];
    let wMean = 0, wSum = 0, highFrac = 0;
    for (const r of s.by_site_codon) {
      const f = lanlFor(r.codon) || 0;
      wMean += f * r.mu_aa; wSum += f;
      if (r.mu_aa >= median) highFrac += f;
    }
    intro.innerHTML += ` <strong>LANL weighting (${lanl.alignmentVersion || "unversioned"}):</strong> ` +
      `circulating-sequence-weighted mean rate ${wSum > 0 ? formatRate(wMean / wSum) : "\u2014"}; ` +
      `${(highFrac * 100).toFixed(0)}% of sequences use an at-or-above-median-supply codon.`;
  } else {
    intro.innerHTML += ` <span class="caveat">LANL codon-frequency weighting is not loaded ` +
      `for this release, so how common each codon is in circulating virus is not shown ` +
      `(the fold range above is the achievable range, not the realised one).</span>`;
  }
  document.getElementById("drill-caption").innerHTML =
    "Highlighted bar is NL4-3's codon. Codons tagged with a step count reach the target amino " +
    "acid only via more nucleotide changes \u2014 a categorical accessibility difference, not codon " +
    "optimization, so they are not in the fold range above.";
  document.getElementById("drill-heading").scrollIntoView({ behavior: "smooth", block: "nearest" });
}

function showSequonDrill(row, loss) {
  document.getElementById("drill-heading").hidden = false;
  const intro = document.getElementById("drill-intro");
  intro.hidden = false;
  const chart = document.getElementById("drill-chart");
  const detail = document.getElementById("drill-detail");
  chart.innerHTML = ""; detail.innerHTML = "";
  intro.innerHTML =
    `<strong>${row.antibody}: loss of the N${row.hxb2_position} glycan (${loss.motif}).</strong> ` +
    `Any change that destroys the N-X-S/T sequon removes the glycan, so the escape rate is the ` +
    `summed mutability over all three codons: <strong>${formatRate(loss.mu_sequon_loss)}</strong>. ` +
    `The single largest route (${loss.largest_contributor.role} \u2192 ${loss.largest_contributor.alt_aa}) ` +
    `is ${(loss.largest_fraction * 100).toFixed(0)}% of that; a view that counted only one change ` +
    `would miss the rest.`;
  barChart(chart, loss.per_codon.map((c) => ({
    label: `${c.role} codon ${c.codon} (${c.ref_aa})`,
    value: c.mu_loss, lo: c.lo95, hi: c.hi95,
    highlight: c.role === loss.largest_contributor.role,
    note: `${c.n_paths} routes`,
  })), { title: `Sequon-loss supply by codon at N${row.hxb2_position}`, width: 720,
         xLabel: "summed mu_aa to any sequon-destroying residue" });
  document.getElementById("drill-caption").textContent =
    "Each bar is the summed rate of all changes at that codon that destroy the glycan sequon. " +
    "N = the asparagine, X = the middle position (loss via proline), S/T = the hydroxyl position.";
}

// --- open on URL params -----------------------------------------------------------
const qs = new URLSearchParams(window.location.search);
populateStrata(classSel.value);
if (qs.get("class") && drm.taxonomy.classes.some((c) => c.class === qs.get("class"))) {
  classSel.value = qs.get("class");
  populateStrata(classSel.value);
  if (qs.get("rx")) stratSel.value = qs.get("rx");
}
renderDrm(); renderBnab(); refreshShare();
