// SPDX-License-Identifier: MIT
//
// data.js -- loading and indexing the shipped atlas tables.
//
// PAYLOAD DISCIPLINE
//   The model bundle is ~170 KB gzipped and the browse tables are several megabytes.
//   A user who only wants the sequence tool must not pay for the browse tables, so
//   nothing here is loaded until a view actually asks for it, and each table is fetched
//   at most once per page. loadSmallAssets() is the sequence tool's whole requirement.
//
// Every published number the portal displays comes from one of these files. Nothing
// published is recomputed in JavaScript.

import { assetUrl, loadModelBundle } from "./version.js";

const cache = new Map();

async function once(key, fn) {
  if (!cache.has(key)) cache.set(key, fn());
  return cache.get(key);
}

async function getJson(version, name) {
  return once(version + "/" + name, async () => {
    const resp = await fetch(assetUrl(version, name));
    if (!resp.ok) throw new Error("Could not load " + name + " (HTTP " + resp.status + ")");
    return resp.json();
  });
}

// --- a small, correct CSV parser ------------------------------------------------
// Our own files have no embedded commas, but an exported file a user re-uploads might,
// and a parser that splits on "," is the kind of thing that fails silently years later.
export function parseCsv(text) {
  const rows = [];
  let row = [], field = "", inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (inQuotes) {
      if (ch === '"') {
        if (text[i + 1] === '"') { field += '"'; i++; } else inQuotes = false;
      } else field += ch;
    } else if (ch === '"') inQuotes = true;
    else if (ch === ",") { row.push(field); field = ""; }
    else if (ch === "\n") { row.push(field); field = ""; rows.push(row); row = []; }
    else if (ch !== "\r") field += ch;
  }
  if (field !== "" || row.length) { row.push(field); rows.push(row); }
  if (!rows.length) return { columns: [], records: [] };
  const columns = rows[0];
  const records = [];
  for (let r = 1; r < rows.length; r++) {
    if (rows[r].length === 1 && rows[r][0] === "") continue;     // trailing blank line
    const rec = {};
    for (let c = 0; c < columns.length; c++) rec[columns[c]] = rows[r][c];
    records.push(rec);
  }
  return { columns, records };
}

async function getCsv(version, name, coerce) {
  return once(version + "/" + name, async () => {
    const resp = await fetch(assetUrl(version, name));
    if (!resp.ok) throw new Error("Could not load " + name + " (HTTP " + resp.status + ")");
    const { columns, records } = parseCsv(await resp.text());
    if (coerce) for (const rec of records) coerce(rec);
    return { columns, records };
  });
}

const num = (v) => (v === "" || v === undefined ? null : Number(v));

// --------------------------------------------------------------------------------
/** Everything the sequence tool needs, and nothing else. ~360 KB gzipped in total. */
export async function loadSmallAssets(version) {
  const [bundle, reference, intervals, hxb2, citation] = await Promise.all([
    once(version + "/bundle", () => loadModelBundle(version)),
    getJson(version, "nl43_reference.json"),
    getJson(version, "interval_ratios.json"),
    getJson(version, "hxb2_align.json"),
    getJson(version, "citation.json"),
  ]);
  return { bundle, reference, intervals, hxb2, citation };
}

/** Published nucleotide layer, indexed by reference position. */
export async function loadNucleotideLayer(version) {
  return once(version + "/nt-index", async () => {
    const { records } = await getCsv(version, "nt_level_supply.csv", (r) => {
      r.ref_position = Number(r.ref_position);
      r.bio_rate_median = num(r.bio_rate_median);
      r.bio_rate_lo95 = num(r.bio_rate_lo95);
      r.bio_rate_hi95 = num(r.bio_rate_hi95);
      r.bio_rate_mean = num(r.bio_rate_mean);
      r.wt_count = num(r.wt_count);
      r.wt_coverage = num(r.wt_coverage);
    });
    const byPosition = new Map();
    for (const r of records) {
      let list = byPosition.get(r.ref_position);
      if (!list) byPosition.set(r.ref_position, (list = []));
      list.push(r);
    }
    return { records, byPosition };
  });
}

/** Published amino-acid layer, indexed by protein and position. */
export async function loadAminoAcidLayer(version) {
  return once(version + "/aa-index", async () => {
    const { records } = await getCsv(version, "aa_level_supply.csv", (r) => {
      r.position = Number(r.position);
      r.mu_aa = num(r.mu_aa);
      r.mu_aa_lo95 = num(r.mu_aa_lo95);
      r.mu_aa_hi95 = num(r.mu_aa_hi95);
      r.min_nt_changes = Number(r.min_nt_changes);
      r.n_paths = Number(r.n_paths);
      r.synonymous = r.synonymous === "True" || r.synonymous === "true";
    });
    const byPosition = new Map();           // "Protein|position" -> rows
    const proteins = new Map();             // protein -> {min, max}
    for (const r of records) {
      const key = r.protein + "|" + r.position;
      let list = byPosition.get(key);
      if (!list) byPosition.set(key, (list = []));
      list.push(r);
      const span = proteins.get(r.protein);
      if (!span) proteins.set(r.protein, { min: r.position, max: r.position });
      else { span.min = Math.min(span.min, r.position); span.max = Math.max(span.max, r.position); }
    }
    return { records, byPosition, proteins };
  });
}

/** The DRM catalog and its selector taxonomy. ~0.1 MB gzipped. */
export async function loadDrmCatalog(version) {
  return once(version + "/drm", async () => {
    const [cat, taxonomy] = await Promise.all([
      getCsv(version, "drm_catalog.csv", (r) => {
        r.protein_position = Number(r.protein_position);
        r.gene_position = Number(r.gene_position);
        r.mu_aa = num(r.mu_aa); r.mu_aa_lo95 = num(r.mu_aa_lo95); r.mu_aa_hi95 = num(r.mu_aa_hi95);
        r.min_nt_changes = num(r.min_nt_changes); r.n_paths = num(r.n_paths);
        r.all_seqs = num(r.all_seqs); r.all_muts = num(r.all_muts); r.freq = num(r.freq);
        r.nl43_matches_consensus = r.nl43_matches_consensus === "True";
      }),
      getJson(version, "drm_taxonomy.json"),
    ]);
    return { records: cat.records, taxonomy };
  });
}

/** Optional LANL codon-frequency weighting. Returns null when the asset is absent --
 *  the portal is fully functional without it, and the UI shows an em dash instead of a
 *  silent zero. Keyed "gene|protein_position|codon" -> {fraction, n_sequences, version}. */
export async function loadLanlCodonFrequencies(version) {
  return once(version + "/lanl", async () => {
    try {
      const resp = await fetch(assetUrl(version, "lanl_codon_frequencies.csv"));
      if (!resp.ok) return null;                 // 404: asset not shipped in this release
      const { records } = parseCsv(await resp.text());
      const byKey = new Map();
      let alignmentVersion = null;
      for (const r of records) {
        byKey.set(r.gene + "|" + Number(r.protein_position) + "|" + r.codon,
                  { fraction: num(r.lanl_fraction), n: num(r.n_sequences) });
        alignmentVersion = r.alignment_version;
      }
      return { byKey, alignmentVersion, present: true };
    } catch {
      return null;
    }
  });
}

/** The curated bNAb escape catalog and the glycan sequon list. */
export async function loadBnabCatalog(version) {
  return once(version + "/bnab", async () => {
    const [cat, sequons] = await Promise.all([
      getCsv(version, "bnab_catalog.csv", (r) => {
        r.hxb2_position = Number(r.hxb2_position);
        r.nl43_env_position = num(r.nl43_env_position);
        // applies_to_nl43 is False where NL4-3 does not carry the published wild-type residue.
        // Those rows keep their annotation and citation but must not show a rate: the atlas
        // can only rate a change that starts from the NL4-3 codon.
        r.applies_to_nl43 = String(r.applies_to_nl43) === "True";
        r.mu_aa = num(r.mu_aa);
        r.mu_aa_lo95 = num(r.mu_aa_lo95);
        r.mu_aa_hi95 = num(r.mu_aa_hi95);
        r.min_nt_changes = num(r.min_nt_changes);
        r.n_star = num(r.n_star);
      }),
      getJson(version, "glycan_sequons.json"),
    ]);
    return { records: cat.records, sequons: sequons.sequons };
  });
}

/** One protein's codon-path shard. Loaded only when a dashboard asks for that protein. */
export async function loadPathShard(version, protein) {
  const { records } = await getCsv(version, "aa_paths/" + protein + ".csv", (r) => {
    r.position = Number(r.position);
    r.n_nt_changes = Number(r.n_nt_changes);
    r.mu_path = num(r.mu_path);
    r.mu_path_lo95 = num(r.mu_path_lo95);
    r.mu_path_hi95 = num(r.mu_path_hi95);
  });
  const byPosition = new Map();
  for (const r of records) {
    let list = byPosition.get(r.position);
    if (!list) byPosition.set(r.position, (list = []));
    list.push(r);
  }
  return { records, byPosition };
}
