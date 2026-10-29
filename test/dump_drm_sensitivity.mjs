// SPDX-License-Identifier: MIT
// dump_drm_sensitivity.mjs -- run codonSensitivity over the whole DRM catalog so the
// figure and any analysis use the same engine the portal uses. Writes
// test/out/drm_sensitivity.json. Substitution DRMs only; indel-coded rows have no codon.
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { MraEngine } from "../js/mra_engine.js";
import { codonSensitivity } from "../js/codon_sensitivity.js";
import { parseCsv } from "../js/data.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const DATA = path.join(HERE, "..", "data", "v0.1");
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));
const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference: J("nl43_reference.json"),
                               intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });

const cat = parseCsv(fs.readFileSync(path.join(DATA, "drm_catalog.csv"), "utf8")).records;
const seen = new Map();
for (const r of cat) {
  if (r.mechanism !== "substitution") continue;
  const key = r.mutation + "|" + r.alt_aa;
  if (seen.has(key)) continue;
  const s = codonSensitivity(engine, {
    gene: r.gene, genePosition: Number(r.gene_position), refAa: r.ref_aa, altAa: r.alt_aa });
  if (!s.applicable) continue;
  seen.set(key, {
    mutation: r.mutation, drug_class: r.drug_class, aa_change: r.aa_change,
    min_nt_changes: Number(r.min_nt_changes),
    nl43_matches_consensus: r.nl43_matches_consensus === "True",
    nl43Codon: s.nl43Codon, n_syn: s.n_synonymous_codons,
    min_steps: s.min_nt_changes, n_at_min: s.n_codons_at_min, n_higher: s.n_codons_higher_step,
    site_fold: s.site_fold_range, full_fold: s.full_fold_range, step_count_varies: s.step_count_varies,
    flank_fold: s.flank_varied ? s.flank_fold_range : null,
    nl43_rate: s.nl43_rate, nl43_rank: s.nl43_rank_among_peers, n_peers: s.n_peers_at_nl43_step,
    site_min: s.site_min, site_max: s.site_max, flank_coupled: s.flank_coupled,
  });
}
const out = [...seen.values()];
fs.mkdirSync(path.join(HERE, "out"), { recursive: true });
fs.writeFileSync(path.join(HERE, "out", "drm_sensitivity.json"), JSON.stringify(out, null, 1));
const multi = out.filter((r) => r.n_at_min > 1);
console.log(`${out.length} substitution DRMs; ${multi.length} have >1 synonymous WT codon at the minimal step count`);
const big = out.filter((r) => r.site_fold >= 2).sort((a, b) => b.site_fold - a.site_fold);
console.log(`${big.length} have a >=2-fold CODON-OPTIMIZATION range (equal accessibility). Top 10:`);
for (const r of big.slice(0, 10))
  console.log(`  ${r.mutation.padEnd(10)} ${r.drug_class.padEnd(5)} ${r.nl43Codon} ` +
              `codon-opt ${r.site_fold.toFixed(2)}x  rank ${r.nl43_rank}/${r.n_peers}`);
const stepv = out.filter((r) => r.step_count_varies);
console.log(`\n${stepv.length} DRMs have synonymous WT codons in DIFFERENT accessibility classes ` +
            `(a change in the minimum number of nt steps).`);
