// SPDX-License-Identifier: MIT
//
// The curated bNAb escape catalog. Run:
//   node test/bnab_cases.mjs --version v0.1
//
// Two things must hold. First, the rate shipped in the CSV must equal what the engine computes
// live from the NL4-3 codon -- the table recomputes it, so a disagreement would mean the
// published file and the page disagree. Second, the evidence and applicability rules the
// catalog claims must actually hold row by row, because those rules are what make it citable.

import fs from "node:fs";
import { MraEngine } from "../js/mra_engine.js";
import { alignToReference } from "../js/align.js";
import { parseCsv } from "../js/data.js";

const version = process.argv.includes("--version")
  ? process.argv[process.argv.indexOf("--version") + 1] : "v0.1";
const D = `data/${version}/`;
const J = (f) => JSON.parse(fs.readFileSync(D + f, "utf8"));

let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

// --- load the catalog -------------------------------------------------------------
// Parsed with the site's own parser, so the test exercises the same quoting rules the
// browser will -- note fields contain commas and would defeat a naive split.
const parsed = parseCsv(fs.readFileSync(D + "bnab_catalog.csv", "utf8"));
const head = parsed.columns, rows = parsed.records;
console.log(`bNAb catalog: ${rows.length} rows, ${new Set(rows.map((r) => r.antibody)).size} antibodies\n`);

const ref = J("nl43_reference.json");
const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference: ref,
  intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });
const out = engine.customAtlas({ name: "NL4-3", sequence: ref.sequence,
  alignment: alignToReference(ref.sequence, ref.sequence) });
const idx = new Map();
for (const r of out.aa) idx.set(`${r.protein}|${r.position}|${r.alt_aa}`, r);

// --- 1. schema and required fields --------------------------------------------------
const REQUIRED = ["antibody","epitope_class","hxb2_position","wt_aa","escape_aa",
                  "applies_to_nl43","evidence_type","evidence_strength","citation","citation_doi"];
check("every required column is present",
      REQUIRED.every((c) => head.includes(c)),
      REQUIRED.filter((c) => !head.includes(c)).join(",") || "all present");

check("every row carries a citation and a resolvable-looking DOI",
      rows.every((r) => r.citation && /^10\.\d{4,}\//.test(r.citation_doi)),
      rows.filter((r) => !/^10\.\d{4,}\//.test(r.citation_doi)).map((r) => r.antibody).join(",") || "all cited");

const STRENGTHS = new Set(["primary", "secondary", "clinical"]);
check("evidence_strength uses only the three agreed tiers",
      rows.every((r) => STRENGTHS.has(r.evidence_strength)),
      [...new Set(rows.map((r) => r.evidence_strength))].join(","));

check("no row is left flagged provisional",
      !head.includes("provisional") && !rows.some((r) => String(r.provisional) === "True"));

// --- 2. the NL4-3 applicability rule ------------------------------------------------
// A row that claims to apply must have a rate; a row that does not must have NO rate and must
// say why. This is the rule that keeps a misleading number off the page.
const applies = rows.filter((r) => r.applies_to_nl43 === "True");
const notApplies = rows.filter((r) => r.applies_to_nl43 !== "True");
check("rows that apply to NL4-3 all carry a rate",
      applies.every((r) => Number(r.mu_aa) > 0), `${applies.length} applicable`);
check("rows that do not apply carry NO rate",
      notApplies.every((r) => !r.mu_aa), `${notApplies.length} not applicable`);
check("every non-applicable row explains why in its note",
      notApplies.every((r) => /NOT APPLICABLE/i.test(r.note || "")),
      notApplies.map((r) => `${r.antibody} ${r.wt_aa}${r.hxb2_position}${r.escape_aa}`).join("; "));

// --- 3. the shipped rate equals what the engine computes live ------------------------
let mismatches = [];
for (const r of applies) {
  const nl = engine.nl43PositionFromHxb2("Env", Number(r.hxb2_position));
  const rec = idx.get(`Env|${nl}|${r.escape_aa}`);
  if (!rec) { mismatches.push(`${r.antibody} ${r.wt_aa}${r.hxb2_position}${r.escape_aa}: no atlas row`); continue; }
  const rel = Math.abs(rec.mu_aa - Number(r.mu_aa)) / rec.mu_aa;
  if (rel > 1e-3) mismatches.push(`${r.wt_aa}${r.hxb2_position}${r.escape_aa}: csv ${r.mu_aa} vs engine ${rec.mu_aa.toExponential(4)}`);
  if (Number(r.min_nt_changes) !== rec.min_nt_changes)
    mismatches.push(`${r.wt_aa}${r.hxb2_position}${r.escape_aa}: steps csv ${r.min_nt_changes} vs engine ${rec.min_nt_changes}`);
}
check("every shipped rate matches the engine recomputation", mismatches.length === 0,
      mismatches.slice(0, 3).join(" | ") || `${applies.length} rows agree`);

// --- 4. NL4-3 really does carry the published wild type where claimed ----------------
const CODON = {TTT:"F",TTC:"F",TTA:"L",TTG:"L",CTT:"L",CTC:"L",CTA:"L",CTG:"L",ATT:"I",ATC:"I",
 ATA:"I",ATG:"M",GTT:"V",GTC:"V",GTA:"V",GTG:"V",TCT:"S",TCC:"S",TCA:"S",TCG:"S",CCT:"P",CCC:"P",
 CCA:"P",CCG:"P",ACT:"T",ACC:"T",ACA:"T",ACG:"T",GCT:"A",GCC:"A",GCA:"A",GCG:"A",TAT:"Y",TAC:"Y",
 CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",AAT:"N",AAC:"N",AAA:"K",AAG:"K",GAT:"D",GAC:"D",GAA:"E",GAG:"E",
 TGT:"C",TGC:"C",TGG:"W",CGT:"R",CGC:"R",CGA:"R",CGG:"R",AGT:"S",AGC:"S",AGA:"R",AGG:"R",GGT:"G",
 GGC:"G",GGA:"G",GGG:"G"};
const Env = ref.proteins.Env;
const aaAt = (nl) => CODON[(Env.codon_seq[nl - Env.first_aa] || "").toUpperCase()];
const wrongWt = applies.filter((r) => {
  const nl = engine.nl43PositionFromHxb2("Env", Number(r.hxb2_position));
  return aaAt(nl) !== r.wt_aa;
});
check("NL4-3 carries the published wild type on every applicable row",
      wrongWt.length === 0,
      wrongWt.map((r) => `${r.wt_aa}${r.hxb2_position}: NL4-3 has ${aaAt(engine.nl43PositionFromHxb2("Env", Number(r.hxb2_position)))}`).join("; "));

// --- 5. the supply claim the catalog is for ------------------------------------------
// A281K escapes N6 as effectively as A281D but needs three nucleotide changes. If the atlas
// ever stopped reproducing that gap, the central example in the write-up would be wrong.
const a281d = applies.find((r) => r.antibody === "N6" && r.escape_aa === "D");
const a281k = applies.find((r) => r.antibody === "N6" && r.escape_aa === "K");
check("N6 A281D is a single nucleotide step", Number(a281d.min_nt_changes) === 1);
check("N6 A281K needs three nucleotide steps, as the source reports",
      Number(a281k.min_nt_changes) === 3, `${a281k.min_nt_changes}`);
check("the one-step N6 escape is orders of magnitude more accessible",
      Number(a281d.mu_aa) / Number(a281k.mu_aa) > 1e6,
      `${(Number(a281d.mu_aa) / Number(a281k.mu_aa)).toExponential(1)}-fold`);

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL bNAb CATALOG CASES PASSED");
