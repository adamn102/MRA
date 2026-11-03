// SPDX-License-Identifier: MIT
//
// The "compare to NL4-3" analyses now folded into the Your-sequence detail view. Run:
//   node test/compare_core_cases.mjs --version v0.1
//
// These were a separate tab ("Compare & genome track") before the merge. The point of the
// cases is that the baseline is the reference's OWN custom atlas, not the published CSV: both
// sides of the fold change then come from the same engine, so an unchanged codon gives exactly
// zero rather than the published file's storage-rounding noise.

import fs from "node:fs";
import { MraEngine, GENETIC_CODE } from "../js/mra_engine.js";
import { alignToReference } from "../js/align.js";
import { codonDifferences, accessibilityRows, oneStepAminoAcids } from "../js/explore_core.js";

const version = process.argv.includes("--version")
  ? process.argv[process.argv.indexOf("--version") + 1] : "v0.1";
const D = `data/${version}/`;
const J = (f) => JSON.parse(fs.readFileSync(D + f, "utf8"));

let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

const ref = J("nl43_reference.json");
const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference: ref,
  intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });
const seq = ref.sequence;

const refAtlas = engine.customAtlas({ name: "NL4-3", sequence: seq,
  alignment: alignToReference(seq, seq) });
const refByPos = new Map();
for (const r of refAtlas.aa) {
  const k = r.protein + "|" + r.position;
  let list = refByPos.get(k); if (!list) refByPos.set(k, (list = [])); list.push(r);
}

// --- 1. oneStepAminoAcids ------------------------------------------------------------
// ATG (Met) reaches these in a single substitution; its own residue is excluded.
const atg = oneStepAminoAcids("ATG", GENETIC_CODE);
check("oneStepAminoAcids excludes the codon's own residue", !atg.has("M"));
check("oneStepAminoAcids(ATG) includes I, L, V, K, T, R", ["I","L","V","K","T","R"].every((a) => atg.has(a)),
      [...atg].sort().join(""));

// --- 2. NL4-3 against itself: nothing differs ----------------------------------------
const cdSelf = codonDifferences(refAtlas.aa, ref.proteins, GENETIC_CODE);
check("NL4-3 vs itself has no codon differences", cdSelf.length === 0, `${cdSelf.length}`);
const accSelf = accessibilityRows(refAtlas.aa, refByPos, ref.proteins, GENETIC_CODE);
const nonzero = accSelf.filter((r) => r.log2fc !== null && Math.abs(r.log2fc) > 1e-12);
check("NL4-3 vs itself has exactly zero fold change everywhere", nonzero.length === 0,
      `${nonzero.length} nonzero of ${accSelf.length}`);

// --- 3. a planted mutation is detected at its position -------------------------------
// Flip NL4-3 Env codon 460 (nt 6804 is its third base) and confirm the codon difference and
// the reshaped one-step reachability come out of the genetic code.
const a = seq.split("");
a[6804 - 1] = "G";                                   // AAT (N) -> AGT (S) at Env 460
const mut = a.join("");
const out = engine.customAtlas({ name: "mut", sequence: mut, alignment: alignToReference(mut, seq) });
const cd = codonDifferences(out.aa, ref.proteins, GENETIC_CODE);
// position is NL4-3 internal numbering (HXB2 460 == NL4-3 458 in Env); the hxb2_label carries
// the clinical number. Key on the codon change itself, which is unambiguous.
const env460 = cd.find((r) => r.protein === "Env" && r.ref_codon === "AAT" && r.your_codon === "AGT");
check("a planted Env codon change (AAT->AGT) is reported", !!env460,
      cd.map((r) => `${r.hxb2_label} ${r.ref_codon}->${r.your_codon}`).join(", "));
check("the HXB2 label names Env 460", env460 && /460/.test(env460.hxb2_label || ""),
      env460 ? env460.hxb2_label : "n/a");
check("it is called an amino-acid difference (N -> S)",
      env460 && env460.ref_aa === "N" && env460.your_aa === "S" && env460.kind === "amino-acid difference");
// gained/lost must match the genetic code directly
const gainTruth = [...oneStepAminoAcids("AGT", GENETIC_CODE)].filter((x) => !oneStepAminoAcids("AAT", GENETIC_CODE).has(x)).sort().join(" ");
check("gained reachability matches the genetic code", env460 && env460.gained === gainTruth,
      env460 ? `${env460.gained} vs ${gainTruth}` : "n/a");

// --- 4. accessibility picks up the shift at the mutated codon ------------------------
const acc = accessibilityRows(out.aa, refByPos, ref.proteins, GENETIC_CODE);
const shifted = acc.filter((r) => r.log2fc !== null && Math.abs(r.log2fc) > 1e-9);
check("the planted change produces accessibility shifts", shifted.length > 0, `${shifted.length} rows`);
// every shifted row must sit at a codon that actually differs (codon_differs true)
check("every shifted accessibility row is at a differing codon",
      shifted.every((r) => r.codon_differs), `${shifted.filter((r) => !r.codon_differs).length} at identical codons`);
// a row absent from the NL4-3 table is only possible where the codon differs
const novel = acc.filter((r) => r.mu_nl43 === null);
check("rows absent from the NL4-3 table only occur at differing codons",
      novel.every((r) => r.codon_differs), `${novel.filter((r) => !r.codon_differs).length} at identical codons`);

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL COMPARE-CORE CASES PASSED");
