// SPDX-License-Identifier: MIT
//
// oneStepAminoAcids, the genetic-code helper that used to feed the codon-difference
// table. That table is no longer on Your Sequence. Run:
//   node test/compare_core_cases.mjs
//
import { GENETIC_CODE } from "../js/mra_engine.js";
import { oneStepAminoAcids } from "../js/explore_core.js";

let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

// ATG (Met) reaches these in a single substitution; its own residue is excluded.
const atg = oneStepAminoAcids("ATG", GENETIC_CODE);
check("oneStepAminoAcids excludes the codon's own residue", !atg.has("M"));
check("oneStepAminoAcids(ATG) includes I, L, V, K, T, R", ["I","L","V","K","T","R"].every((a) => atg.has(a)),
      [...atg].sort().join(""));

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL COMPARE-CORE CASES PASSED");
