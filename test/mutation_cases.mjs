// SPDX-License-Identifier: MIT
//
// mutation_cases.mjs -- the parsing and coordinate logic behind the table-analysis view
// and every HXB2 label the portal prints.
//
//   node test/mutation_cases.mjs [--version v0.1]
//
// These are the cases where a silent error would be worst: a subunit offset that is off
// by one returns the mutability of the wrong position, with no outward sign that
// anything is wrong. The clinical positions below (RT 103, RT 184, PR 46, PR 90,
// IN 148, IN 155) are ones a virologist will recognise on sight, so a wrong answer here
// is at least visible to the reader.

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { MraEngine } from "../js/mra_engine.js";
import { parseMutationLine, resolveMutation } from "../js/mutations.js";
import { parseCsv } from "../js/data.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const i = process.argv.indexOf("--version");
const VERSION = i >= 0 ? process.argv[i + 1] : "v0.1";
const DATA = path.join(HERE, "..", "data", VERSION);
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));

const engine = new MraEngine({ bundle: J("mra_web_model.json"), reference: J("nl43_reference.json"),
                               intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });

let failures = 0;
const check = (label, cond, detail) => {
  if (cond) console.log(`  OK   ${label}`);
  else { console.log(`  FAIL ${label}${detail ? " -- " + detail : ""}`); failures++; }
};

// --- build the amino-acid index the browse page builds ---------------------------
const aa = parseCsv(fs.readFileSync(path.join(DATA, "aa_level_supply.csv"), "utf8"));
const byPosition = new Map();
const proteins = new Set();
for (const r of aa.records) {
  r.position = Number(r.position);
  r.mu_aa = Number(r.mu_aa);
  r.mu_aa_lo95 = Number(r.mu_aa_lo95);
  r.mu_aa_hi95 = Number(r.mu_aa_hi95);
  r.min_nt_changes = Number(r.min_nt_changes);
  r.n_paths = Number(r.n_paths);
  proteins.add(r.protein);
  const key = r.protein + "|" + r.position;
  let list = byPosition.get(key);
  if (!list) byPosition.set(key, (list = []));
  list.push(r);
}
const ctx = { engine, proteins, lookup: (p, q) => byPosition.get(p + "|" + q) };

console.log(`\n0. the index the browse page builds (${aa.records.length} rows)`);
check("all 62,518 amino-acid rows parsed", aa.records.length === 62518, String(aa.records.length));
check("nine gene products", proteins.size === 9, [...proteins].join(","));

// --- CSV parser ------------------------------------------------------------------
console.log("\n1. CSV parsing");
{
  const { columns, records } = parseCsv('a,b,c\n1,"x,y",3\n4,"say ""hi""",6\n');
  check("header parsed", columns.join("|") === "a|b|c", columns.join("|"));
  check("quoted comma kept inside the field", records[0].b === "x,y", records[0].b);
  check("escaped quotes unescaped", records[1].b === 'say "hi"', records[1].b);
  check("trailing newline does not make a blank record", records.length === 2, String(records.length));
}

// --- the parser -------------------------------------------------------------------
console.log("\n2. mutation-string forms");
{
  const forms = [["RT K65R", "RT", "K", 65, "R"], ["RT 65R", "RT", "", 65, "R"],
                 ["PR:M46I", "PR", "M", 46, "I"], ["  pol m184v  ", "pol", "M", 184, "V"],
                 ["Env N160K", "Env", "N", 160, "K"], ["Gag M1*", "Gag", "M", 1, "*"]];
  for (const [line, name, refAa, pos, altAa] of forms) {
    const p = parseMutationLine(line);
    check(`"${line}"`, p && p.name === name && p.refAa === refAa && p.position === pos && p.altAa === altAa,
          JSON.stringify(p));
  }
  for (const bad of ["RT", "K65R", "RT K65", "", "RT K65RR"]) {
    check(`"${bad}" is rejected rather than guessed at`, parseMutationLine(bad) === null,
          JSON.stringify(parseMutationLine(bad)));
  }
}

// --- subunit offsets ----------------------------------------------------------------
// From the workbook alignment sheet: protease is pol AA 52-150, RT 151-590,
// RNase 591-710, integrase 711-999.
console.log("\n3. subunit offsets and HXB2 labels");
{
  const cases = [["PR", 1, 52], ["PR", 46, 97], ["PR", 90, 141], ["PR", 99, 150],
                 ["RT", 1, 151], ["RT", 65, 215], ["RT", 103, 253], ["RT", 184, 334], ["RT", 440, 590],
                 ["RNH", 1, 591], ["IN", 1, 711], ["IN", 148, 858], ["IN", 155, 865], ["IN", 289, 999]];
  for (const [abbrev, own, expectedPol] of cases) {
    const hit = engine.resolveSubProtein(abbrev, own);
    check(`${abbrev} ${own} -> Pol ${expectedPol}`,
          hit && hit.protein === "Pol" && hit.position === expectedPol, JSON.stringify(hit));
    check(`Pol ${expectedPol} labels back as "${abbrev} ${own}"`,
          engine.hxb2Label("Pol", expectedPol) === `${abbrev} ${own}`,
          engine.hxb2Label("Pol", expectedPol));
  }
  check("Gag subunits resolve too: p24 1 -> Gag 133",
        engine.resolveSubProtein("p24", 1)?.position === 133,
        JSON.stringify(engine.resolveSubProtein("p24", 1)));
}

// --- Env and Vpu need the position-by-position lookup --------------------------------
console.log("\n4. HXB2 numbering for Env and Vpu");
{
  const envMap = engine.nl43ToHxb2.Env;
  check("an Env lookup exists", envMap && envMap.size > 0, String(envMap && envMap.size));
  // The mapping must NOT be a constant offset -- that is exactly the mistake the brief
  // warns about. Check that at least two different offsets occur.
  const offsets = new Set([...envMap.entries()].map(([nl43, hxb2]) => hxb2 - nl43));
  check("Env mapping is not a constant offset", offsets.size > 1, `${offsets.size} distinct offsets`);
  check("Vpu mapping is not a constant offset",
        new Set([...engine.nl43ToHxb2.Vpu.entries()].map(([n, h]) => h - n)).size > 1);
  // Every other protein maps identity.
  for (const p of ["Gag", "Pol", "Vif", "Vpr", "Tat", "Rev", "Nef"]) {
    check(`${p} maps identity`, engine.hxb2Position(p, 42) === 42, String(engine.hxb2Position(p, 42)));
  }
}

// --- resolving against the real atlas -------------------------------------------------
console.log("\n5. resolving real mutations against the atlas");
{
  for (const line of ["RT K103N", "RT M184V", "PR M46I", "PR L90M", "IN N155H", "IN Q148H"]) {
    const r = resolveMutation(line, ctx);
    const ok = !r.error && r.mu_aa > 0 && r.min_nt_changes >= 1 && r.note === "";
    check(`${line} -> ${r.hxb2} ${r.aa_change}, mu_aa ${r.mu_aa ? r.mu_aa.toExponential(2) : "?"}, ` +
          `${r.min_nt_changes} nt, ${r.n_paths} path(s)`, ok, r.error || r.note);
  }
  // A wrong reference amino acid must be reported, not ignored.
  const wrong = resolveMutation("RT A103N", ctx);
  check("a wrong reference amino acid is reported", wrong.note && wrong.note.includes("atlas has K"),
        JSON.stringify(wrong.note));
  // Out of range, with an explanation rather than an empty table.
  const nef = resolveMutation("Nef G100A", ctx);
  check("nef beyond coverage explains itself",
        nef.error && nef.error.includes("beyond the modelled region"), nef.error);
  const bad = resolveMutation("Pol M9999V", ctx);
  check("a position past the protein is refused", !!bad.error, JSON.stringify(bad));
  const unknown = resolveMutation("XYZ M1V", ctx);
  check("an unknown protein is named in the error",
        unknown.error && unknown.error.includes("XYZ"), unknown.error);
}

console.log("\n" + "=".repeat(60));
console.log(failures === 0 ? "ALL MUTATION/COORDINATE CASES PASSED" : `${failures} CASE(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
