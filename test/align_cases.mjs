// SPDX-License-Identifier: MIT
//
// align_cases.mjs -- the alignment cases that matter, with their expected behaviour.
//
//   node test/align_cases.mjs [--version v0.1]
//
// The NL4-3 round-trip covers the identity case exhaustively. These cases cover what a
// real submission looks like: a partial amplicon, an insertion relative to the
// reference, a sequence straddling the end of the modelled region, a sequence entirely
// outside it, and a reverse-complemented paste.
//
// The case that matters most is the insertion. If a large indel shifted the band and
// the aligner silently clipped it, every downstream codon would be read out of frame
// and the portal would return plausible, wrong numbers. So that case asserts not just
// that the insertion is reported, but that the reference positions AFTER it still map
// correctly.

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { alignToReference, parseFasta, reverseComplement } from "../js/align.js";
import { MraEngine } from "../js/mra_engine.js";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, "..");
const i = process.argv.indexOf("--version");
const VERSION = i >= 0 ? process.argv[i + 1] : "v0.1";
const DATA = path.join(REPO, "data", VERSION);
const J = (f) => JSON.parse(fs.readFileSync(path.join(DATA, f), "utf8"));

const bundle = J("mra_web_model.json");
const reference = J("nl43_reference.json");
const engine = new MraEngine({ bundle, reference, intervals: J("interval_ratios.json"), hxb2: J("hxb2_align.json") });
const REF = reference.sequence;
const FIRST = bundle.modelled_region.first_position;
const LAST = bundle.modelled_region.last_position;

let failures = 0;
function check(label, condition, detail) {
  if (condition) console.log(`  OK   ${label}`);
  else { console.log(`  FAIL ${label}${detail ? " -- " + detail : ""}`); failures++; }
}
function section(name) { console.log(`\n${name}`); }

// Deterministic pseudo-random mutagenesis so cases are reproducible run to run.
function mutate(seq, rate, seed) {
  let s = seed >>> 0;
  const rnd = () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
  const bases = ["A", "C", "G", "T"];
  return [...seq].map((b) => (rnd() < rate ? bases[Math.floor(rnd() * 4)] : b)).join("");
}

// ---------------------------------------------------------------------------------
section("1. the reference itself");
{
  const a = alignToReference(REF, REF);
  check("identity is 1", a.identity === 1, `got ${a.identity}`);
  check("every base maps", a.alignedCount === REF.length, `${a.alignedCount}/${REF.length}`);
  check("mapping is the identity permutation",
        [...a.queryToRef].every((v, k) => v === k + 1));
  check("forward orientation chosen", a.orientation === "forward", a.orientation);
}

// ---------------------------------------------------------------------------------
section("2. a clinical pol fragment (PR + start of RT), 5% divergent");
{
  const start = 1311, end = 2600;                 // 1-based reference coordinates
  const frag = mutate(REF.slice(start - 1, end), 0.05, 7);
  const a = alignToReference(frag, REF);
  check("placed at the right reference start", a.refStart === start, `refStart=${a.refStart}`);
  check("placed at the right reference end", a.refEnd === end, `refEnd=${a.refEnd}`);
  check("whole query covered", a.coverageOfQuery === 1, `${a.coverageOfQuery}`);
  check("identity near 0.95", a.identity > 0.9 && a.identity < 1, `${a.identity.toFixed(3)}`);
  const out = engine.customAtlas({ name: "pol fragment", sequence: frag, alignment: a });
  check("nucleotide rows only inside the fragment",
        out.nt.every((r) => r.ref_position >= start && r.ref_position <= end));
  // Reference 1311-1503 is the real gag-pol overlap, so Gag codons appear here too and
  // that is correct, not a leak. HIV-1 genes overlap; the published per-position
  // assignment is what decides, and the engine honours it rather than re-deriving frames.
  const analysed = Object.fromEntries(Object.entries(out.summary.proteins)
    .filter(([, v]) => v.codonsAnalysed > 0).map(([k, v]) => [k, v.codonsAnalysed]));
  check("only Gag and Pol codons analysed, Gag only via the gag-pol overlap",
        Object.keys(analysed).sort().join(",") === "Gag,Pol",
        JSON.stringify(analysed));
  check("every analysed Gag codon lies inside the gag-pol overlap",
        out.aa.filter((r) => r.protein === "Gag")
           .every((r) => r.query_nt_start >= 1 && reference.proteins.Gag.codon_nt[r.position - 1][2] <= 1503));
  check("no spline fallback inside an aligned fragment", out.summary.fallbackPositions === 0);
  console.log(`       ${out.nt.length} nt rows, ${out.aa.length} aa rows, ` +
              `${out.summary.codonsAnalysed} codons`);
}

// ---------------------------------------------------------------------------------
section("3. a 9 bp insertion relative to the reference");
{
  const start = 2000, end = 3000, insertAt = 500;   // insertion 500 bp into the fragment
  const base = REF.slice(start - 1, end);
  const frag = base.slice(0, insertAt) + "ACGTACGTA" + base.slice(insertAt);
  const a = alignToReference(frag, REF);
  const inserted = a.unalignedSpans.reduce((n, s) => n + (s.queryEnd - s.queryStart + 1), 0);
  check("the inserted bases are reported as unaligned", inserted === 9,
        `${inserted} unaligned bases in ${a.unalignedSpans.length} span(s)`);
  // The decisive assertion: everything downstream of the insertion still maps in frame.
  // Query base q (1-based) past the insertion corresponds to fragment base q-9, i.e.
  // reference position start + (q - 9) - 1 = start + q - 10.
  let downstreamOk = true, firstBad = null;
  for (let q = insertAt + 10; q <= frag.length; q++) {
    const expected = start + q - 10;
    if (a.queryToRef[q - 1] !== expected) { downstreamOk = false; firstBad = q; break; }
  }
  check("reference positions after the insertion are unshifted", downstreamOk,
        firstBad ? `query base ${firstBad} maps to ${a.queryToRef[firstBad - 1]}` : "");
  check("band was widened or did not need to be", a.band >= 64, `band=${a.band}`);
}

// ---------------------------------------------------------------------------------
section("4. a nef amplicon straddling the end of the modelled region");
{
  const start = 8100, end = REF.length;             // 8100-8618, crosses 8213
  const frag = REF.slice(start - 1, end);
  const a = alignToReference(frag, REF);
  const out = engine.customAtlas({ name: "nef amplicon", sequence: frag, alignment: a });
  check("positions beyond the modelled region are dropped and counted",
        out.summary.droppedOutsideRegion === REF.length - LAST,
        `dropped ${out.summary.droppedOutsideRegion}, expected ${REF.length - LAST}`);
  check("no nucleotide row past the last modelled position",
        out.nt.every((r) => r.ref_position <= LAST));
  const nef = out.summary.proteins.Nef;
  check("nef codons are analysed", nef.codonsAnalysed > 0, JSON.stringify(nef));
  check("nef analysis stops at the covered range",
        nef.lastAnalysed <= 72, `lastAnalysed=${nef.lastAnalysed}`);
  console.log(`       nef codons ${nef.firstAnalysed}-${nef.lastAnalysed} of ` +
              `${nef.codonsInReference} modelled; ${out.summary.droppedOutsideRegion} positions dropped`);
}

// ---------------------------------------------------------------------------------
section("5. a sequence lying entirely outside the modelled region");
{
  const frag = REF.slice(8300 - 1);                 // 8300-8618, all past 8213
  const a = alignToReference(frag, REF);
  const out = engine.customAtlas({ name: "past the boundary", sequence: frag, alignment: a });
  check("it aligns (the reference is not trimmed before alignment)", a.refStart === 8300,
        `refStart=${a.refStart}`);
  check("but yields no rows at all", out.nt.length === 0 && out.aa.length === 0,
        `${out.nt.length} nt, ${out.aa.length} aa`);
  check("and the drop is counted so the UI can explain itself",
        out.summary.droppedOutsideRegion === frag.length,
        `${out.summary.droppedOutsideRegion}/${frag.length}`);
}

// ---------------------------------------------------------------------------------
section("6. a reverse-complemented paste");
{
  const frag = reverseComplement(REF.slice(1311 - 1, 2600));
  const a = alignToReference(frag, REF);
  check("orientation is detected", a.orientation === "reverse_complement", a.orientation);
  check("and it aligns as well as the forward strand does", a.identity === 1,
        `${a.identity}`);
}

// ---------------------------------------------------------------------------------
section("7. junk input");
{
  const a = alignToReference("ACGT".repeat(50), REF, { tryReverseComplement: true });
  console.log(`       unalignable=${!!a.unalignable}, identity=${a.identity.toFixed(3)}, ` +
              `aligned=${a.alignedCount}`);
  check("a non-HIV sequence is either refused or aligns poorly",
        a.unalignable || a.identity < 0.8, `identity ${a.identity}`);
}

// ---------------------------------------------------------------------------------
section("8. multi-FASTA parsing");
{
  const text = ">one\nACGTACGT\nACGT\n>two desc here\nTTTT\n";
  const recs = parseFasta(text);
  check("two records", recs.length === 2);
  check("lines are joined", recs[0].sequence === "ACGTACGTACGT", recs[0].sequence);
  check("header kept", recs[1].name === "two desc here", recs[1].name);
  check("a bare sequence with no header is accepted",
        parseFasta("ACGTACGT")[0].sequence === "ACGTACGT");
}

// ---------------------------------------------------------------------------------
section("9. a full LTR-to-LTR genome (query LONGER than the reference)");
{
  // 5' LTR + leader before gag, and 3' LTR after nef, neither in the gag-nef reference.
  const rnd = (n, seed) => { let s = seed >>> 0, o = ""; const b = "ACGT";
    for (let i = 0; i < n; i++) { s = (s * 1103515245 + 12345) >>> 0; o += b[s % 4]; } return o; };
  const fivePrime = rnd(634, 7), threePrime = rnd(450, 99);
  const genome = fivePrime + REF + threePrime;
  const a = alignToReference(genome, REF);
  check("a query longer than the reference still aligns", !a.unalignable,
        a.unalignable ? a.reason : "");
  check("the whole reference is covered (refSpan 1..n)",
        a.refStart === 1 && a.refEnd === REF.length, `refSpan ${a.refStart}-${a.refEnd}`);
  // the gag-nef core must map 1:1 onto the reference, in frame
  let bad = 0;
  for (let k = 1; k <= REF.length; k++) if (a.queryToRef[634 + k - 1] !== k) bad++;
  check("the gag-nef core maps 1:1 to the reference", bad === 0, `${bad} positions off`);
  check("the two LTR overhangs are reported as unaligned, not forced in",
        a.unalignedSpans.length === 2 &&
        a.unalignedSpans[0].queryStart === 1 && a.unalignedSpans[0].queryEnd === 634 &&
        a.unalignedSpans[1].queryStart === 635 + REF.length,
        JSON.stringify(a.unalignedSpans));
  const out = engine.customAtlas({ name: "full genome", sequence: genome, alignment: a });
  check("the engine returns the full covered atlas",
        out.summary.codonsAnalysed === 2984 && out.nt.length === 8213 * 3,
        `${out.summary.codonsAnalysed} codons, ${out.nt.length} nt rows`);
  check("positions past 8,213 are dropped and counted",
        out.summary.droppedOutsideRegion === REF.length - LAST,
        `dropped ${out.summary.droppedOutsideRegion}`);
  console.log(`       genome ${genome.length} bp -> aligned ${a.alignedCount}, ` +
              `${out.summary.codonsAnalysed} codons, ${out.summary.droppedOutsideRegion} dropped past 8,213`);
}

// --- RNA spelling -----------------------------------------------------------------------------
// A pasted RNA sequence used to align at ~78% identity and silently lose a third of its codons,
// because U is not one of the four bases the model is defined on. It must now behave exactly
// like the DNA spelling of the same sequence.
{
  const frag = reference.sequence.slice(2000, 3000);
  const dna = parseFasta(">x\n" + frag)[0].sequence;
  const rna = parseFasta(">x\n" + frag.replace(/T/g, "U"))[0].sequence;
  check("U is mapped to T on input", rna === dna);
  check("no U survives parsing", !/U/.test(rna));
  const aDna = alignToReference(dna, reference.sequence);
  const aRna = alignToReference(rna, reference.sequence);
  check("an RNA sequence aligns identically to its DNA spelling",
        Math.abs(aDna.identity - aRna.identity) < 1e-12 && aRna.identity > 0.999,
        `DNA ${(100 * aDna.identity).toFixed(1)}% vs RNA ${(100 * aRna.identity).toFixed(1)}%`);
  // Lower case and alignment gap characters were already handled; pin them so they stay that way.
  check("lower case and gap characters are normalised away",
        parseFasta(">x\n" + frag.slice(0, 60).toLowerCase())[0].sequence === frag.slice(0, 60) &&
        parseFasta(">x\n" + frag.slice(0, 30) + "---" + frag.slice(30, 60))[0].sequence
          === frag.slice(0, 60));
}

console.log("\n" + "=".repeat(60));
console.log(failures === 0 ? "ALL ALIGNMENT CASES PASSED" : `${failures} ALIGNMENT CASE(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
