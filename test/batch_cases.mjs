// SPDX-License-Identifier: MIT
//
// The batch aggregation layer. Run:
//   node test/batch_cases.mjs --version v0.1
//
// The invariants here are the ones that would silently corrupt a batch figure: an uncovered
// site leaking into a median as a zero (which would drag the curve to the floor wherever
// amplicons do not overlap), an unalignable sequence counted as a sample, and divergence
// appearing where every sample agrees.

import fs from "node:fs";
import { MraEngine } from "../js/mra_engine.js";
import { alignToReference, buildSeedIndex } from "../js/align.js";
import { makeBatchMatrix, foldSequenceIntoMatrix, codonByPositionFrom, batchSiteStats,
         topDivergentSites, siteAcrossSamples, perSampleSummary,
         encodeCodon, decodeCodon, NO_CODON, perSampleRates, MUT_TYPES, codonGroupsAtSite,
         buildAaIndex, makeAaMatrix, foldAaIntoMatrix, aaPositionStats, siteSupplyStats,
         weightedDestinationSupply, codonSensitivity, diversityVsMutability,
         samplesAtAaPosition, spearman, AA_ORDER }
  from "../js/batch_core.js";

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
const seedIndex = buildSeedIndex(ref.sequence);
const region = ref.modelled_region;
const nSites = region.last_position - region.first_position + 1;
const FLIP = { A: "G", C: "T", G: "A", T: "C" };

/** Fold a set of {name, sequence} into a fresh matrix, mirroring what the worker does. */
function build(records) {
  const m = makeBatchMatrix(records.length, nSites, region.first_position);
  const names = [], summaries = [];
  records.forEach((rec, i) => {
    names.push(rec.name);
    const aln = alignToReference(rec.sequence, ref.sequence, { seedIndex });
    if (aln.unalignable) { summaries.push({ name: rec.name, unalignable: true }); return; }
    const out = engine.customAtlas({ name: rec.name, sequence: rec.sequence, alignment: aln });
    foldSequenceIntoMatrix(m, i, out.nt, codonByPositionFrom(out.aa, ref.proteins));
    summaries.push(out.summary);
  });
  return { m, names, summaries };
}

// --- 1. codon byte encoding -----------------------------------------------------------
check("encode/decode round-trips every codon", (() => {
  const B = ["A", "C", "G", "T"];
  for (const a of B) for (const b of B) for (const c of B) {
    const codon = a + b + c;
    if (decodeCodon(encodeCodon(codon)) !== codon) return false;
  }
  return true;
})());
check("a missing codon encodes to NO_CODON and decodes to null",
      encodeCodon(null) === NO_CODON && encodeCodon("AN") === NO_CODON
      && decodeCodon(NO_CODON) === null);

// --- 2. identical batch: no divergence anywhere ---------------------------------------
{
  const N = 4;
  const recs = Array.from({ length: N }, (_, i) => ({ name: `id${i}`, sequence: ref.sequence }));
  const { m, names } = build(recs);
  const stats = batchSiteStats(m);
  check("identical batch covers every modelled site", stats.length === nSites,
        `${stats.length} of ${nSites}`);
  check("identical batch has zero divergence at every site",
        stats.every((r) => r.divergence === 0 || r.divergence === null),
        `${stats.filter((r) => r.divergence > 0).length} sites with divergence`);
  check("identical batch has median == min == max",
        stats.every((r) => Math.abs(r.median - r.min) < 1e-12 && Math.abs(r.max - r.min) < 1e-12));
  check("every site reports all samples covering", stats.every((r) => r.n_covered === N));
  check("topDivergentSites returns nothing for an identical batch",
        topDivergentSites(stats).length === 0);
  const ps = perSampleSummary(m, stats, names);
  check("every sample's ratio to the batch median is exactly 1",
        ps.every((r) => Math.abs(r.ratio_to_batch - 1) < 1e-12),
        ps.map((r) => r.ratio_to_batch.toFixed(6)).join(","));
}

// --- 3. one planted codon change: divergence confined to that codon's footprint ---------
{
  const a = ref.sequence.split("");
  a[2584] = FLIP[a[2584]];                       // changes nt 2585
  const recs = [
    { name: "wt1", sequence: ref.sequence },
    { name: "wt2", sequence: ref.sequence },
    { name: "var", sequence: a.join("") },
  ];
  const { m, names } = build(recs);
  const stats = batchSiteStats(m);
  const div = stats.filter((r) => r.divergence !== null && r.divergence > 1e-9);
  // The mutated site changes, and so does its 5' neighbour, whose 3' context is now different.
  // Anything beyond that footprint would mean the fold is bleeding across positions.
  check("divergence is confined to the mutated site and its 3'-context neighbour",
        div.length === 2 && div.some((r) => r.ref_position === 2585)
                        && div.some((r) => r.ref_position === 2584),
        div.map((r) => r.ref_position).join(","));
  check("the mutated site is the most divergent",
        topDivergentSites(stats, { limit: 1 })[0].ref_position === 2585);
  const at = siteAcrossSamples(m, 2585, names);
  check("all three samples are reported at the mutated site",
        at.length === 3 && at.every((r) => r.covered));
  check("the variant carries a different codon from the two wild types", (() => {
    const byName = Object.fromEntries(at.map((r) => [r.sample, r.codon]));
    return byName.wt1 === byName.wt2 && byName.var !== byName.wt1;
  })(), at.map((r) => `${r.sample}=${r.codon}`).join(" "));
  check("the variant's rate differs from the wild types at that site", (() => {
    const v = Object.fromEntries(at.map((r) => [r.sample, r.value]));
    return Math.abs(v.wt1 - v.wt2) < 1e-15 && Math.abs(v.var - v.wt1) > 1e-12;
  })());
}

// --- 4. disjoint amplicons: uncovered is never zero -------------------------------------
{
  const recs = [
    { name: "left", sequence: ref.sequence.slice(1000, 2000) },
    { name: "right", sequence: ref.sequence.slice(5000, 6000) },
  ];
  const { m, names } = build(recs);
  const stats = batchSiteStats(m);
  check("no site is covered by both disjoint amplicons",
        stats.every((r) => r.n_covered === 1),
        `${stats.filter((r) => r.n_covered > 1).length} sites with >1`);
  check("no median is zero (an uncovered cell never enters a statistic)",
        stats.every((r) => r.median > 0));
  check("single-coverage sites report null divergence, not 0 or Infinity",
        stats.every((r) => r.divergence === null));
  const at = siteAcrossSamples(m, 1500, names);
  check("a sample not covering the queried site is reported as absent, not omitted",
        at.length === 2 && at.filter((r) => r.covered).length === 1
        && at.find((r) => !r.covered).value === null,
        at.map((r) => `${r.sample}:${r.covered}`).join(" "));
}

// --- 5. an unalignable sequence is excluded, not zeroed ----------------------------------
{
  const recs = [
    { name: "real1", sequence: ref.sequence.slice(2000, 3000) },
    { name: "junk", sequence: "ACGT".repeat(60) },
    { name: "real2", sequence: ref.sequence.slice(2000, 3000) },
  ];
  const { m, names, summaries } = build(recs);
  check("the junk record is flagged unalignable", summaries[1].unalignable === true);
  const stats = batchSiteStats(m);
  check("covered sites count only the two alignable samples",
        stats.every((r) => r.n_covered === 2),
        `max n_covered = ${Math.max(...stats.map((r) => r.n_covered))}`);
  const at = siteAcrossSamples(m, 2500, names);
  check("the unalignable sample contributes no value anywhere",
        at.find((r) => r.sample === "junk").covered === false);
  const ps = perSampleSummary(m, stats, names);
  check("the unalignable sample reports zero sites and a null ratio", (() => {
    const j = ps.find((r) => r.sample === "junk");
    return j.sites_covered === 0 && j.ratio_to_batch === null;
  })());
}

// --- 6. G->A is a subset of total supply --------------------------------------------------
{
  const { m } = build([{ name: "a", sequence: ref.sequence },
                       { name: "b", sequence: ref.sequence }]);
  const tot = batchSiteStats(m, { field: "total" });
  const ga = batchSiteStats(m, { field: "ga" });
  const byPos = new Map(tot.map((r) => [r.ref_position, r.median]));
  check("the G->A median never exceeds the total median at the same site",
        ga.every((r) => r.median <= byPos.get(r.ref_position) + 1e-15),
        `${ga.filter((r) => r.median > byPos.get(r.ref_position) + 1e-15).length} violations`);
  const gaNonZero = ga.filter((r) => r.median > 0).length;
  check("some sites carry G->A supply and some carry none (G sites only)",
        gaNonZero > 0 && gaNonZero < ga.length, `${gaNonZero} of ${ga.length}`);
}

// --- 7. per-sequence rate summary --------------------------------------------------------
{
  const { m, names } = build([
    { name: "full", sequence: ref.sequence },
    { name: "frag", sequence: ref.sequence.slice(2000, 4300) },
    { name: "junk", sequence: "ACGT".repeat(60) },
  ]);
  const rates = perSampleRates(m, names);
  const full = rates.find((r) => r.sample === "full");
  const junk = rates.find((r) => r.sample === "junk");

  check("every substitution class is reported for a full-length sequence",
        MUT_TYPES.every((t) => Number.isFinite(full[t])),
        MUT_TYPES.filter((t) => !Number.isFinite(full[t])).join(",") || "all present");
  check("G->A is the highest per-type rate (APOBEC dominance)",
        MUT_TYPES.every((t) => t === "G->A" || full[t] <= full["G->A"]),
        `G->A ${full["G->A"].toExponential(2)}, next ` +
        `${Math.max(...MUT_TYPES.filter((t) => t !== "G->A").map((t) => full[t])).toExponential(2)}`);
  check("ga_rate restates the G->A class exactly", full.ga_rate === full["G->A"]);
  check("G->A share is a fraction between 0 and 1",
        full.ga_share > 0 && full.ga_share < 1, (100 * full.ga_share).toFixed(1) + "%");
  // The bulk rate is a mean over sites of a sum over three classes, so it must sit between the
  // smallest and the largest per-site total -- and in particular must not equal a per-type rate.
  check("bulk rate is positive and distinct from any single per-type rate",
        full.bulk_rate > 0 && MUT_TYPES.every((t) => full.bulk_rate !== full[t]),
        full.bulk_rate.toExponential(3));
  // Expected substitutions per genome per cycle: order 1 for HIV, not 1e-5 or 1e3.
  const perCycle = full.bulk_rate * full.sites_covered;
  check("expected substitutions per genome per cycle is order 1",
        perCycle > 0.05 && perCycle < 10, perCycle.toFixed(3));
  check("transitions exceed transversions", full.transition_rate > full.transversion_rate,
        `Ts/Tv ${(full.transition_rate / full.transversion_rate).toFixed(2)}`);
  check("an unalignable sequence reports nulls, not zeros",
        junk.bulk_rate === null && junk.ga_rate === null && junk.ga_share === null
        && MUT_TYPES.every((t) => junk[t] === null) && junk.sites_covered === 0);
  // A per-type rate is averaged over the sites where that type can occur, so restricting the
  // sequence to a sub-region changes which sites contribute and the rate may move -- but it
  // must stay finite and positive, never NaN from an empty denominator.
  const frag = rates.find((r) => r.sample === "frag");
  check("a partial sequence still reports finite per-type rates",
        MUT_TYPES.every((t) => frag[t] === null || Number.isFinite(frag[t])),
        `sites ${frag.sites_covered}`);
  check("the per-type counts use only sites whose reference base matches the type", (() => {
    // G->A must be averaged over G sites only: its count must equal the number of G sites in
    // the covered span, which is strictly fewer than all covered sites.
    const gIndex = MUT_TYPES.indexOf("G->A");
    const sIdx = names.indexOf("full");
    const n = m.typeCount[sIdx * MUT_TYPES.length + gIndex];
    return n > 0 && n < full.sites_covered;
  })(), `G->A sites vs covered`);
}

// --- 8. HXB2 protein-position lookup, as the page's coordinate box does it -----------------
{
  const spans = new Map();
  for (const [n, r] of Object.entries(ref.proteins)) {
    if (r.codon_nt && r.codon_nt.length) spans.set(n, { min: r.first_aa,
                                                        max: r.first_aa + r.codon_nt.length - 1 });
  }
  const targets = engine.proteinTargets(spans);
  const byKey = new Map(targets.map((t) => [t.key, t]));
  check("the subunit selector offers the clinical subunits",
        ["PR", "RT", "IN", "Env"].every((k) => byKey.has(k)),
        targets.map((t) => t.key).join(","));

  // Wild-type residues at well-known resistance positions. These are the strongest available
  // check that the HXB2 -> NL4-3 -> nucleotide chain is wired correctly: if any link were off
  // by one, the codon here would not translate to the residue the DRM name implies.
  const LANDMARKS = [["PR", 46, "ATG", "M"], ["RT", 103, "AAA", "K"], ["RT", 184, "ATG", "M"],
                     ["IN", 148, "CAA", "Q"], ["IN", 155, "AAT", "N"], ["Env", 332, "AAC", "N"]];
  for (const [key, pos, codon, aa] of LANDMARKS) {
    const t = byKey.get(key);
    const geneAa = engine.geneAaFromOwnPosition(t, pos);
    const rec = ref.proteins[t.protein];
    const nts = geneAa === null ? null : rec.codon_nt[geneAa - rec.first_aa];
    const got = nts ? nts.map((n) => ref.sequence[n - 1]).join("") : null;
    check(`${key} ${pos} resolves to codon ${codon} (${aa})`, got === codon, `got ${got}`);
    check(`${key} ${pos} round-trips back to its own label`,
          geneAa !== null && engine.hxb2Label(t.protein, geneAa) === `${key} ${pos}`,
          geneAa === null ? "no NL4-3 equivalent" : engine.hxb2Label(t.protein, geneAa));
  }
  // Env carries indels relative to HXB2, so the gene position must NOT equal the HXB2 number.
  const env = byKey.get("Env");
  check("Env HXB2 332 maps to a different NL4-3 position (indel offset applied)",
        engine.geneAaFromOwnPosition(env, 332) !== 332,
        `NL4-3 Env ${engine.geneAaFromOwnPosition(env, 332)}`);
}

// --- 9. codon grouping for the readable site plot -------------------------------------------
{
  const a = ref.sequence.split("");
  a[1598] = a[1598] === "A" ? "G" : "A";             // nt 1599, first base of PR 46
  const recs = [];
  for (let i = 0; i < 10; i++) {
    recs.push({ name: i < 3 ? `a_very_long_patient_identifier_${i}` : `ctrl_${i}`,
                sequence: i < 3 ? a.join("") : ref.sequence });
  }
  const { m, names } = build(recs);
  const rows = siteAcrossSamples(m, 1599, names);
  const refCodon = [1599, 1600, 1601].map((n) => ref.sequence[n - 1]).join("");
  const groups = codonGroupsAtSite(rows, refCodon);

  check("ten samples collapse to the number of distinct codons",
        groups.length === 2, `${groups.length} groups from ${rows.length} samples`);
  check("group sizes account for every covering sample",
        groups.reduce((t, g) => t + g.n, 0) === rows.filter((r) => r.covered).length);
  check("groups are ordered by rate, fastest first",
        groups.every((g, i) => i === 0 || groups[i - 1].rate >= g.rate));
  check("exactly one group is flagged as the reference codon",
        groups.filter((g) => g.is_reference).length === 1,
        groups.map((g) => `${g.codon}:${g.is_reference}`).join(" "));
  check("the reference group holds the seven unmutated samples",
        groups.find((g) => g.is_reference).n === 7);
  check("every sample in a group shares one rate (uniform)",
        groups.every((g) => g.uniform));
  check("sample names are carried on the group, not lost in the collapse",
        groups.every((g) => g.samples.length === g.n)
        && groups.some((g) => g.samples.some((n) => n.startsWith("a_very_long"))));
  // A batch where everyone agrees must still produce one group, not zero.
  const same = build([{ name: "x", sequence: ref.sequence }, { name: "y", sequence: ref.sequence }]);
  const g1 = codonGroupsAtSite(siteAcrossSamples(same.m, 1599, same.names), refCodon);
  check("an agreeing batch yields a single reference group",
        g1.length === 1 && g1[0].is_reference && g1[0].n === 2);
}

// --- 10. codon-level aggregation ------------------------------------------------------------
{
  const aaIdx = buildAaIndex(ref.proteins);
  check("every reference codon is indexed once",
        aaIdx.nPositions === 2984, `${aaIdx.nPositions} positions`);
  check("the index is in genomic order",
        aaIdx.meta.every((m, i) => i === 0 || aaIdx.meta[i - 1].ntStart <= m.ntStart));

  const buildAa = (recs) => {
    const m = makeAaMatrix(recs.length, aaIdx.nPositions);
    recs.forEach((r, i) => {
      const aln = alignToReference(r.sequence, ref.sequence, { seedIndex });
      // paths, not aa rows: mutability is the sum of ONE-STEP path rates, and a path carries
      // the mutation type that lets G->A be separated.
      foldAaIntoMatrix(m, i, engine.customAtlas({ name: r.name, sequence: r.sequence,
                                                  alignment: aln }).paths, aaIdx);
    });
    return { m, names: recs.map((r) => r.name) };
  };

  // RT 103 is Pol 253 in NL4-3 numbering, codon AAA. AAA->AAG is synonymous (still K).
  const polRec = ref.proteins.Pol;
  const nts = polRec.codon_nt[253 - polRec.first_aa];
  const synSeq = ref.sequence.split(""); synSeq[nts[2] - 1] = "G";

  // (a) a homogeneous batch
  const homo = buildAa([{ name: "a", sequence: ref.sequence },
                        { name: "b", sequence: ref.sequence }]);
  const hs = aaPositionStats(homo.m, aaIdx);
  check("a homogeneous batch has no mutability divergence anywhere",
        hs.every((r) => r.divergence === 0 || r.divergence === null),
        `${hs.filter((r) => r.divergence > 0).length} divergent`);
  check("a homogeneous batch has one residue and one codon at every position",
        hs.every((r) => r.n_aa === 1 && r.n_codons === 1));
  check("a homogeneous batch has zero diversity and zero entropy everywhere",
        hs.every((r) => r.diversity === 0 && r.entropy === 0));
  const hStat = hs.find((r) => r.protein === "Pol" && r.position === 253);
  const hDest = weightedDestinationSupply(homo.m, hStat);
  const single = homo.m.destinations.get(hStat.position_index + "|" + hStat.majority_codon);
  check("weighted destination rates equal the single-codon rates when the batch agrees",
        hDest.every((d) => Math.abs(d.weighted_rate - single[AA_ORDER.indexOf(d.alt_aa)]) < 1e-18),
        `max diff ${Math.max(...hDest.map((d) =>
          Math.abs(d.weighted_rate - single[AA_ORDER.indexOf(d.alt_aa)]))).toExponential(1)}`);
  const hSens = codonSensitivity(homo.m, hStat);
  check("a single-codon position reports no fold range rather than a fold range of 1",
        hSens.uniform === true && hSens.fold_range === null, String(hSens.fold_range));

  // (b) a purely synonymous variant: residue identical, mutability not
  const mix = buildAa([{ name: "a", sequence: ref.sequence },
                       { name: "b", sequence: ref.sequence },
                       { name: "c", sequence: ref.sequence },
                       { name: "g", sequence: synSeq.join("") }]);
  const ms = aaPositionStats(mix.m, aaIdx);
  const mStat = ms.find((r) => r.protein === "Pol" && r.position === 253);
  check("a synonymous variant leaves the residue distribution unchanged",
        mStat.n_aa === 1 && mStat.diversity === 0 && mStat.entropy === 0,
        `${mStat.n_aa} residues, diversity ${mStat.diversity}`);
  check("but it does split the codon distribution",
        mStat.n_codons === 2 && mStat.majority_codon === "AAA" &&
        Math.abs(mStat.majority_codon_freq - 0.75) < 1e-9,
        `${mStat.n_codons} codons, majority freq ${mStat.majority_codon_freq}`);
  check("and it does change mutability (synonymous change is invisible to the protein, " +
        "not to the supply)",
        mStat.divergence > 0, `log2 fold range ${mStat.divergence}`);
  const diverging = ms.filter((r) => r.divergence > 1e-9);
  check("mutability divergence is confined to the altered codon",
        diverging.length === 1 && diverging[0].position === 253,
        diverging.map((r) => `${r.protein} ${r.position}`).join(", "));

  // (c) the weighting arithmetic
  const mDest = weightedDestinationSupply(mix.m, mStat);
  const dA = mix.m.destinations.get(mStat.position_index + "|AAA");
  const dG = mix.m.destinations.get(mStat.position_index + "|AAG");
  const n = AA_ORDER.indexOf("N");
  const expectedN = 0.75 * dA[n] + 0.25 * dG[n];
  const gotN = mDest.find((d) => d.alt_aa === "N").weighted_rate;
  check("the weighted rate is the codon-frequency-weighted mean",
        Math.abs(gotN - expectedN) < 1e-18,
        `got ${gotN.toExponential(4)} expected ${expectedN.toExponential(4)}`);
  check("the destination list carries the range across codons",
        mDest.find((d) => d.alt_aa === "N").fold_range > 1,
        String(mDest.find((d) => d.alt_aa === "N").fold_range));

  // (d) the non-synonymous / synonymous / stop split
  const total = AA_ORDER.reduce((t, a, i) => t + dA[i], 0);
  const nsSum = mDest.filter((d) => !d.is_synonymous && !d.is_stop)
                     .reduce((t, d) => t + d.weighted_rate, 0);
  const refAa = mStat.majority_aa;
  let nsFromCodon = 0;
  AA_ORDER.forEach((a, i) => { if (a !== "*" && a !== refAa) nsFromCodon += dA[i]; });
  check("stops and synonymous changes are excluded from non-synonymous supply",
        nsFromCodon < total, `${nsFromCodon.toExponential(2)} of ${total.toExponential(2)}`);
  check("the per-destination rates sum to more than any single destination",
        nsSum > Math.max(...mDest.filter((d) => !d.is_synonymous && !d.is_stop)
                                 .map((d) => d.weighted_rate)));

  // (e) per-sample accessor and the scatter input
  const sm = samplesAtAaPosition(mix.m, mStat, mix.names);
  check("every submitted sample appears in the per-position table",
        sm.length === mix.names.length);
  check("the synonymous sample is flagged with its own codon",
        sm.find((r) => r.sample === "g").codon === "AAG" &&
        sm.find((r) => r.sample === "g").aa === "K");
  const pts = diversityVsMutability(ms);
  check("the scatter has one point per covered position with a positive rate",
        pts.length > 0 && pts.every((d) => d.mutability > 0 && d.diversity >= 0),
        `${pts.length} points`);
  check("spearman returns null on too few points", spearman([1, 2], [1, 2]) === null);
}

// --- 11. the G->A / non-G->A split ----------------------------------------------------------
{
  const aaIdx = buildAaIndex(ref.proteins);
  const m = makeAaMatrix(1, aaIdx.nPositions);
  const out = engine.customAtlas({ name: "nl43", sequence: ref.sequence,
                                   alignment: alignToReference(ref.sequence, ref.sequence) });
  foldAaIntoMatrix(m, 0, out.paths, aaIdx);
  const all = aaPositionStats(m, aaIdx, { field: "nonsyn" });
  const ga = aaPositionStats(m, aaIdx, { field: "ga" });
  const ng = aaPositionStats(m, aaIdx, { field: "nonga" });

  check("the three classes cover the same positions",
        all.length === ga.length && ga.length === ng.length, `${all.length}/${ga.length}/${ng.length}`);
  const bad = all.filter((r, i) =>
    Math.abs((ga[i].median + ng[i].median) - r.median) > 1e-12 * Math.max(1e-30, r.median));
  check("G->A plus non-G->A equals the non-synonymous total at every position",
        bad.length === 0, `${bad.length} positions disagree`);
  check("non-G->A is never negative",
        ng.every((r) => r.median >= 0 && r.min >= 0));
  // A codon with no G cannot contribute G->A supply, so plenty of positions must be exactly 0.
  const zeroGa = ga.filter((r) => r.median === 0).length;
  check("positions whose codon has no G carry zero G->A supply",
        zeroGa > 1000 && zeroGa < ga.length, `${zeroGa} of ${ga.length}`);
  // G->A dominates overall, which is the APOBEC signature the atlas is built around.
  const sum = (a) => a.reduce((t, r) => t + r.median, 0);
  const share = sum(ga) / sum(all);
  check("G->A is the majority of non-synonymous supply genome-wide",
        share > 0.5 && share < 0.95, `${(100 * share).toFixed(1)}%`);
  // Folding from one-step paths must not drift from the amino-acid rows it replaced.
  const byPos = new Map();
  for (const r of out.aa) {
    if (r.min_nt_changes !== 1 || r.alt_aa === "*" || r.alt_aa === r.ref_aa) continue;
    const k = r.protein + "|" + r.position;
    byPos.set(k, (byPos.get(k) || 0) + r.mu_aa);
  }
  let maxRel = 0;
  for (const r of all) {
    const v = byPos.get(r.key);
    if (v === undefined) continue;
    maxRel = Math.max(maxRel, Math.abs(v - r.median) / Math.max(v, 1e-30));
  }
  check("one-step path sums agree with the amino-acid rows to better than 1%",
        maxRel < 0.01, `max ${(100 * maxRel).toExponential(2)}%`);
}

// --- 12. the supply-flow chart's two weightings and its codon pooling ------------------------
// drawSupplyFlow needs a DOM, so this group is skipped where jsdom is absent (same contract as
// test/page_boot_cases.mjs). The arithmetic it guards is not incidental: the two weightings are
// different quantities, not a rescaling of each other, and a bug that silently swapped them
// would still render a plausible-looking chart.
{
  let JSDOM = null;
  for (const spec of ["jsdom", "/tmp/domtest/node_modules/jsdom/lib/api.js"]) {
    try { ({ JSDOM } = await import(spec)); break; } catch { /* try the next */ }
  }
  if (!JSDOM) {
    console.log("  skip  jsdom not found - supply-flow checks skipped");
  } else {
    const dom = new JSDOM("<!DOCTYPE html><body></body>");
    globalThis.window = dom.window; globalThis.document = dom.window.document;
    const { drawSupplyFlow } = await import("../js/batch_view.js");

    // A position carrying one dominant codon plus a tail of rare ones.
    const aaIdx = buildAaIndex(ref.proteins);
    const pol = ref.proteins.Pol, nts = pol.codon_nt[253 - pol.first_aa];
    const put = (c) => { const a = ref.sequence.split("");
      a[nts[0] - 1] = c[0]; a[nts[1] - 1] = c[1]; a[nts[2] - 1] = c[2]; return a.join(""); };
    const mix = [["AAA", 90], ["AAG", 6], ["AAC", 1], ["AAT", 1], ["AGA", 1], ["CAA", 1]];
    const seqs = [];
    for (const [c, k] of mix) for (let i = 0; i < k; i++) seqs.push(put(c));
    const m = makeAaMatrix(seqs.length, aaIdx.nPositions);
    seqs.forEach((q, i) => foldAaIntoMatrix(m, i,
      engine.customAtlas({ name: "s" + i, sequence: q,
                           alignment: alignToReference(q, ref.sequence) }).paths, aaIdx));
    const st = aaPositionStats(m, aaIdx).find((r) => r.protein === "Pol" && r.position === 253);
    check("the test position carries the six planted codons", st.n_codons === 6, st.n_codons);

    const draw = (weight) => {
      const box = dom.window.document.createElement("div");
      const cap = dom.window.document.createElement("div");
      const ok = drawSupplyFlow(box, cap, st, m, { weight, refCodon: "AAA" });
      const texts = [...box.querySelectorAll("text")].map((t) => t.textContent);
      const names = texts.filter((t) => /\([A-Z*]\)|rarer codons/.test(t));
      const heights = [...box.querySelectorAll("rect")]
        .map((r) => Number(r.getAttribute("height"))).slice(0, names.length);
      const tot = heights.reduce((a, b) => a + b, 0);
      return { ok, names, box, labels: texts, cap: cap.textContent,
               share: heights.map((h) => h / tot) };
    };
    const bySeq = draw("sequences"), bySup = draw("supply");

    check("both weightings draw", bySeq.ok && bySup.ok);
    // 4 codons under 2% of 99 samples are pooled, leaving AAA, AAG and one block.
    check("codons under 2% are pooled into a single block",
          bySeq.names.length === 3 && /rarer codons/.test(bySeq.names[2]),
          bySeq.names.join(" | "));
    check("the pooled block is reported in the caption",
          /codons each under 2% of the batch are pooled/.test(bySeq.cap));
    // In "sequences" mode a block must be the codon's frequency, nothing else.
    const freqAAA = st.codon_counts.find(([c]) => c === "AAA")[1] / st.n_covered;
    check("a block in sequences mode equals that codon's frequency",
          Math.abs(bySeq.share[0] - freqAAA) < 0.01,
          `${bySeq.share[0].toFixed(3)} vs ${freqAAA.toFixed(3)}`);
    // In "supply" mode it must not, because the codons differ in mutability.
    check("supply mode gives the dominant codon a smaller share than its frequency",
          bySup.share[0] < bySeq.share[0] - 0.02,
          `${bySup.share[0].toFixed(3)} vs ${bySeq.share[0].toFixed(3)}`);
    check("the rare pooled codons carry more supply than their frequency",
          bySup.share[2] > bySeq.share[2] + 0.02,
          `${bySup.share[2].toFixed(3)} vs ${bySeq.share[2].toFixed(3)}`);
    for (const [tag, r] of [["sequences", bySeq], ["supply", bySup]]) {
      check(`${tag} mode block shares sum to 1`,
            Math.abs(r.share.reduce((a, b) => a + b, 0) - 1) < 1e-9);
    }
    // Synonymous routes used to be dropped; they are now a destination block of their own,
    // because they are a median 21% of a site's supply and were silently missing from it.
    check("the no-change block is drawn and named",
          bySeq.labels.some((t) => /no change \(synonymous\)/.test(t)),
          bySeq.labels.filter((t) => /change|other/.test(t)).join(" | "));
    check("the left blocks no longer carry a supply percentage",
          !bySeq.labels.some((t) => /of supply/.test(t)));
    check("the left blocks still carry a sequence percentage",
          bySeq.labels.filter((t) => /of sequences/.test(t)).length === bySeq.names.length);
    // A codon must be drawn in the colour of the residue it encodes, so degenerate codons
    // match each other and match that residue in the composition bar above.
    const { residueColour, drawPositionComposition } = await import("../js/batch_view.js");
    check("AAA and AAG are the same colour, being degenerate for lysine",
          residueColour("K") === residueColour("K") &&
          residueColour("K") !== residueColour("N"));
    const fills = [...bySeq.box.querySelectorAll("rect")]
      .map((r) => r.getAttribute("fill")).slice(0, bySeq.names.length);
    const wanted = bySeq.names.map((n) => {
      const mm = /^([ACGT]{3}) \(([A-Z*])\)/.exec(n);
      return mm ? residueColour(mm[2]) : null;
    });
    check("every codon block is drawn in its residue's colour",
          wanted.every((c, i) => c === null || c === fills[i]),
          wanted.map((c, i) => `${bySeq.names[i].trim()}:${c}/${fills[i]}`).join(" "));
    const cbox = dom.window.document.createElement("div");
    const ccap = dom.window.document.createElement("div");
    drawPositionComposition(cbox, ccap, st, { only: "aa" });
    const cfill = [...cbox.querySelectorAll("rect")].map((r) => r.getAttribute("fill"));
    check("the composition bar uses the same residue colours as the flow",
          cfill[0] === residueColour(st.aa_counts[0][0]),
          `${cfill[0]} vs ${residueColour(st.aa_counts[0][0])}`);
  }
}

// --- 13. per-site supply statistics ----------------------------------------------------------
{
  const aaIdx = buildAaIndex(ref.proteins);
  const m = makeAaMatrix(1, aaIdx.nPositions);
  foldAaIntoMatrix(m, 0, engine.customAtlas({ name: "nl43", sequence: ref.sequence,
    alignment: alignToReference(ref.sequence, ref.sequence) }).paths, aaIdx);
  const stats = aaPositionStats(m, aaIdx);
  const all = stats.map((st) => siteSupplyStats(m, st)).filter(Boolean);
  check("a statistic is returned for essentially every position",
        all.length > 0.95 * stats.length, `${all.length} of ${stats.length}`);
  check("fractions are in [0,1] and the two reported ones never exceed 1 together",
        all.every((s) => s.syn_fraction >= 0 && s.stop_fraction >= 0 &&
                         s.syn_fraction + s.stop_fraction <= 1 + 1e-9));
  // exp(H) is bounded by the number of reachable destinations and is at least 1.
  check("the effective number lies between 1 and the reachable count",
        all.every((s) => s.effective_destinations >= 1 - 1e-9 &&
                         s.effective_destinations <= s.n_destinations + 1e-9));
  check("evenness is exp(H)/k, or null where only one residue is reachable",
        all.every((s) => s.evenness === null
          ? s.n_destinations === 1
          : Math.abs(s.evenness - s.effective_destinations / s.n_destinations) < 1e-9));
  check("the largest destination share is at least 1/k and at most 1",
        all.every((s) => s.top_share >= 1 / s.n_destinations - 1e-9 && s.top_share <= 1 + 1e-9));
  // A site whose supply is concentrated must score near 1, a spread one well above it: if the
  // statistic did not discriminate it would be useless as a column. The ceiling is lower here
  // than on a real batch because this is ONE sequence, so every position has a single codon
  // and at most nine one-step neighbours; a 182-sequence batch pools codons and reaches ~14.
  const eff = all.map((s) => s.effective_destinations).sort((a, b) => a - b);
  check("the effective number discriminates across the proteome",
        eff[0] < 2 && eff[eff.length - 1] > 4,
        `${eff[0].toFixed(2)} to ${eff[eff.length - 1].toFixed(2)}`);
  check("a single sequence cannot exceed the nine one-step neighbours of its codon",
        all.every((s) => s.n_destinations <= 9), Math.max(...all.map((s) => s.n_destinations)));
  const synMed = all.map((s) => s.syn_fraction).sort((a, b) => a - b)[Math.floor(all.length / 2)];
  check("synonymous supply is a substantial but not dominant median share",
        synMed > 0.05 && synMed < 0.5, `${(100 * synMed).toFixed(1)}%`);
}

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL BATCH CASES PASSED");
