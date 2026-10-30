// SPDX-License-Identifier: MIT
//
// mutations.js -- parse and resolve mutation strings like "RT K65R".
//
// Kept out of the page files and free of the DOM so it can be tested directly
// (test/mutation_cases.mjs). The table-analysis view is a primary user surface and a
// parser that silently mis-resolves a subunit offset would quietly return the mutability
// of the wrong position.
//
// Accepted forms, case-insensitive on the protein name:
//     RT K65R      subunit abbreviation, reference amino acid given
//     RT 65R       reference amino acid omitted
//     Pol M184V    gene product directly
//     Env N160K    reference amino acid checked against the Env lookup
//     PR:M46I      colon separator
//     PR M46I      space separator
// The reference amino acid, when given, is CHECKED and a mismatch is reported.

const MUTATION_RE = /^\s*([A-Za-z0-9]+)\s*[:\s]\s*([A-Za-z*]?)(\d+)([A-Za-z*])\s*$/;

/** Split a line into its parts without resolving it against the atlas. */
export function parseMutationLine(line) {
  const m = MUTATION_RE.exec(line);
  if (!m) return null;
  const [, name, refAa, position, altAa] = m;
  return {
    name,
    refAa: refAa ? refAa.toUpperCase() : "",
    position: Number(position),
    altAa: altAa.toUpperCase(),
  };
}

/**
 * Resolve a mutation line to an atlas row.
 *
 * @param {string} line
 * @param {object} ctx
 *   ctx.engine    an MraEngine, for subunit offsets and HXB2 labels
 *   ctx.proteins  iterable of gene-product names present in the atlas
 *   ctx.lookup    (protein, genePosition) -> array of amino-acid rows, or undefined
 * @returns {object} either {error} or the resolved row
 */
export function resolveMutation(line, { engine, proteins, lookup }) {
  const trimmed = line.trim();
  const parsed = parseMutationLine(trimmed);
  if (!parsed) {
    return { line: trimmed, error: 'could not parse; expected something like "RT K65R"' };
  }
  const { name, refAa, position, altAa } = parsed;

  let protein = null, genePosition = null, label = name;
  const direct = [...proteins].find((p) => p.toUpperCase() === name.toUpperCase());
  if (direct) {
    protein = direct; genePosition = position; label = direct;
  } else {
    const hit = engine.resolveSubProtein(name, position);
    if (hit) { protein = hit.protein; genePosition = hit.position; label = hit.subprotein; }
  }
  if (!protein) {
    return { line: trimmed, error: `unknown protein or subunit "${name}"` };
  }

  const rows = lookup(protein, genePosition);
  if (!rows || !rows.length) {
    return {
      line: trimmed,
      error: protein === "Nef"
        ? `nef ${position} is beyond the modelled region (72 of about 206 nef codons are covered)`
        : `${label} ${position} is outside the modelled region`,
      protein, gene_position: genePosition,
    };
  }
  const hit = rows.find((r) => r.alt_aa === altAa);
  if (!hit) {
    return { line: trimmed, error: `no amino-acid change to "${altAa}" at ${label} ${position}` };
  }
  const observedRef = rows[0].ref_aa;
  return {
    line: trimmed,
    label: `${label} ${position}`,
    hxb2: engine.hxb2Label(protein, genePosition),
    hxb2_change: engine.hxb2Mutation(protein, genePosition, observedRef, altAa),
    protein,
    gene_position: genePosition,
    aa_change: hit.aa_change,
    ref_aa: observedRef,
    alt_aa: altAa,
    mu_aa: hit.mu_aa,
    mu_aa_lo95: hit.mu_aa_lo95,
    mu_aa_hi95: hit.mu_aa_hi95,
    min_nt_changes: hit.min_nt_changes,
    n_paths: hit.n_paths,
    note: refAa && refAa !== observedRef
      ? `reference amino acid given as ${refAa} but the atlas has ${observedRef} here`
      : "",
  };
}
