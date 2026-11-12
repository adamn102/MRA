// SPDX-License-Identifier: MIT
//
// batch_view.js -- residue colours and the supply-flow chart.
//
// The mutability tracks, the site-comparison plot, the diversity scatter and the
// composition bar belonged to the Batch analysis page, which is no longer part of the
// site. What remains is what Browse and Your Sequence draw.

import { LAYOUT } from "./chart.js";
import { AA_ORDER } from "./batch_core.js";
import { formatRate } from "./table.js";

// ---------------------------------------------------------------------------------
/**
 * One colour per amino acid, used by the supply-flow chart.
 *
 * Keyed on the RESIDUE, so a codon takes the colour of what it encodes and degenerate codons
 * are necessarily the same colour -- AAA and AAG are both lysine and both read as lysine at a
 * glance. The flow chart uses this map, so a codon and the residue it encodes match.
 *
 * The assignment is fixed rather than order-of-appearance: a residue keeps its colour between
 * positions, so moving the marker does not recolour the figure under the reader. Chemically
 * similar residues are deliberately NOT grouped -- the question here is identity, not
 * property, and property-grouped schemes make the common case (telling two destinations apart)
 * harder.
 */
const AA_COLOUR = {
  A: "#0b5c8a", R: "#c2561a", N: "#2e8b9a", D: "#8a4f7d", C: "#6aaa64",
  E: "#b5485d", Q: "#3f7fa8", G: "#c2973a", H: "#7a5ea8", I: "#2f8f6f",
  L: "#a8553a", K: "#1f6fb2", M: "#9c6b2f", F: "#4a7d3a", P: "#b03a6b",
  S: "#2a9d8f", T: "#846a3c", V: "#5566b5", W: "#7d4a2a", Y: "#3d8ab5",
};
const STOP_COLOUR = "#8c1d33";      // a stop is an outcome, not a residue: kept visually apart
const NEUTRAL_COLOUR = "#9aa4ae";   // pooled blocks and the no-change block

/** The colour for a residue, falling back to neutral for stops, pools and unknowns. */
export function residueColour(aa) {
  if (aa === "*") return STOP_COLOUR;
  return AA_COLOUR[aa] || NEUTRAL_COLOUR;
}

/**
 * Where this position's residue change comes from and where it goes: codons on the left,
 * destination residues on the right, bands between them.
 *
 * TWO WEIGHTINGS, AND WHY BOTH EXIST. The chart can apportion its bands by either quantity,
 * and the difference between them is the whole argument of this page:
 *
 *   weight = "sequences"  band = freq(codon) x rate(codon->residue) / total rate of that codon
 *                         A left block is simply how common the codon is. Reading: if one
 *                         sequence drawn from this batch underwent a residue change here,
 *                         where would it land? It deliberately ignores that some codons
 *                         mutate far more readily than others.
 *
 *   weight = "supply"     band = freq(codon) x rate(codon->residue)
 *                         A left block is the codon's share of the residue-change supply.
 *                         Reading: across the whole batch, where does change actually get
 *                         generated? A codon carried by a tenth of the sequences can dominate.
 *
 * Switching between them moves the left column by exactly the per-codon mutability, so the
 * toggle is a direct readout of how much codon identity matters at this position. Both are
 * mass-conserving; neither is a rescaling of the other.
 *
 * RARE CODONS ARE POOLED. A real batch carries a long tail of codons at around one percent
 * each -- in the LANL compendium, over 20 of them at some env positions. Drawn individually
 * they are unreadable hairlines that squeeze the codons worth seeing into a thin band.
 * Anything below `minCodonFreq` (2% by default) merges into one "other codons" block. Pooling
 * happens AFTER each codon's synonymous routes are removed, because the pooled codons do not
 * all encode the same residue and "synonymous" would otherwise be ill-defined.
 *
 * SYNONYMOUS ROUTES ARE INCLUDED, as one destination block rather than scattered among the
 * residues. A synonymous change from AAA and one from GGT have nothing in common except that
 * neither alters the protein, so pooling them is semantically exact in a way that pooling two
 * different residues would not be. They are a median 21% of a site's supply across the LANL
 * compendium and range from 4% to 53%, so leaving them out understated how much of what a
 * position generates never reaches the protein.
 */
export function drawSupplyFlow(container, captionEl, stat, matrix, opts = {}) {
  container.innerHTML = "";
  const mode = opts.weight === "supply" ? "supply" : "sequences";
  const codons = stat.codon_counts.map(([c]) => c);
  const freqOf = new Map(stat.codon_counts.map(([c, n]) => [c, n / stat.n_covered]));

  // Raw per-codon rates first; the weighting is applied once, below, so the two modes cannot
  // drift apart in their handling of synonymous or unreachable routes.
  const raw = new Map();
  for (const codon of codons) {
    const d = matrix.destinations.get(stat.position_index + "|" + codon);
    if (!d) continue;
    const own = stat.codon_aa.get(codon);
    const m = new Map();
    let tot = 0;
    for (let i = 0; i < AA_ORDER.length; i++) {
      const alt = AA_ORDER[i];
      if (!(d[i] > 0)) continue;
      // Every synonymous route lands in one block: what they share is that the protein does
      // not change, not the residue they happen to preserve.
      const key = alt === own ? SYNONYMOUS : alt;
      m.set(key, (m.get(key) || 0) + d[i]);
      tot += d[i];
    }
    if (m.size && tot > 0) raw.set(codon, { m, tot });
  }
  if (!raw.size) {
    const p = document.createElement("p");
    p.className = "caveat";
    p.textContent = "No residue change is reachable from the codons in this batch in one step.";
    container.appendChild(p);
    captionEl.textContent = "";
    return false;
  }

  const flows = new Map();          // codon -> Map(destination -> band value)
  const supplyOf = new Map();       // codon -> share basis for the "% of supply" label
  for (const [codon, { m, tot }] of raw) {
    const f = freqOf.get(codon);
    const out = new Map();
    for (const [alt, rate] of m) out.set(alt, mode === "supply" ? f * rate : f * rate / tot);
    flows.set(codon, out);
    supplyOf.set(codon, f * tot);   // always the supply basis, whichever mode is drawn
  }

  // Pool the rare tail. The most common codon is kept whatever its frequency, so a batch of
  // uniformly rare codons still shows one named block rather than a single anonymous one.
  const minFreq = opts.minCodonFreq == null ? 0.02 : opts.minCodonFreq;
  const keepAnyway = codons[0];
  const rare = [...flows.keys()].filter((c) => freqOf.get(c) < minFreq && c !== keepAnyway);
  let pooledCodons = 0, pooledFreq = 0;
  if (rare.length > 1) {
    const merged = new Map();
    let mergedSupply = 0;
    for (const c of rare) {
      for (const [a, v] of flows.get(c)) merged.set(a, (merged.get(a) || 0) + v);
      pooledFreq += freqOf.get(c);
      mergedSupply += supplyOf.get(c);
      flows.delete(c); supplyOf.delete(c);
    }
    pooledCodons = rare.length;
    flows.set(OTHER_CODONS, merged);
    freqOf.set(OTHER_CODONS, pooledFreq);
    supplyOf.set(OTHER_CODONS, mergedSupply);
  }

  const dTotal = new Map();
  for (const m of flows.values()) {
    for (const [a, v] of m) dTotal.set(a, (dTotal.get(a) || 0) + v);
  }
  // Keep the largest destinations and pool the tail, or the bands become hairlines.
  const limit = opts.limit || 8;
  let order = [...dTotal.entries()].sort((a, b) => b[1] - a[1]).map(([a]) => a);
  // The no-change block always keeps its own row: it is the one destination a reader needs to
  // see whatever its size, and pooling it into "other" would merge it with residue changes.
  const hasSyn = order.includes(SYNONYMOUS);
  order = order.filter((a) => a !== SYNONYMOUS);
  const tail = order.slice(limit - (hasSyn ? 1 : 0));
  order = order.slice(0, limit - (hasSyn ? 1 : 0));
  if (hasSyn) order.unshift(SYNONYMOUS);
  if (tail.length) {
    const pooled = new Set(tail);
    for (const m of flows.values()) {
      let o = 0;
      for (const a of pooled) { if (m.has(a)) { o += m.get(a); m.delete(a); } }
      if (o > 0) m.set(OTHER, (m.get(OTHER) || 0) + o);
    }
    dTotal.set(OTHER, tail.reduce((t, a) => t + dTotal.get(a), 0));
    for (const a of tail) dTotal.delete(a);
    order.push(OTHER);
  }
  const total = [...dTotal.values()].reduce((t, v) => t + v, 0);
  const cTotal = new Map([...flows].map(([c, m]) =>
    [c, [...m.values()].reduce((t, v) => t + v, 0)]));
  const supplyTotal = [...supplyOf.values()].reduce((t, v) => t + v, 0);
  const live = codons.filter((c) => flows.has(c) && cTotal.get(c) > 0)
    .concat(flows.has(OTHER_CODONS) && cTotal.get(OTHER_CODONS) > 0 ? [OTHER_CODONS] : []);

  // A codon is drawn in the colour of the residue it encodes, so degenerate codons share one
  // colour. The pooled block spans several residues and takes the neutral colour.
  const bandColour = (c) =>
    c === OTHER_CODONS ? NEUTRAL_COLOUR : residueColour(stat.codon_aa.get(c));

  const width = opts.width || LAYOUT.width;
  // A drawn title needs room above the first block. The dashboard keeps the title in the
  // card heading (drawTitle: false) and only needs enough inset that the top label is not
  // clipped by the viewBox.
  const showTitle = !!(opts.title && opts.drawTitle !== false);
  const padTop = (showTitle && opts.subtitle !== false) ? 34
               : showTitle ? 22
               : (opts.title ? 20 : 12);
  const padBottom = 14;
  const height = opts.height || Math.max(260, 42 * Math.max(live.length, order.length) + 70);
  const plotH = height - padTop - padBottom;
  // Label gutters, sized to what the labels actually need rather than to a round number. The
  // widest left label over the whole compendium is "ATG (M) <" at about 60 px and the widest
  // right one "synonymous  49%" at about 90 px; the gutters were 320 and 290, roughly three
  // times that, which is what squeezed the bands out of the chart when it was given less than
  // the full page. Everything the labels do not need goes to the bands, because the bands are
  // the chart -- the columns on their own are a table.
  const BW = 16;
  const GUT_L = 104, GUT_R = 116;      // label width plus its 10 px standoff, rounded up
  const MIN_BAND = 150;
  // Below about 400 px the gutters have to give way too, or they would meet in the middle and
  // the bands would be drawn with negative width -- two columns and nothing between them.
  const need = GUT_L + GUT_R + BW * 2 + MIN_BAND;
  const squeeze = width < need ? (width - BW * 2 - MIN_BAND) / (GUT_L + GUT_R) : 1;
  const gutL = Math.max(0, Math.round(GUT_L * Math.min(1, squeeze)));
  const gutR = Math.max(0, Math.round(GUT_R * Math.min(1, squeeze)));
  const xL = gutL;
  const xR = width - gutR;
  const NS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("width", "100%");
  svg.setAttribute("role", "img");
  if (opts.title) svg.setAttribute("aria-label", opts.title);
  const mk = (n, a, t) => {
    const e = document.createElementNS(NS, n);
    for (const [k, v] of Object.entries(a)) e.setAttribute(k, v);
    if (t != null) e.textContent = t;
    return e;
  };
  const stack = (keys, totals) => {
    const gap = 9;
    const usable = plotH - gap * Math.max(0, keys.length - 1);
    const out = new Map();
    let y = padTop;
    for (const k of keys) {
      const h = Math.max(2, (totals.get(k) / total) * usable);
      out.set(k, [y, y + h]);
      y += h + gap;
    }
    return out;
  };
  const cy = stack(live, cTotal), dy = stack(order, dTotal);

  // Title, subtitle and the two column headings are drawn only where the caller asks for a
  // title. A figure with its own heading and caption does not need them repeated inside the
  // plot, and on a narrow figure they are the first thing to collide; the caption carries
  // them instead. A caller that passes a title and leaves subtitle on gets the full set.
  if (showTitle) {
    svg.appendChild(mk("text", { x: 0, y: 13, "font-size": 12 }, opts.title));
  }
  // The subtitle and the two column headings are the verbose half of the header. A figure
  // that sits beside another one, under a shared caption, wants only enough to tell the two
  // apart; the caption says what the blocks mean.
  if (showTitle && opts.subtitle !== false) {
    svg.appendChild(mk("text", { x: 0, y: 27, "font-size": 9.5, fill: "#5c6670" },
                       mode === "supply"
                         ? "Sized by mutational supply: a block is its share of the change this "
                           + "position generates."
                         : "Sized by how common each codon is: a block is its share of the "
                           + "sequences."));
    svg.appendChild(mk("text", { x: xL + BW, y: padTop - 5, "text-anchor": "end",
                                 "font-size": 9, fill: "#5c6670" }, "codon carried"));
    svg.appendChild(mk("text", { x: xR - BW, y: padTop - 5, "font-size": 9, fill: "#5c6670" },
                       "residue it becomes"));
  }

  const cCur = new Map(live.map((c) => [c, cy.get(c)[0]]));
  const dCur = new Map(order.map((a) => [a, dy.get(a)[0]]));
  live.forEach((c, ci) => {
    const m = flows.get(c);
    for (const a of order) {
      const v = m.get(a);
      if (!(v > 0)) continue;
      const hL = (v / cTotal.get(c)) * (cy.get(c)[1] - cy.get(c)[0]);
      const hR = (v / dTotal.get(a)) * (dy.get(a)[1] - dy.get(a)[0]);
      const y0 = cCur.get(c), y1 = dCur.get(a);
      cCur.set(c, y0 + hL); dCur.set(a, y1 + hR);
      const xa = xL + BW, xb = xR - BW, xm = (xa + xb) / 2;
      const d = `M${xa},${y0.toFixed(1)} C${xm},${y0.toFixed(1)} ${xm},${y1.toFixed(1)} ` +
                `${xb},${y1.toFixed(1)} L${xb},${(y1 + hR).toFixed(1)} ` +
                `C${xm},${(y1 + hR).toFixed(1)} ${xm},${(y0 + hL).toFixed(1)} ` +
                `${xa},${(y0 + hL).toFixed(1)} Z`;
      svg.appendChild(mk("path", { d, fill: bandColour(c), "fill-opacity": 0.42 }));
    }
  });

  live.forEach((c, ci) => {
    const [y0, y1] = cy.get(c);
    svg.appendChild(mk("rect", { x: xL, y: y0.toFixed(1), width: BW,
                                 height: (y1 - y0).toFixed(1),
                                 fill: bandColour(c), "fill-opacity": 0.9 }));
    const mid = (y0 + y1) / 2;
    const isPool = c === OTHER_CODONS;
    svg.appendChild(mk("text", { x: xL - 10, y: mid - 2, "text-anchor": "end", "font-size": 11,
                                 "font-family": isPool ? "inherit" : "monospace" },
                       isPool ? `${pooledCodons} rarer codons`
                              : `${c} (${stat.codon_aa.get(c)})` +
                                // The marker alone: a reader who needs telling what the
                                // triangle means has the caption, and spelling it out on
                                // every chart cost 50 px of gutter on every one of them.
                                (c === opts.refCodon ? " \u25c0" : "")));
    svg.appendChild(mk("text", { x: xL - 10, y: mid + 11, "text-anchor": "end", "font-size": 9,
                                 fill: "#5c6670" },
                       `${(100 * freqOf.get(c)).toFixed(0)}%`));
  });
  // Thin blocks sit closer together than a line of type. Keep each label on its block
  // when there is room, and step it down only far enough that two lines never share pixels.
  let labelY = padTop + 9;
  for (const a of order) {
    const [y0, y1] = dy.get(a);
    const dCol = (a === OTHER || a === SYNONYMOUS) ? NEUTRAL_COLOUR : residueColour(a);
    svg.appendChild(mk("rect", { x: xR - BW, y: y0.toFixed(1), width: BW,
                                 height: (y1 - y0).toFixed(1), fill: dCol,
                                 "fill-opacity": 0.75 }));
    const label = a === OTHER ? "other" : a === SYNONYMOUS ? "synonymous"
                : a === "*" ? "stop" : a;
    const y = Math.max((y0 + y1) / 2 + 4, labelY);
    // A label that cannot clear the one above it without leaving the box is omitted.
    // Drawing it anyway stacks the smallest percentages on one line.
    if (y > height - 4) continue;
    svg.appendChild(mk("text", { x: xR + 10, y: y.toFixed(1), "font-size": 11 },
                       `${label}  ${(100 * dTotal.get(a) / total).toFixed(0)}%`));
    labelY = y + 13;
  }
  container.appendChild(svg);
  if (typeof opts.onLegend === "function") {
    const refAa = (opts.refCodon && stat.codon_aa.get(opts.refCodon)) || opts.refAa || "";
    opts.onLegend(order.map((a) => ({
      label: a === OTHER ? "other"
           : a === SYNONYMOUS ? (refAa ? `synonymous (${refAa})` : "synonymous")
           : a === "*" ? "stop" : a,
      colour: (a === OTHER || a === SYNONYMOUS) ? NEUTRAL_COLOUR : residueColour(a),
      percent: total > 0 ? 100 * dTotal.get(a) / total : 0,
    })));
  }

  let lead = null;
  for (const c of live) {
    if (c === OTHER_CODONS) continue;      // a pooled block is not a finding about one codon
    const f = freqOf.get(c), share = supplyOf.get(c) / supplyTotal;
    if (!lead || share - f > lead.gap) lead = { codon: c, f, share, gap: share - f };
  }
  // describeColumns is Browse and Your Sequence. Those pages size the two columns
  // differently: the left is the LANL frequency, the right is that frequency times each
  // codon's pie. Callers that omit it keep the older supply-versus-frequency wording.
  const opening = opts.describeColumns
    ? `Left blocks are the LANL codon frequency only. Right blocks are relative supply: ` +
      `codon frequency multiplied by that codon's fraction of its one-step changes. `
    : (mode === "supply"
      ? `Band width is a codon's frequency multiplied by its rate to that residue, so a left ` +
        `block is its share of the <em>supply</em> rather than of the sequences. `
      : `Band width is how common the codon is, split by where that codon's single ` +
        `substitutions lead. Left blocks are shares of <em>sequences</em>; switch the ` +
        `weighting to size them by supply instead. `);
  const comparison = opts.describeColumns ? ""
    : (lead && lead.gap > 0.05
      ? `<strong>${lead.codon}</strong> is carried by ${(100 * lead.f).toFixed(0)}% of these ` +
        `sequences but accounts for <strong>${(100 * lead.share).toFixed(0)}%</strong> of what ` +
        `this position generates. `
      : `Supply is spread across the codons roughly in proportion to how common they are. `);
  captionEl.innerHTML = opening + comparison +
    // Deliberately the supply fraction from siteSupplyStats, not this chart's own synonymous
    // band share: in "sequences" weighting a band is a per-codon proportion, so its share is a
    // mean of proportions and would read differently from the same quantity in the table and
    // the readout. One number for one thing, whichever weighting is on screen.
    (opts.stats
      ? `<strong>${(100 * opts.stats.syn_fraction).toFixed(0)}%</strong> of the supply here is ` +
        `synonymous and never reaches the protein. `
      : "") +
    (opts.stats
      ? `Across the residue changes, the supply is spread over the equivalent of ` +
        `<strong>${opts.stats.effective_destinations.toFixed(1)}</strong> residues ` +
        `(exp of the entropy; ${opts.stats.n_destinations} are reachable at all), with ` +
        `${(100 * opts.stats.top_share).toFixed(0)}% of it going to ` +
        `<strong>${opts.stats.top_destination}</strong>. `
      : "") +
    (tail.length ? ` The ${tail.length} smallest destinations are pooled as "other".` : "") +
    (pooledCodons
      ? ` ${pooledCodons} codons each under ${(100 * minFreq).toFixed(0)}% of the batch are ` +
        `pooled into one block (${(100 * pooledFreq).toFixed(0)}% of sequences between them); ` +
        `the per-position table lists them separately.`
      : "");
  return true;
}

const OTHER = "\u2014other\u2014";
const SYNONYMOUS = "\u2014synonymous\u2014";
const OTHER_CODONS = "\u2014other codons\u2014";
