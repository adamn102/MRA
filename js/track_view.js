// SPDX-License-Identifier: MIT
//
// track_view.js -- the genome supply track and the position pie, shared by the published
// atlas (browse.html) and a submitted sequence (sequence.html).
//
// Both pages show the same two things: a smoothed per-site supply track along the genome,
// and a pie of where one substitution's supply goes at the selected codon. The only
// differences are the source of the numbers (published atlas vs a custom atlas built from
// the submitted sequence) and, for a submitted sequence, the variant ticks. Keeping one
// renderer here is what stops the two views drifting apart.

import { genomeTrack, pieChart } from "./chart.js";
import { smoothSeries, oneStepSlices } from "./explore_core.js";
import { formatRate } from "./table.js";

export const CLASS_LABEL = {
  total: "total substitution supply",
  nonga: "non-G\u2192A supply",
  ga: "G\u2192A supply",
};
const CLASS_COLOUR = { total: "#2f4b7c", nonga: "#3f8f7a", ga: "#c2561a" };

/**
 * Draw the supply track and return the caption HTML describing exactly what is plotted.
 *
 * @param {HTMLElement} container
 * @param {object} o
 *   perSite    Map<ntPosition, {total, ga, nonga}>
 *   sites      sorted array of the positions in perSite
 *   win        smoothing window, in sites
 *   cls        "total" | "nonga" | "ga"
 *   segments   gene annotation band
 *   first,last x-axis bounds (the modelled region, not the covered span)
 *   marker     selected nucleotide
 *   markerLabel
 *   variants   optional array of positions differing from NL4-3
 *   onSelect   callback(position)
 *   source     short phrase naming where the rates came from
 */
export function drawSupplyTrack(container, o) {
  const pts = smoothSeries(o.sites, o.perSite, o.cls, o.win);
  genomeTrack(container, [
    { label: CLASS_LABEL[o.cls], colour: CLASS_COLOUR[o.cls], points: pts, width: 1.1 },
  ], {
    title: o.title || `Smoothed ${CLASS_LABEL[o.cls]} along the genome`,
    xMin: o.first, xMax: o.last, geneBand: o.segments,
    width: o.width || 1020, height: o.height || 300,
    yLabel: "supply per site per replication",
    marker: o.marker, markerLabel: o.markerLabel,
    variants: o.variants, onSelect: o.onSelect,
  });

  const vals = pts.map((p) => p[1]).sort((a, b) => a - b);
  const med = vals.length ? vals[Math.floor(vals.length / 2)] : null;
  const at = pts.find((p) => p[0] >= o.marker);
  const site = o.perSite.get(o.marker);
  // The covered span matters on a submitted sequence: a partial amplicon leaves most of the
  // axis empty, and saying so is clearer than letting the reader infer it from a gap.
  const covered = o.sites.length;
  const spanNote = covered < (o.last - o.first + 1)
    ? ` Only ${covered.toLocaleString()} of ${(o.last - o.first + 1).toLocaleString()} ` +
      `reference positions are covered by this sequence; the rest of the axis is blank.`
    : "";
  return `Centred running mean of the per-site rate over a ${o.win} nt window, summed across ` +
    `substitution classes${o.cls === "total" ? "" : " and restricted to " + CLASS_LABEL[o.cls]}. ` +
    `Rates are ${o.source}. Insertions and deletions are excluded. The band above the plot ` +
    `shows coding segments; overlapping reading frames sit on separate rows, and spliced ` +
    `genes (tat, rev) appear as one box per exon, with the gene holding the marker filled ` +
    `solid.` + spanNote +
    (med ? ` Median over the plotted window ${formatRate(med)}` : "") +
    (at ? `, value at the marker ${formatRate(at[1])}` : "") +
    (site ? `, of which G\u2192A is ${(100 * site.ga / site.total).toFixed(0)}% at this exact site` : "") +
    `. <strong>G\u2192A carries most of the total supply</strong>, so the total track partly ` +
    `traces G density \u2014 use the class selector to separate the two.`;
}

/**
 * Pie of the one-nucleotide-accessible amino-acid changes at a codon, with its caption.
 * Returns true if a pie was drawn, false if no single-step change exists from this codon.
 */
export function drawPositionPie(chartEl, captionEl, rows, label, opts = {}) {
  const oneStep = oneStepSlices(rows);
  chartEl.innerHTML = "";
  if (!oneStep.length) {
    captionEl.innerHTML =
      "No amino-acid change is reachable from this codon by a single nucleotide substitution; " +
      "see the table below for the 2- and 3-step changes.";
    return false;
  }
  pieChart(chartEl, oneStep.map((r) => ({
    label: `${r.alt_aa}${r.synonymous ? "*" : ""}`, value: r.mu_aa, highlight: false,
  })), { title: `Single-nucleotide amino-acid changes at ${label}`, width: opts.width || 680 });
  captionEl.innerHTML =
    "Each slice is an amino-acid change reachable from " +
    (opts.custom ? "<strong>your</strong> codon" : "the reference codon") +
    " by a <strong>single</strong> nucleotide substitution, sized by its share of the one-step " +
    "mutational supply at this position. <code>*</code> marks a synonymous change. The shares " +
    "are a property of the whole codon, not of the single nucleotide marked on the track. " +
    "Two- and three-step changes are in the table below but not in the pie, because they are " +
    "generated far more rarely. " +
    (opts.custom
      ? "Values are computed for the codon in your sequence, so they differ from the published " +
        "NL4-3 atlas wherever your codon differs."
      : "Values are read from the atlas, not recomputed.");
  return true;
}
