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
 *   refBase    optional base at the marker, so a 0% G→A row can explain itself
 *   statsEl    optional element to receive the quantity table
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
  // The three figures that used to run on as a sentence. A table because they are read as
  // values, not as prose: a reader comparing the marker against the window median wants them
  // aligned under each other, not separated by commas.
  if (o.statsEl) {
    renderTrackStats(o.statsEl, [
      { label: `Median over the plotted window`,
        value: med === null ? "\u2014" : formatRate(med),
        note: `${CLASS_LABEL[o.cls]}, ${o.win} nt window` },
      { label: `Value at the marker`,
        value: at ? formatRate(at[1]) : "\u2014",
        note: `smoothed window centred on nt ${o.marker.toLocaleString()}` },
      { label: `G\u2192A share at this exact site`,
        value: site ? `${(100 * site.ga / site.total).toFixed(0)}%` : "\u2014",
        // A non-G base cannot have G->A supply, so 0% there is a property of the base rather
        // than a small measurement. Saying which base it is stops that reading as a bug.
        note: !site ? "site not covered"
          : o.refBase && o.refBase !== "G"
            ? `0 by definition — the base here is ${o.refBase}, not G`
            : `unsmoothed, nt ${o.marker.toLocaleString()} only` },
    ]);
  }

  return `Centred running mean of the per-site rate over a ${o.win} nt window, summed across ` +
    `substitution classes${o.cls === "total" ? "" : " and restricted to " + CLASS_LABEL[o.cls]}. ` +
    `Rates are ${o.source}. Insertions and deletions are excluded. The band above the plot ` +
    `shows coding segments; overlapping reading frames sit on separate rows, and spliced ` +
    `genes (tat, rev) appear as one box per exon, with the gene holding the marker filled ` +
    `solid.` + spanNote;
}

/** The small quantity table under the track. */
function renderTrackStats(el, rows) {
  el.innerHTML =
    '<table class="data track-stats"><thead><tr>' +
    "<th>Quantity</th><th>Value</th><th>Basis</th></tr></thead><tbody>" +
    rows.map((r) => `<tr><td>${r.label}</td><td class="num">${r.value}</td>` +
                    `<td>${r.note}</td></tr>`).join("") +
    "</tbody></table>";
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
