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

import { genomeTrack, pieChart, emergenceCurve, LAYOUT } from "./chart.js";
import { smoothSeries, oneStepSlices, emergenceStats, emergenceProbability }
  from "./explore_core.js";
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
    width: o.width || LAYOUT.width, height: o.height || 320,
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
  if (o.statsEl) renderTrackStats(o.statsEl, trackStats({
    vals, med, at, site, cls: o.cls, win: o.win, marker: o.marker, refBase: o.refBase,
  }));

  return `Centred running mean of the per-site rate over a ${o.win} nt window, summed across ` +
    `substitution classes${o.cls === "total" ? "" : " and restricted to " + CLASS_LABEL[o.cls]}. ` +
    `Rates are ${o.source}. Insertions and deletions are excluded. The band above the plot ` +
    `shows coding segments; overlapping reading frames sit on separate rows, and spliced ` +
    `genes (tat, rev) appear as one box per exon, with the gene holding the marker filled ` +
    `solid.` + spanNote;
}

/**
 * The quantity rows shown under a track. Pure, so test/explore_cases.mjs can assert what each
 * one is supposed to do -- in particular that the median is a FIXED reference across the whole
 * track while the relative row MOVES with the marker. An earlier version labelled the median
 * "over the plotted window", which read as if it were local to the marker; it never was.
 *
 * @param {object} a
 *   vals     every plotted value, sorted ascending
 *   med      the median of vals
 *   at       [position, value] at the marker, or undefined
 *   site     {total, ga, nonga} at the exact marker position, or undefined
 *   cls, win, marker, refBase
 */
export function trackStats(a) {
  let rank = null;
  if (a.at && a.vals.length) {
    let below = 0;
    for (const v of a.vals) { if (v <= a.at[1]) below++; else break; }   // vals is sorted
    rank = 100 * below / a.vals.length;
  }
  return [
    // Deliberately NOT a median of the window around the marker: it is taken over every
    // plotted point, so it is a fixed reference line for the chosen class and window and
    // does not move when the marker does. That is what makes the row below meaningful.
    { label: "Median across the whole track",
      value: a.med === null || a.med === undefined ? "\u2014" : formatRate(a.med),
      note: `all ${a.vals.length.toLocaleString()} plotted positions, ${CLASS_LABEL[a.cls]}, ` +
            `${a.win} nt window \u2014 fixed while you move the marker` },
    { label: "Value at the marker",
      value: a.at ? formatRate(a.at[1]) : "\u2014",
      note: `smoothed window centred on nt ${a.marker.toLocaleString()}` },
    { label: "Marker relative to the track",
      value: rank === null ? "\u2014" : `${(a.at[1] / a.med).toFixed(2)}\u00d7`,
      note: rank === null ? "no value at this position"
        : `${rank.toFixed(0)}th percentile of the plotted track` },
    { label: "G\u2192A share at this exact site",
      value: a.site ? `${(100 * a.site.ga / a.site.total).toFixed(0)}%` : "\u2014",
      // A non-G base cannot have G->A supply, so 0% there is a property of the base rather
      // than a small measurement. Saying which base it is stops that reading as a bug.
      note: !a.site ? "site not covered"
        : a.refBase && a.refBase !== "G"
          ? `0 by definition \u2014 the base here is ${a.refBase}, not G`
          : `unsmoothed, nt ${a.marker.toLocaleString()} only` },
  ];
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
 * The codon panel: a pie of the one-nucleotide-accessible amino-acid changes, a selector for
 * one of them, and -- for whichever is selected -- the expected number of de novo copies as a
 * function of the replicating population size.
 *
 * The pie answers "where does supply go at this codon?"; the curve answers "how big does the
 * population have to be before I expect to see this particular change at all?". The second
 * question is the one that connects a rate to something observable, which is why the two sit
 * side by side and share a selection.
 *
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

  // The pie and the curve together span exactly the track's width, so the three figures on
  // the page share one left and right edge. The pie takes the larger share because its legend
  // is text that must not wrap.
  const PIE_W = Math.round((LAYOUT.width - LAYOUT.gap) * 0.56);
  const CURVE_W = LAYOUT.width - LAYOUT.gap - PIE_W;
  const FIG_H = 460;

  const keyOf = (r) => `${r.alt_aa}${r.synonymous ? "*" : ""}`;
  // Two rows could share a target amino acid only via different codons, which cannot happen
  // within one codon; still, keep the first so the map is total.
  const byKey = new Map();
  for (const r of oneStep) if (!byKey.has(keyOf(r))) byKey.set(keyOf(r), r);
  let selectedKey = keyOf(oneStep[0]);            // the largest slice, already sorted

  // --- controls ---------------------------------------------------------------------
  const ctrl = document.createElement("div");
  ctrl.className = "controls";
  const lab = document.createElement("label");
  lab.textContent = "Amino-acid change";
  const sel = document.createElement("select");
  for (const r of oneStep) {
    const o = document.createElement("option");
    o.value = keyOf(r);
    o.textContent = `${r.ref_aa}\u2192${r.alt_aa}${r.synonymous ? " (synonymous)" : ""}` +
                    `  \u2014  ${r.mu_aa.toExponential(2)}`;
    sel.appendChild(o);
  }
  const wrap = document.createElement("div");
  wrap.append(lab, sel);
  ctrl.appendChild(wrap);
  chartEl.appendChild(ctrl);

  // --- side-by-side figures -----------------------------------------------------------
  const row = document.createElement("div");
  row.className = "panel-row";
  const pieSlot = document.createElement("div");
  const curveSlot = document.createElement("div");
  // Explicit bases: a flex item holding an SVG at width:100% has no intrinsic width and
  // would otherwise collapse, which is what made these figures render far smaller than asked.
  pieSlot.style.flex = `1 1 ${PIE_W}px`;
  pieSlot.style.maxWidth = `${PIE_W}px`;
  curveSlot.style.flex = `1 1 ${CURVE_W}px`;
  curveSlot.style.maxWidth = `${CURVE_W}px`;
  row.append(pieSlot, curveSlot);
  chartEl.appendChild(row);
  const readout = document.createElement("div");
  chartEl.appendChild(readout);

  function pick(key) {
    if (!byKey.has(key)) return;
    selectedKey = key;
    sel.value = key;
    redraw();
  }

  function redraw() {
    pieChart(pieSlot, oneStep.map((r) => ({
      label: keyOf(r), value: r.mu_aa, highlight: false,
    })), { title: `Single-nucleotide amino-acid changes at ${label}`,
           width: opts.width || PIE_W, height: opts.height || FIG_H,
           selected: selectedKey, onSelect: pick });

    const r = byKey.get(selectedKey);
    emergenceCurve(curveSlot, {
      mu: r.mu_aa, lo95: r.mu_aa_lo95, hi95: r.mu_aa_hi95,
      change: `${label} ${r.ref_aa}\u2192${r.alt_aa}`,
      // Out to 1e8 rather than 1e6: the transition sits at 1e3-1e6, but the saturation
      // regime is where a real infection lives, and cutting the axis at 1e6 hides the
      // conclusion that one-step changes are not supply-limited in vivo.
      nMin: 1, nMax: 1e8, width: CURVE_W, height: FIG_H,
    });

    const { nStar, nLo, nHi } = emergenceStats(r);
    const pAt = (n) => emergenceProbability(r.mu_aa, n);
    const pct = (x) => x >= 0.9995 ? "&gt;99.9%" : `${(100 * x).toFixed(1)}%`;
    readout.innerHTML =
      '<table class="data track-stats"><thead><tr><th>Quantity</th><th>Value</th>' +
      "<th>Basis</th></tr></thead><tbody>" +
      `<tr><td>Rate of ${r.ref_aa}\u2192${r.alt_aa}</td>` +
      `<td class="num">${r.mu_aa.toExponential(2)}</td>` +
      `<td>per replication event; 95% CI ${r.mu_aa_lo95.toExponential(2)}\u2013` +
      `${r.mu_aa_hi95.toExponential(2)}</td></tr>` +
      `<tr><td>N* \u2014 replication events for one expected copy</td>` +
      `<td class="num">${fmt(nStar)}</td>` +
      `<td>1/rate; 95% CI ${fmt(nLo)}\u2013${fmt(nHi)}. At N* the chance of at least one ` +
      `is 63%, not certainty</td></tr>` +
      `<tr><td>P(generated) at N = 10\u2075</td><td class="num">${pct(pAt(1e5))}</td>` +
      `<td>a small or recently founded population</td></tr>` +
      `<tr><td>P(generated) at N = 10\u2077</td><td class="num">${pct(pAt(1e7))}</td>` +
      `<td>order of a chronic untreated infection, per cycle</td></tr>` +
      "</tbody></table>";
  }

  sel.addEventListener("change", () => pick(sel.value));
  redraw();

  captionEl.innerHTML =
    "<strong>Left:</strong> each slice is an amino-acid change reachable from " +
    (opts.custom ? "<strong>your</strong> codon" : "the reference codon") +
    " by a <strong>single</strong> nucleotide substitution, sized by its share of the one-step " +
    "mutational supply at this position. <code>*</code> marks a synonymous change. Click a " +
    "slice, or use the selector, to choose one. " +
    "<strong>Right:</strong> the probability that the selected change is <em>generated at " +
    "least once</em> in a single replication cycle, in a population performing N replication " +
    "events \u2014 P = 1 \u2212 (1 \u2212 rate)<sup>N</sup>, treating events as independent " +
    "trials. The band comes from the rate's 95% interval; N* marks where one copy is expected, " +
    "which is a 63% chance, not a certainty. <strong>N counts replication events</strong> " +
    "(population size \u00d7 the fraction replicating), so it is not a fraction and is not " +
    "bounded by 1. " +
    (opts.custom
      ? "Rates are computed for the codon in your sequence, so they differ from the published " +
        "NL4-3 atlas wherever your codon differs. "
      : "Rates are read from the atlas, not recomputed. ") +
    "<strong>This is generation, not establishment:</strong> a variant generated here may " +
    "still be lost to drift, purged by selection, or never sampled \u2014 the curve speaks " +
    "only to mutational supply. It is also per cycle, so it understates supply accumulated " +
    "over an infection, and it applies only to the one-step changes in the pie; a multi-step " +
    "change is not a single-cycle event and must not be read off it. The shaded band marks " +
    "the order of magnitude of replication in a chronic untreated infection, taken from the " +
    "HIV dynamics literature rather than measured here.";
  return true;
}
