// SPDX-License-Identifier: MIT
//
// track_view.js -- the genome supply track and the position pie, shared by the published
// atlas (the landing page, index.html) and a submitted sequence (sequence.html).
//
// Both pages show the same two things: a smoothed per-site supply track along the genome,
// and a pie of where one substitution's supply goes at the selected codon. The only
// differences are the source of the numbers (published atlas vs a custom atlas built from
// the submitted sequence) and, for a submitted sequence, the variant ticks. Keeping one
// renderer here is what stops the two views drifting apart.

import { genomeTrack, pieChart, LAYOUT } from "./chart.js";
import { residueColour, drawSupplyFlow } from "./batch_view.js";
import { smoothSeries, oneStepSlices }
  from "./explore_core.js";
import { formatRate } from "./table.js";
import { setDetails } from "./figure_details.js";

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
  const series = [
    { label: o.label || CLASS_LABEL[o.cls], colour: CLASS_COLOUR[o.cls], points: pts, width: 1.1 },
  ];
  // NL4-3 OVERLAY. Drawn on exactly the sites the submission covers, not on the whole
  // reference: the smoothing window must contain the same positions in both series or the
  // two lines are not comparable point for point, and a partial amplicon would otherwise be
  // set against a whole-genome curve. Dashed and grey so the submission reads as the subject
  // and the reference as the comparator.
  let refPts = null;
  if (o.refPerSite) {
    const shared = o.sites.filter((p) => o.refPerSite.has(p));
    if (shared.length >= 3) {
      refPts = smoothSeries(shared, o.refPerSite, o.cls, o.win);
      series.push({ label: "NL4-3", colour: "#8a94a0", points: refPts, width: 1.0,
                    dash: "4 3", opacity: 0.95 });
    }
  }
  genomeTrack(container, series, {
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
  if (o.statsEl) {
    const rows = trackStats({
      vals, med, at, site, cls: o.cls, win: o.win, marker: o.marker, refBase: o.refBase,
    });
    // Where the two lines sit close together the eye cannot separate them, so the comparison
    // is also given as a number at the marker.
    if (refPts) {
      const refAt = refPts.find((p) => p[0] >= o.marker);
      if (at && refAt && refAt[1] > 0) {
        const ratio = at[1] / refAt[1];
        rows.push({ label: "NL4-3 at the marker", value: formatRate(refAt[1]),
                    note: "the reference's own rate over the same positions and window" });
        rows.push({
          label: "Your sequence \u00f7 NL4-3",
          value: (ratio >= 1 ? "\u00d7" : "\u00f7") +
                 (ratio >= 1 ? ratio : 1 / ratio).toFixed(2),
          note: Math.abs(Math.log2(ratio)) < 0.05
            ? "within 3% \u2014 the codons and contexts here are effectively the reference's"
            : (ratio > 1 ? "this window generates more change than NL4-3 does"
                         : "this window generates less change than NL4-3 does"),
        });
      }
    }
    renderTrackStats(o.statsEl, rows);
  }

  return `Centred running mean of the per-site rate over a ${o.win} nt window, summed across ` +
    `substitution classes${o.cls === "total" ? "" : " and restricted to " + CLASS_LABEL[o.cls]}. ` +
    `Rates are ${o.source}. Insertions and deletions are excluded. ` +
    (refPts
      ? `The dashed grey line is <strong>NL4-3</strong> over the same positions and the same ` +
        `window, so the two are comparable point for point. They coincide wherever your ` +
        `codons and 3\u2032 contexts match the reference, which for a close relative is most ` +
        `of the genome \u2014 the gaps between them are the whole of what your sequence ` +
        `changes about mutational supply. `
      : "") +
    `The band above the plot shows coding segments; overlapping reading frames sit on ` +
    `separate rows, and spliced genes (tat, rev) appear as one box per exon, with the gene ` +
    `holding the marker filled solid.` + spanNote;
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

/**
 * The flow card's Details text: what the two columns are, then whatever the chart itself
 * reported about this position. The chart writes a caption describing its own contents --
 * pooling, totals, how many sequences -- and that is worth keeping rather than rewriting
 * here, where it would drift from what was actually drawn.
 */
function flowDetailHtml(label, chartCaption) {
  return "<p>The same position across circulating virus, from a batch over the 182 LANL " +
    "group M compendium genomes. <strong>Left column:</strong> the codons those sequences " +
    "carry at " + label + ", with the share carrying each; \u25c0 marks the codon NL4-3 " +
    "has. <strong>Right column:</strong> the residue each becomes. <strong>Band width</strong> " +
    "is that route's share of the residue change this position generates.</p>" +
    "<p>The pie beside it is one codon's draw; this is the population's. Both are coloured " +
    "by destination residue, so a slice and a band for the same residue match.</p>" +
    (chartCaption ? "<p>" + chartCaption + "</p>" : "");
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
 * The codon panel: a pie of the one-nucleotide-accessible amino-acid changes at this codon,
 * and beside it -- where a flow source is supplied -- where residue change at this position
 * is actually generated across circulating virus.
 *
 * WHAT REPLACED THE EMERGENCE CURVE. This slot used to hold "probability of generating this
 * change de novo" against replicating population size. That curve plotted
 * 1 - (1 - mu)^N, which is a restatement of the rate the pie already shows against an N the
 * atlas cannot measure -- so its shape came from the reader's assumption about population
 * size, not from the data. The flow chart answers a question the atlas can actually speak
 * to: given the codons circulating virus carries here, where does residue change get
 * generated? It is the same chart the batch tab draws, from a batch over the LANL group M
 * compendium precomputed at build time.
 *
 * The two figures share a colour key -- both are coloured by destination residue, and
 * degenerate codons for one residue take that residue's colour -- so a slice and a band can
 * be matched by eye.
 *
 * @param {object} [opts.flow]  {stat, matrix} from lanl_flow.js, or null for the pie alone
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

  // Side by side, together spanning exactly the track's width, so the dashboard and the
  // genome track above it share one left and one right edge.
  //
  // This only works because the flow chart's label gutters are now sized to the labels
  // (~104 px and ~116 px) instead of the 320/290 they used to reserve. At this split the
  // bands get ~350 px, against a negative width -- literally no bands -- when half the page
  // was tried against the old gutters. The pie takes slightly less than half: it is the
  // simpler figure and its legend is text that must not wrap.
  const hasFlow = !!(opts.flow && opts.flow.stat && opts.flow.matrix);
  // Two equal boxes spanning exactly the track's width. Equal matters here: these are two
  // views of one position, and a wider box reads as the more important figure when neither
  // is. Both are given the same height as well, so the row is a pair of matching panels
  // rather than a tall chart next to a short one.
  const PIE_W = hasFlow ? Math.floor((LAYOUT.width - LAYOUT.gap) / 2) : LAYOUT.width;
  const FLOW_W = hasFlow ? LAYOUT.width - LAYOUT.gap - PIE_W : LAYOUT.width;
  const FIG_H = 500;

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

  // --- side-by-side figures -----------------------------------------------------------
  const row = document.createElement("div");
  row.className = "panel-row";
  // Each plot gets its own bordered card, so two charts standing side by side read as two
  // plots rather than one sprawling graphic. The card is the box; the slot inside it is what
  // the drawing function is handed, so neither chart has to know it is in a card.
  // Each card is a <figure class="chart">, which is the unit js/figure_export.js offers a PNG
  // for. That is the whole reason they are figures and not plain boxes: the pie and the flow
  // chart are separately useful, so each gets its own download rather than sharing one
  // capture of the pair.
  const pieCard = document.createElement("figure");
  pieCard.className = "chart plot-card";
  pieCard.setAttribute("data-export-name", "position_one_step_pie");
  const flowCard = document.createElement("figure");
  flowCard.className = "chart plot-card";
  flowCard.setAttribute("data-export-name", "position_supply_flow");
  const pieSlot = document.createElement("div");
  const flowSlot = document.createElement("div");
  // The change selector belongs to the pie, so it lives in the pie's card rather than
  // floating above both figures as a page-level control.
  ctrl.classList.add("card-controls");
  pieCard.appendChild(ctrl);
  pieCard.appendChild(pieSlot);
  flowCard.appendChild(flowSlot);
  // Explicit bases: a flex item holding an SVG at width:100% has no intrinsic width and
  // would otherwise collapse, which is what made these figures render far smaller than asked.
  pieCard.style.flex = `1 1 ${PIE_W}px`;
  pieCard.style.maxWidth = `${PIE_W}px`;
  if (hasFlow) {
    flowCard.style.flex = `1 1 ${FLOW_W}px`;
    flowCard.style.maxWidth = `${FLOW_W}px`;
  }
  row.append(pieCard);
  if (hasFlow) row.append(flowCard);
  chartEl.appendChild(row);
  // drawSupplyFlow writes its explanation here; it is never added to the page, only read
  // back into the flow card's Details dialog.
  const flowCaption = document.createElement("div");

  function pick(key) {
    if (!byKey.has(key)) return;
    selectedKey = key;
    sel.value = key;
    redraw();
  }

  function redraw() {
    pieChart(pieSlot, oneStep.map((r) => ({
      label: keyOf(r), value: r.mu_aa, highlight: false,
      // Coloured by the destination residue, matching the flow chart beside it. A synonymous
      // change keeps the residue, so it takes that residue's colour too.
      colour: residueColour(r.alt_aa),
    })), { // Short title: the position is in the heading above and the caption below, so the
           // chart only has to say which of the two figures it is. The compact legend drops
           // the rate, which the selector above and the table below both already give, and
           // the width it frees goes into the circle.
           title: hasFlow ? `${label}: from the NL4-3 codon`
                          : `Single-nucleotide changes at ${label}`,
           width: opts.width || PIE_W, height: opts.height || FIG_H,
           legendWidth: hasFlow ? 132 : undefined, compactLegend: hasFlow,
           selected: selectedKey, onSelect: pick });

    if (hasFlow) {
      // No title: the figure has a heading and a caption, and the chart's own header block
      // would be the first thing to collide at this width. Height is left to the chart,
      // which sizes itself from the number of blocks it has to stack.
      drawSupplyFlow(flowSlot, flowCaption, opts.flow.stat, opts.flow.matrix,
                     { weight: "supply", width: FLOW_W, height: FIG_H,
                       title: `${label}: across circulating virus`, subtitle: false,
                       // Pooled harder than the batch tab's 2%, so the block count fits the
                       // shared box height with its labels still legible. The threshold is
                       // named in the chart's own caption, which the Details dialog carries.
                       minCodonFreq: 0.05,
                       refCodon: opts.refCodon || null });
    }
  }

  sel.addEventListener("change", () => pick(sel.value));
  redraw();

  // The explanation now lives behind each card's Details button rather than as prose under
  // the row. The visible caption is one line naming the position, so the figure still says
  // what it is at a glance and a downloaded PNG is not anonymous.
  // Each card now carries its own title inside the plot, so the shared caption under the
  // pair would only say it a third time.
  captionEl.innerHTML = hasFlow
    ? ""
    : `<strong>${label}</strong> \u2014 where a single substitution's supply goes.`;

  setDetails(pieCard,
    "<p>Each slice is an amino-acid change reachable from " +
    (opts.custom ? "<strong>your</strong> codon" : "the reference codon") +
    " in one nucleotide substitution, sized by its share of the one-step supply at this " +
    "position. <code>*</code> marks a synonymous change. Click a slice, or use the " +
    "selector above the chart.</p><p>" +
    (opts.custom
      ? "Rates use the codon in your sequence, so they differ from the NL4-3 atlas wherever "
        + "that codon differs. "
      : "Rates are read from the published NL4-3 atlas. ") +
    "Slices are coloured by the destination residue" +
    (hasFlow ? ", the same key the chart beside it uses" : "") +
    ".</p><p><strong>This is generation, not establishment</strong> \u2014 a variant " +
    "generated here can still be lost to drift, purged by selection, or never sampled.</p>",
    { title: `One-step changes from the codon at ${label}` });

  if (hasFlow) {
    setDetails(flowCard, flowDetailHtml(label, flowCaption.innerHTML),
               { title: `${label} across circulating virus` });
  }

  return true;
}
