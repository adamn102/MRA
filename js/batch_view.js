// SPDX-License-Identifier: MIT
//
// batch_view.js -- the two figures on the batch-analysis tab.
//
// WHY MIN-MAX AND NOT AN INTERQUARTILE BAND. Measured on a 24-sequence batch: the IQR is
// exactly zero at all 8,213 sites. Per-site rates differ only where codon identity or local
// context differs, so at any one site at most a couple of samples deviate and the middle 50%
// always agrees -- an IQR ribbon renders as a flat invisible line. The min-max range is
// non-zero at 1,922 sites and reaches ~19x the median, so that is the envelope drawn here.
//
// For the same reason the batch median alone is close to any single sequence's track and says
// little about the batch. The quantity that is genuinely batch-level is the per-codon
// divergence, log2(max/min), which is why it gets its own panel rather than a tooltip.

import { genomeTrack, dotPlot, barChart, LAYOUT } from "./chart.js";
import { smoothSeries } from "./explore_core.js";
import { codonGroupsAtSite } from "./batch_core.js";
import { formatRate } from "./table.js";

/**
 * Batch supply track: median line, min-max envelope, gene band, clickable.
 *
 * @param {object} o  {stats, segments, first, last, win, marker, markerLabel, onSelect,
 *                     nSamples, statsEl}
 */
/**
 * Amino-acid mutability along the genome: one point per CODON, not per nucleotide.
 *
 * The y value is the summed rate of every one-step change that alters the residue. That is
 * the quantity with evolutionary meaning -- a nucleotide rate says nothing about which amino
 * acid is reached -- and it is what the rest of this page is keyed on.
 *
 * @param {object} o {stats, segments, first, last, win, marker, markerLabel, onSelect,
 *                    nSamples, field}
 */
export function drawMutabilityTrack(container, o) {
  container.innerHTML = "";
  const stats = o.stats;
  if (!stats.length) return false;
  const xs = stats.map((r) => r.nt_start);
  const asMap = (key) => new Map(stats.map((r) => [r.nt_start, { total: r[key] }]));
  const win = o.win || 51;
  const med = smoothSeries(xs, asMap("median"), "total", win);
  const lo = smoothSeries(xs, asMap("min"), "total", win);
  const hi = smoothSeries(xs, asMap("max"), "total", win);
  const ribbon = med.map((p, i) => [p[0], lo[i] ? lo[i][1] : p[1], hi[i] ? hi[i][1] : p[1]]);

  genomeTrack(container, [{ points: med, label: "median", colour: "#0b5c8a", width: 1.2 }], {
    geneBand: o.segments, xMin: o.first, xMax: o.last,
    marker: o.marker, markerLabel: o.markerLabel,
    ribbon, ribbonColour: "#0b5c8a", ribbonOpacity: 0.2, ribbonLabel: "min\u2013max",
    title: `Amino-acid mutability across ${o.nSamples} sequences ` +
           `(${win}-codon smoothing)`,
    yLabel: "one-step rate of residue change",
    onSelect: o.onSelect, width: LAYOUT.width, height: 250,
  });
  return true;
}

export function drawMutabilityDivergenceTrack(container, o) {
  container.innerHTML = "";
  const pts = o.stats
    .filter((r) => r.divergence !== null && r.divergence > 0)
    .map((r) => [r.nt_start, r.divergence]);
  if (!pts.length) {
    const p = document.createElement("p");
    p.className = "caveat";
    p.textContent = "Every sample carries the same codon at every shared position, so the " +
                    "mutability is identical across the batch and there is nothing to plot.";
    container.appendChild(p);
    return false;
  }
  // Reuse genomeTrack's frame by passing the divergence as a series on a linear scale: the
  // track's own y-scale is logarithmic, so stems are drawn here instead.
  const width = LAYOUT.width, height = 150;
  const padLeft = LAYOUT.padLeft, padRight = LAYOUT.padRight, padTop = 18, padBottom = 34;
  const plotW = width - padLeft - padRight, plotH = height - padTop - padBottom;
  const x0 = o.first, x1 = o.last;
  const yMax = Math.max(...pts.map((p) => p[1])) * 1.08;
  const NS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("width", "100%");
  svg.setAttribute("role", "img");
  const mk = (n, a) => { const e = document.createElementNS(NS, n);
    for (const [k, v] of Object.entries(a)) e.setAttribute(k, v); return e; };
  const sx = (p) => padLeft + ((p - x0) / (x1 - x0)) * plotW;
  const sy = (v) => padTop + plotH - (v / yMax) * plotH;

  svg.appendChild(mk("text", { x: 0, y: 12, "font-size": 12 })).textContent =
    "Where mutability differs between samples";
  for (const frac of [0, 0.5, 1]) {
    const y = padTop + plotH - frac * plotH;
    svg.appendChild(mk("line", { x1: padLeft, x2: padLeft + plotW, y1: y, y2: y,
                                 stroke: "#e3e8ec", "stroke-width": 0.8 }));
    const t = mk("text", { x: padLeft - 6, y: y + 3.5, "text-anchor": "end",
                           "font-size": 10, fill: "#5c6670" });
    t.textContent = (frac * yMax).toFixed(1);
    svg.appendChild(t);
  }
  for (const [pos, v] of pts) {
    svg.appendChild(mk("line", { x1: sx(pos), x2: sx(pos), y1: sy(0), y2: sy(v),
                                 stroke: "#b5485d", "stroke-width": 0.6, "stroke-opacity": 0.8 }));
  }
  if (Number.isFinite(o.marker)) {
    svg.appendChild(mk("line", { x1: sx(o.marker), x2: sx(o.marker), y1: padTop,
                                 y2: padTop + plotH, stroke: "#b5485d", "stroke-width": 1.2,
                                 "stroke-dasharray": "4 3" }));
  }
  svg.appendChild(mk("line", { x1: padLeft, x2: padLeft + plotW, y1: sy(0), y2: sy(0),
                               stroke: "#999" }));
  const yl = mk("text", { x: 11, y: padTop + plotH / 2, "font-size": 10, fill: "#5c6670",
                          "text-anchor": "middle",
                          transform: `rotate(-90 11 ${padTop + plotH / 2})` });
  yl.textContent = "log2(max/min)";
  svg.appendChild(yl);
  if (typeof o.onSelect === "function") {
    const hit = mk("rect", { x: padLeft, y: padTop, width: plotW, height: plotH,
                             fill: "transparent", cursor: "crosshair" });
    hit.addEventListener("click", (ev) => {
      const box = svg.getBoundingClientRect();
      const scale = box.width / width;
      const px = (ev.clientX - box.left) / scale;
      const pos = Math.round(x0 + ((px - padLeft) / plotW) * (x1 - x0));
      o.onSelect(Math.min(Math.max(pos, x0), x1));
    });
    svg.appendChild(hit);
  }
  container.appendChild(svg);
  return true;
}

/**
 * The site-specific comparison: every sample's supply at one position, ordered by value, with
 * the codon each sample carries. Samples that do not cover the site are listed as absent rather
 * than omitted -- an omission would read as "no difference here".
 */
/**
 * The site-specific comparison, as one row per distinct codon rather than per sample.
 *
 * The earlier per-sample layout was unreadable with real data: every sample sharing a codon has
 * an identical rate, so the rows were mostly duplicates, and real sequence names are long
 * enough to crowd out the plot. Codons are short, few, and the thing that actually differs.
 * Sample names are listed in the table beneath, where length costs nothing.
 */
export function drawSiteComparison(container, captionEl, rows, refPosition, extra = {}) {
  container.innerHTML = "";
  const covered = rows.filter((r) => r.covered);
  const absent = rows.length - covered.length;
  if (!covered.length) {
    captionEl.textContent =
      `No sample in this batch covers ${extra.label || "NL4-3 position " + refPosition}.`;
    return false;
  }
  const groups = codonGroupsAtSite(covered, extra.refCodon || null);
  const vals = covered.map((r) => r.value).sort((a, b) => a - b);
  const mid = vals.length % 2 ? vals[vals.length >> 1]
            : 0.5 * (vals[(vals.length >> 1) - 1] + vals[vals.length >> 1]);

  // Row label carries the codon, its sample count, and -- only when the group is small enough
  // to name without crowding -- the samples themselves. Everything else is in the table.
  const NAME_LIMIT = 2;
  const plotRows = groups.map((g) => {
    const named = g.n <= NAME_LIMIT ? g.samples.join(", ") : `${g.n} sequences`;
    return {
      label: `${g.codon}  (${named})`,
      value: (mid > 0 && g.rate > 0) ? Math.log2(g.rate / mid) : 0,
      note: `${g.rate.toExponential(2)}${g.is_reference ? "  \u25c0 NL4-3" : ""}`,
    };
  });
  dotPlot(container, plotRows, {
    title: `Supply at ${extra.label || "NL4-3 nt " + refPosition.toLocaleString()} by codon`,
    width: LAYOUT.width,
    xLabel: "log2 of each codon's rate relative to the batch median at this site",
  });

  const max = groups[0].rate, min = groups[groups.length - 1].rate;
  const nonUniform = groups.filter((g) => !g.uniform);
  captionEl.innerHTML =
    `${covered.length} of ${rows.length} samples cover this site` +
    (absent ? `; ${absent} do not and are not plotted` : "") +
    `, carrying <strong>${groups.length} distinct codon${groups.length === 1 ? "" : "s"}</strong>. ` +
    `Batch median ${formatRate(mid)}. ` +
    (groups.length > 1 && min > 0
      ? `The fastest codon generates substitutions <strong>${(max / min).toFixed(1)}\u00d7</strong> ` +
        `more often than the slowest (${formatRate(min)} to ${formatRate(max)}). Every sample ` +
        `with the same codon has the same rate: the model is identical for all of them, so the ` +
        `codon is what differs. `
      : `Every covering sample carries the same codon, so they share one rate. `) +
    (nonUniform.length
      ? `<strong>${nonUniform.length} codon group spans a range</strong> because the base just ` +
        `after the codon differs between those samples, which changes the 3\u2032 context of its ` +
        `last position. `
      : "") +
    `Sample names are in the table below.`;
  return true;
}

/**
 * Observed residue diversity against predicted mutability -- the figure that separates
 * mutational supply from selection.
 *
 * Each point is a codon. The x axis is what the model predicts about how readily the residue
 * changes; the y axis is what the submitted sequences actually show. The informative points
 * are the ones off the diagonal:
 *
 *   high mutability, no diversity -> variants are generated often and are not there. Something
 *                                    removes them, which is what purifying selection looks
 *                                    like from here: the change is mutationally easy and
 *                                    evidently not tolerated.
 *   low mutability, high diversity -> variation that is hard to generate yet common, pointing
 *                                    to positive selection or to shared ancestry rather than
 *                                    to repeated independent mutation.
 *
 * The caption states the limit rather than leaving it implied: submitted sequences are related
 * by descent, so diversity is not an independent count of mutation-and-selection events.
 */
export function drawDiversityVsMutability(container, captionEl, points, opts = {}) {
  container.innerHTML = "";
  if (!points.length) {
    captionEl.textContent = "No position is covered by enough samples to compare.";
    return false;
  }
  const width = opts.width || LAYOUT.width, height = 420;
  const padLeft = 72, padRight = 150, padTop = 26, padBottom = 48;
  const plotW = width - padLeft - padRight, plotH = height - padTop - padBottom;
  const NS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("width", "100%");
  svg.setAttribute("role", "img");
  const mk = (n, a, t) => {
    const e = document.createElementNS(NS, n);
    for (const [k, v] of Object.entries(a)) e.setAttribute(k, v);
    if (t != null) e.textContent = t;
    return e;
  };
  const xs = points.map((p) => Math.log10(p.mutability));
  const x0 = Math.floor(Math.min(...xs)), x1 = Math.ceil(Math.max(...xs));
  const yMax = Math.max(0.05, ...points.map((p) => p.diversity));
  const sx = (v) => padLeft + ((v - x0) / (x1 - x0 || 1)) * plotW;
  const sy = (v) => padTop + plotH - (v / yMax) * plotH;

  svg.appendChild(mk("text", { x: 0, y: 13, "font-size": 12 },
                     opts.title || "Observed diversity against predicted mutability"));
  for (let ex = x0; ex <= x1; ex++) {
    svg.appendChild(mk("line", { x1: sx(ex), x2: sx(ex), y1: padTop, y2: padTop + plotH,
                                 stroke: "#eef2f4", "stroke-width": 0.8 }));
    svg.appendChild(mk("text", { x: sx(ex), y: padTop + plotH + 15, "text-anchor": "middle",
                                 "font-size": 10, fill: "#5c6670" }, "10" + supDigits(ex)));
  }
  for (const f of [0, 0.25, 0.5, 0.75, 1]) {
    const v = f * yMax;
    svg.appendChild(mk("line", { x1: padLeft, x2: padLeft + plotW, y1: sy(v), y2: sy(v),
                                 stroke: "#eef2f4", "stroke-width": 0.8 }));
    svg.appendChild(mk("text", { x: padLeft - 6, y: sy(v) + 3.5, "text-anchor": "end",
                                 "font-size": 10, fill: "#5c6670" }, v.toFixed(2)));
  }
  // The batch's own median is the divider, so "high" and "low" are read off the submitted set
  // rather than from an arbitrary absolute threshold.
  const sorted = xs.slice().sort((a, b) => a - b);
  const medX = sorted[sorted.length >> 1];
  svg.appendChild(mk("line", { x1: sx(medX), x2: sx(medX), y1: padTop, y2: padTop + plotH,
                               stroke: "#9aa5ad", "stroke-width": 1, "stroke-dasharray": "4 3" }));
  svg.appendChild(mk("text", { x: sx(medX) + 4, y: padTop + 10, "font-size": 8.5,
                               fill: "#9aa5ad" }, "median mutability"));

  // Invariant positions are the majority and would hide everything else if drawn opaque.
  for (const p of points) {
    const invariant = p.diversity === 0;
    svg.appendChild(mk("circle", {
      cx: sx(Math.log10(p.mutability)).toFixed(1), cy: sy(p.diversity).toFixed(1),
      r: invariant ? 1.5 : 2.6,
      fill: invariant ? "#b8c2c9" : "#0b5c8a",
      "fill-opacity": invariant ? 0.35 : 0.7 }));
  }
  const constrained = points
    .filter((p) => p.diversity === 0 && Math.log10(p.mutability) > medX)
    .sort((a, b) => b.mutability - a.mutability);
  for (const p of constrained.slice(0, opts.labelCount || 6)) {
    svg.appendChild(mk("circle", { cx: sx(Math.log10(p.mutability)).toFixed(1),
                                   cy: sy(0).toFixed(1), r: 3.2, fill: "none",
                                   stroke: "#b5485d", "stroke-width": 1.2 }));
  }
  let ly = padTop + 14;
  svg.appendChild(mk("text", { x: padLeft + plotW + 8, y: ly, "font-size": 9, fill: "#b5485d" },
                     "highest supply, no variation:"));
  for (const p of constrained.slice(0, opts.labelCount || 6)) {
    ly += 12;
    svg.appendChild(mk("text", { x: padLeft + plotW + 8, y: ly, "font-size": 9, fill: "#5c6670" },
                       opts.labelFor ? opts.labelFor(p) : p.protein + " " + p.position));
  }
  svg.appendChild(mk("line", { x1: padLeft, x2: padLeft + plotW, y1: padTop + plotH,
                               y2: padTop + plotH, stroke: "#999" }));
  svg.appendChild(mk("text", { x: padLeft + plotW / 2, y: padTop + plotH + 34,
                               "text-anchor": "middle", "font-size": 10, fill: "#5c6670" },
                     "predicted one-step rate of residue change"));
  svg.appendChild(mk("text", { x: 14, y: padTop + plotH / 2, "font-size": 10, fill: "#5c6670",
                               "text-anchor": "middle",
                               transform: `rotate(-90 14 ${padTop + plotH / 2})` },
                     "observed diversity (1 \u2212 majority residue frequency)"));
  container.appendChild(svg);

  const nVar = points.filter((p) => p.diversity > 0).length;
  captionEl.innerHTML =
    `${points.length.toLocaleString()} codons. <strong>${nVar.toLocaleString()}</strong> vary in ` +
    `residue across the batch; <strong>${constrained.length.toLocaleString()}</strong> sit above ` +
    `the median mutability yet show no variation at all \u2014 changes that are readily generated ` +
    `and nonetheless absent, which is what purifying selection looks like from here. ` +
    (opts.rho !== null && opts.rho !== undefined
      ? `Spearman &rho; between mutability and diversity is <strong>${opts.rho.toFixed(3)}</strong>. `
      : "") +
    `<strong>One caution that changes how this reads:</strong> the submitted sequences are ` +
    `related by descent, so diversity counts shared ancestry as well as independent mutation, ` +
    `and a consensus sequence has already discarded within-host variation. Read this as a ` +
    `description of the set you submitted, not as an estimate of selection.`;
  return true;
}

/** Unicode superscript for a base-10 exponent, for axis labels without MathML. */
function supDigits(e) {
  const map = { "-": "\u207b", 0: "\u2070", 1: "\u00b9", 2: "\u00b2", 3: "\u00b3", 4: "\u2074",
                5: "\u2075", 6: "\u2076", 7: "\u2077", 8: "\u2078", 9: "\u2079" };
  return String(e).split("").map((c) => map[c] || c).join("");
}

/** A horizontal stacked-proportion bar: what the batch carries at one position. */
function compositionBar(svg, mk, x, y, w, h, entries, total, palette, labelFmt) {
  let cx = x;
  entries.forEach(([label, count], i) => {
    const frac = count / total;
    const bw = frac * w;
    svg.appendChild(mk("rect", { x: cx.toFixed(1), y, width: Math.max(bw, 0.6).toFixed(1),
                                 height: h, fill: palette[i % palette.length],
                                 stroke: "#fff", "stroke-width": 0.6 }));
    // Only label a slice wide enough to hold the text; the rest are in the legend.
    if (bw > 28) {
      svg.appendChild(mk("text", { x: (cx + bw / 2).toFixed(1), y: y + h / 2 + 3.5,
                                   "text-anchor": "middle", "font-size": 9.5, fill: "#fff" },
                         labelFmt(label, frac)));
    }
    cx += bw;
  });
}

const COMP_PALETTE = ["#0b5c8a", "#2e8b9a", "#6aaa64", "#c2973a", "#c2561a", "#8a4f7d",
                      "#5c6670", "#b5485d"];

/**
 * What the batch carries at one codon: the residue composition and the codon composition,
 * stacked so the relationship between them is visible.
 *
 * The two bars answer different questions and are deliberately shown together. A position can
 * be completely invariant in residue while its codons differ -- synonymous variation is
 * invisible to the protein but changes the mutational supply, which is the whole reason this
 * page works at codon level rather than at residue level alone.
 */
export function drawPositionComposition(container, captionEl, stat, opts = {}) {
  container.innerHTML = "";
  const width = opts.width || LAYOUT.width, rowH = 30, gap = 38;
  const height = 2 * rowH + gap + 54;
  const NS = "http://www.w3.org/2000/svg";
  const svg = document.createElementNS(NS, "svg");
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  svg.setAttribute("width", "100%");
  svg.setAttribute("role", "img");
  const mk = (n, a, t) => {
    const e = document.createElementNS(NS, n);
    for (const [k, v] of Object.entries(a)) e.setAttribute(k, v);
    if (t != null) e.textContent = t;
    return e;
  };
  const padLeft = 112, padRight = 16;
  const w = width - padLeft - padRight;
  const n = stat.n_covered;

  svg.appendChild(mk("text", { x: 0, y: 13, "font-size": 12 },
                     opts.title || "What the batch carries here"));
  svg.appendChild(mk("text", { x: padLeft - 8, y: 24 + rowH / 2 + 3.5, "text-anchor": "end",
                               "font-size": 10, fill: "#5c6670" }, "amino acid"));
  compositionBar(svg, mk, padLeft, 24, w, rowH, stat.aa_counts, n, COMP_PALETTE,
                 (l, f) => `${l}  ${(100 * f).toFixed(0)}%`);
  const y2 = 24 + rowH + gap;
  svg.appendChild(mk("text", { x: padLeft - 8, y: y2 + rowH / 2 + 3.5, "text-anchor": "end",
                               "font-size": 10, fill: "#5c6670" }, "codon"));
  compositionBar(svg, mk, padLeft, y2, w, rowH, stat.codon_counts, n, COMP_PALETTE,
                 (l, f) => `${l}  ${(100 * f).toFixed(0)}%`);
  // Full legends under each bar, since narrow slices carry no inline label.
  svg.appendChild(mk("text", { x: padLeft, y: 24 + rowH + 14, "font-size": 9, fill: "#5c6670" },
                     stat.aa_counts.map(([a, c]) => `${a} ${c}`).join("   ")));
  svg.appendChild(mk("text", { x: padLeft, y: y2 + rowH + 14, "font-size": 9, fill: "#5c6670" },
                     stat.codon_counts.map(([c, k]) => `${c} ${k}`).join("   ")));
  container.appendChild(svg);

  const synonymousOnly = stat.n_aa === 1 && stat.n_codons > 1;
  captionEl.innerHTML =
    `${n} sample${n === 1 ? "" : "s"} cover this position. ` +
    `<strong>${stat.majority_aa}</strong> in ${(100 * stat.majority_aa_freq).toFixed(0)}% of them` +
    (stat.n_aa > 1 ? `, across ${stat.n_aa} residues (Shannon entropy ` +
                     `${stat.entropy.toFixed(2)} bits)` : ", with no residue variation") +
    `; <strong>${stat.n_codons}</strong> distinct codon${stat.n_codons === 1 ? "" : "s"}. ` +
    (synonymousOnly
      ? `<strong>The variation here is entirely synonymous</strong> \u2014 every sample encodes ` +
        `the same residue, but not with the same codon, so the protein is identical while the ` +
        `mutational supply is not. `
      : "") +
    (stat.divergence !== null && stat.divergence > 0
      ? `Mutability spans ${Math.pow(2, stat.divergence).toFixed(1)}-fold across the batch.`
      : "Every covering sample has the same mutability here.");
  return true;
}

/**
 * Where this position's supply goes: the batch-weighted rate to each destination residue.
 *
 * Weighted by the codon frequencies observed in the batch, so the bar heights are the supply
 * as distributed in the submitted set rather than the supply from one reference codon. The
 * fold range across codons is shown per destination, because a destination whose rate depends
 * strongly on which synonymous codon a sample carries is one where codon usage is itself an
 * evolutionary variable.
 */
export function drawDestinationSupply(container, captionEl, dest, stat, opts = {}) {
  container.innerHTML = "";
  const rows = dest.filter((d) => !d.is_synonymous && d.weighted_rate > 0);
  if (!rows.length) {
    captionEl.textContent = "No residue change is reachable from this codon in one step.";
    return false;
  }
  const shown = rows.slice(0, opts.limit || 12);
  barChart(container, shown.map((d) => ({
    label: d.alt_aa + (d.is_stop ? " (stop)" : ""),
    value: d.weighted_rate,
    note: d.fold_range && d.fold_range > 1.01
      ? `${d.fold_range.toFixed(1)}\u00d7 across codons` : "",
  })), { title: opts.title || "Where a single substitution goes from here",
  // Linear scale deliberately: the bars are read as shares of this position's supply, and
  // the synonymous destination (which can be orders of magnitude larger) is excluded
  // above, so nothing here compresses the rest. Exact rates print per row regardless.
         width: opts.width || LAYOUT.width,
         xLabel: "batch-weighted rate per replication" });

  const total = rows.reduce((t, d) => t + d.weighted_rate, 0);
  const top = shown[0];
  const sensitive = rows.filter((d) => d.fold_range && d.fold_range > 2);
  captionEl.innerHTML =
    `Rates are weighted by how often each codon occurs in this batch, so they describe the ` +
    `supply in the set you submitted rather than in one reference codon. ` +
    `The most accessible change is <strong>${stat.majority_aa}\u2192${top.alt_aa}</strong> at ` +
    `${formatRate(top.weighted_rate)}, ` +
    `${(100 * top.weighted_rate / total).toFixed(0)}% of the non-synonymous supply here. ` +
    (sensitive.length
      ? `<strong>${sensitive.length} destination${sensitive.length === 1 ? "" : "s"} vary more ` +
        `than two-fold</strong> depending on which synonymous codon a sample carries ` +
        `(${sensitive.slice(0, 4).map((d) => d.alt_aa).join(", ")}) \u2014 at those, codon usage ` +
        `changes how readily the change is generated. `
      : "") +
    `Stops are shown but are lethal, not exploratory; synonymous supply is excluded here and ` +
    `reported separately.`;
  return true;
}
