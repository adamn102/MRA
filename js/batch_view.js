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
import { codonGroupsAtSite, AA_ORDER } from "./batch_core.js";
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
  // SMOOTHING. One stem per diverging codon is unreadable once a few hundred positions
  // disagree, so a window mean is drawn over the stems. The mean is taken over EVERY covered
  // codon in the window, counting non-diverging ones as the zero they are -- averaging only
  // the diverging ones would make a region with one big outlier look like a region where
  // everything disagrees. At win = 1 the mean is the raw value and only the stems show.
  const win = Math.max(1, o.win || 1);
  let smooth = null;
  if (win > 1) {
    const xs = o.stats.map((r) => r.nt_start);
    const byPos = new Map(o.stats.map((r) => [r.nt_start,
      { total: r.divergence === null ? 0 : r.divergence }]));
    smooth = smoothSeries(xs, byPos, "total", win);
  }
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
  // Stems fade back when a smoothed line is drawn over them, so the line reads first while
  // individual spikes stay locatable.
  const stemOpacity = smooth ? 0.22 : 0.8;
  for (const [pos, v] of pts) {
    svg.appendChild(mk("line", { x1: sx(pos), x2: sx(pos), y1: sy(0), y2: sy(v),
                                 stroke: "#b5485d", "stroke-width": 0.6,
                                 "stroke-opacity": stemOpacity }));
  }
  if (smooth && smooth.length) {
    const d = smooth.map((p, i) => `${i ? "L" : "M"}${sx(p[0]).toFixed(1)},${sy(p[1]).toFixed(1)}`)
                    .join(" ");
    svg.appendChild(mk("path", { d, fill: "none", stroke: "#8c1d33", "stroke-width": 1.6,
                                 "stroke-linejoin": "round" }));
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
  // `only: "aa"` renders just the residue bar. The flow chart below already carries the codon
  // composition in its left column, and showing it twice is what made this section hard to read.
  const aaOnly = opts.only === "aa";
  const width = opts.width || LAYOUT.width, rowH = 30, gap = 38;
  const height = aaOnly ? rowH + 44 : 2 * rowH + gap + 54;
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
  if (aaOnly) {
    svg.appendChild(mk("text", { x: padLeft, y: 24 + rowH + 14, "font-size": 9, fill: "#5c6670" },
                       stat.aa_counts.map(([a, c]) => `${a} ${c}`).join("   ")));
    container.appendChild(svg);
    captionEl.innerHTML =
      `${n} sample${n === 1 ? "" : "s"} cover this position, carrying ` +
      `<strong>${stat.majority_aa}</strong> in ${(100 * stat.majority_aa_freq).toFixed(0)}%` +
      (stat.n_aa > 1
        ? ` and ${stat.n_aa - 1} other residue${stat.n_aa === 2 ? "" : "s"} ` +
          `(Shannon entropy ${stat.entropy.toFixed(2)} bits)`
        : " and nothing else") +
      `, on <strong>${stat.n_codons}</strong> distinct codon` +
      `${stat.n_codons === 1 ? "" : "s"}.` +
      (stat.n_aa === 1 && stat.n_codons > 1
        ? ` <strong>Every sample encodes the same residue</strong> \u2014 the variation here is ` +
          `entirely synonymous, so the protein is identical while the supply is not.`
        : "");
    return true;
  }
  svg.appendChild(mk("text", { x: padLeft - 8, y: y2 + rowH / 2 + 3.5, "text-anchor": "end",
                               "font-size": 10, fill: "#5c6670" }, "codon"));
  // The residue each codon encodes, so a reader does not have to translate in their head and
  // can see at a glance which codon differences are synonymous.
  const withAa = (c) => stat.codon_aa && stat.codon_aa.get(c)
    ? `${c} (${stat.codon_aa.get(c)})` : c;
  compositionBar(svg, mk, padLeft, y2, w, rowH, stat.codon_counts, n, COMP_PALETTE,
                 (l, f) => `${withAa(l)}  ${(100 * f).toFixed(0)}%`);
  // Full legends under each bar, since narrow slices carry no inline label.
  svg.appendChild(mk("text", { x: padLeft, y: 24 + rowH + 14, "font-size": 9, fill: "#5c6670" },
                     stat.aa_counts.map(([a, c]) => `${a} ${c}`).join("   ")));
  svg.appendChild(mk("text", { x: padLeft, y: y2 + rowH + 14, "font-size": 9, fill: "#5c6670" },
                     stat.codon_counts.map(([c, k]) => `${withAa(c)} ${k}`).join("   ")));
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
 * SYNONYMOUS ROUTES ARE EXCLUDED throughout. They are often the largest rate at a position
 * (AAG->AAA is a G->A transition at ~1e-3) and would swamp every residue-changing band while
 * saying nothing about where the protein can go. The per-position table carries them.
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
      if (alt === own || !(d[i] > 0)) continue;
      m.set(alt, d[i]);
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
  const tail = order.slice(limit);
  order = order.slice(0, limit);
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

  const width = opts.width || LAYOUT.width;
  const padTop = 34, padBottom = 14;
  const height = opts.height || Math.max(260, 42 * Math.max(live.length, order.length) + 70);
  const plotH = height - padTop - padBottom;
  const BW = 16, xL = 320, xR = width - 290;
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

  svg.appendChild(mk("text", { x: 0, y: 13, "font-size": 12 },
                     opts.title || "Where residue change at this position comes from, " +
                                   "and where it goes"));
  svg.appendChild(mk("text", { x: 0, y: 27, "font-size": 9.5, fill: "#5c6670" },
                     mode === "supply"
                       ? "Sized by mutational supply: a block is its share of the change this "
                         + "position generates."
                       : "Sized by how common each codon is: a block is its share of the "
                         + "sequences."));
  svg.appendChild(mk("text", { x: xL + BW, y: padTop - 5, "text-anchor": "end", "font-size": 9,
                               fill: "#5c6670" }, "codon carried"));
  svg.appendChild(mk("text", { x: xR - BW, y: padTop - 5, "font-size": 9, fill: "#5c6670" },
                     "residue it becomes"));

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
      svg.appendChild(mk("path", { d, fill: FLOW_PALETTE[ci % FLOW_PALETTE.length],
                                   "fill-opacity": 0.42 }));
    }
  });

  live.forEach((c, ci) => {
    const [y0, y1] = cy.get(c);
    svg.appendChild(mk("rect", { x: xL, y: y0.toFixed(1), width: BW,
                                 height: (y1 - y0).toFixed(1),
                                 fill: FLOW_PALETTE[ci % FLOW_PALETTE.length],
                                 "fill-opacity": 0.9 }));
    const mid = (y0 + y1) / 2;
    const isPool = c === OTHER_CODONS;
    svg.appendChild(mk("text", { x: xL - 10, y: mid - 2, "text-anchor": "end", "font-size": 11,
                                 "font-family": isPool ? "inherit" : "monospace" },
                       isPool ? `${pooledCodons} rarer codons`
                              : `${c} (${stat.codon_aa.get(c)})` +
                                (c === opts.refCodon ? "  \u25c0 NL4-3" : "")));
    // Both numbers on every block, in both modes: their gap is the thing worth seeing.
    svg.appendChild(mk("text", { x: xL - 10, y: mid + 11, "text-anchor": "end", "font-size": 9,
                                 fill: "#5c6670" },
                       `${(100 * freqOf.get(c)).toFixed(0)}% of sequences \u2192 ` +
                       `${(100 * supplyOf.get(c) / supplyTotal).toFixed(0)}% of supply`));
  });
  for (const a of order) {
    const [y0, y1] = dy.get(a);
    svg.appendChild(mk("rect", { x: xR - BW, y: y0.toFixed(1), width: BW,
                                 height: (y1 - y0).toFixed(1), fill: "#5c6670",
                                 "fill-opacity": 0.55 }));
    const label = a === OTHER ? "other" : a === "*" ? "stop" : a;
    svg.appendChild(mk("text", { x: xR + 10, y: (y0 + y1) / 2 + 4, "font-size": 11 },
                       `${label}   ${(100 * dTotal.get(a) / total).toFixed(0)}%`));
  }
  container.appendChild(svg);

  let lead = null;
  for (const c of live) {
    if (c === OTHER_CODONS) continue;      // a pooled block is not a finding about one codon
    const f = freqOf.get(c), share = supplyOf.get(c) / supplyTotal;
    if (!lead || share - f > lead.gap) lead = { codon: c, f, share, gap: share - f };
  }
  const topDest = order[0];
  const dName = (a) => a === OTHER ? "the pooled remainder" : a === "*" ? "a stop codon" : a;
  captionEl.innerHTML =
    (mode === "supply"
      ? `Band width is a codon's frequency multiplied by its rate to that residue, so a left ` +
        `block is its share of the <em>supply</em> rather than of the sequences. `
      : `Band width is how common the codon is, split by where that codon's single ` +
        `substitutions lead. Left blocks are shares of <em>sequences</em>; switch the ` +
        `weighting to size them by supply instead. `) +
    (lead && lead.gap > 0.05
      ? `<strong>${lead.codon}</strong> is carried by ${(100 * lead.f).toFixed(0)}% of these ` +
        `sequences but accounts for <strong>${(100 * lead.share).toFixed(0)}%</strong> of what ` +
        `this position generates. `
      : `Supply is spread across the codons roughly in proportion to how common they are. `) +
    `Most goes to <strong>${dName(topDest)}</strong> ` +
    `(${(100 * dTotal.get(topDest) / total).toFixed(0)}%). ` +
    `Synonymous routes are excluded: they are often the largest rate at a position and would ` +
    `swamp the rest without saying anything about where the protein can go.` +
    (tail.length ? ` The ${tail.length} smallest destinations are pooled as "other".` : "") +
    (pooledCodons
      ? ` ${pooledCodons} codons each under ${(100 * minFreq).toFixed(0)}% of the batch are ` +
        `pooled into one block (${(100 * pooledFreq).toFixed(0)}% of sequences between them); ` +
        `the per-position table lists them separately.`
      : "");
  return true;
}

const OTHER = "\u2014other\u2014";
const OTHER_CODONS = "\u2014other codons\u2014";
const FLOW_PALETTE = ["#0b5c8a", "#2e8b9a", "#c2973a", "#c2561a", "#6aaa64", "#8a4f7d",
                      "#5c6670", "#b5485d"];
