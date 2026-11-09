// SPDX-License-Identifier: MIT
//
// chart.js -- the three charts this portal needs, in plain SVG. No dependencies.
//
//   barChart     relative mutational supply across the alternatives at one position
//   genomeTrack  rate along the genome, optionally with a second series overlaid
//   dotPlot      log2 fold change per amino-acid change, for the comparison view
//
// SVG rather than canvas so the output is inspectable, scales on a high-DPI screen,
// and can be saved by right-click. Rates span more than three orders of magnitude, so
// the value axis is log10 wherever a rate is plotted -- and bars are never drawn on a
// log axis, because a bar's length would then encode a ratio to an arbitrary floor.

const NS = "http://www.w3.org/2000/svg";

function el(name, attrs = {}, text = null) {
  const node = document.createElementNS(NS, name);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  if (text !== null) node.textContent = text;
  return node;
}

function svgRoot(width, height, title) {
  const svg = el("svg", { viewBox: `0 0 ${width} ${height}`, width: "100%",
                          role: "img", "aria-label": title || "chart" });
  svg.style.maxWidth = width + "px";
  svg.style.font = "11px -apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif";
  return svg;
}

function niceLogTicks(lo, hi) {
  const ticks = [];
  for (let e = Math.floor(lo); e <= Math.ceil(hi); e++) ticks.push(e);
  return ticks;
}

const supLabel = (e) => "10" + String(e).replace(/-/g, "\u2212").split("").map(
  (c) => ({ "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3", "4": "\u2074",
            "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079",
            "\u2212": "\u207b" }[c] || c)).join("");

// ---------------------------------------------------------------------------------
/**
 * Horizontal bars of mutational supply across alternatives at one position.
 * Bars are drawn on a LINEAR axis of rate relative to the largest alternative, which
 * is what "relative supply" means and keeps bar length honest. The absolute rate and
 * its interval are printed beside each bar.
 *
 * rows: [{label, value, lo, hi, highlight?, note?}]
 */
// ---------------------------------------------------------------------------------
/**
 * One layout width for every full-width figure, so the genome track and the codon panel
 * below it line up instead of drifting apart whenever one of them is resized.
 *
 * 1140 is the content box of `main` (max-width 1180 less 20px padding each side). The SVGs
 * carry width:100% with this as their max-width, so they fill it on a wide screen and scale
 * down on a narrow one.
 */
export const LAYOUT = {
  width: 1140,          // full figure width, shared by the track and the codon panel
  padLeft: 62,          // track plot-area insets, needed to align the position slider
  padRight: 118,
  gap: 16,              // gap between the pie and the supply-flow chart beside it
};

/** Percentage offsets that place a slider exactly under the track's plot area. */
export function trackSliderGeometry(width = LAYOUT.width) {
  const plotW = width - LAYOUT.padLeft - LAYOUT.padRight;
  return { marginLeft: `${(100 * LAYOUT.padLeft / width).toFixed(4)}%`,
           width: `${(100 * plotW / width).toFixed(4)}%`,
           maxWidth: `${width}px` };
}

export function barChart(container, rows, options = {}) {
  container.innerHTML = "";
  if (!rows.length) {
    container.appendChild(Object.assign(document.createElement("p"),
      { className: "caveat", textContent: options.emptyMessage || "Nothing to plot." }));
    return;
  }
  const rowH = 17, padTop = 26, padBottom = 34, padLeft = 54, padRight = 128;
  const width = options.width || 720;
  const height = padTop + rows.length * rowH + padBottom;
  const plotW = width - padLeft - padRight;
  const max = Math.max(...rows.map((r) => r.value || 0));
  // LOG OPTION. Rates at one position can span two orders of magnitude, where a linear bar
  // for the smallest is a sliver against the largest and the row is unreadable. On a log
  // scale the bar length is log10(value) measured up from a floor two decades below the
  // largest value, so every non-zero row has visible length and the ratios stay legible.
  const logScale = !!options.log && max > 0;
  const logFloor = logScale ? Math.log10(max) - (options.logDecades || 2.5) : 0;
  const frac10 = (v) => {
    if (!(v > 0)) return 0;
    const f = (Math.log10(v) - logFloor) / (Math.log10(max) - logFloor);
    return Math.max(0, Math.min(1, f));
  };
  const svg = svgRoot(width, height, options.title);

  svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title || ""));

  rows.forEach((r, i) => {
    const y = padTop + i * rowH;
    svg.appendChild(el("text", { x: padLeft - 6, y: y + 11, "text-anchor": "end",
                                 "font-weight": r.highlight ? 600 : 400 }, r.label));
    const frac = logScale ? frac10(r.value) : (max > 0 ? (r.value || 0) / max : 0);
    const w = Math.max(frac * plotW, r.value > 0 ? 1.5 : 0);
    svg.appendChild(el("rect", { x: padLeft, y: y + 3, width: w, height: rowH - 7,
                                 fill: r.highlight ? "#c2561a" : "#0b5c8a",
                                 "fill-opacity": r.highlight ? 0.85 : 0.6 }));
    // Interval as a whisker on the same linear scale.
    if (r.lo != null && r.hi != null && max > 0) {
      const x1 = padLeft + (logScale ? frac10(r.lo) : r.lo / max) * plotW;
      const x2 = padLeft + (logScale ? frac10(r.hi) : r.hi / max) * plotW;
      const yc = y + rowH / 2 - 0.5;
      svg.appendChild(el("line", { x1, x2, y1: yc, y2: yc, stroke: "#333", "stroke-width": 0.9 }));
      for (const x of [x1, x2]) {
        svg.appendChild(el("line", { x1: x, x2: x, y1: yc - 3, y2: yc + 3,
                                     stroke: "#333", "stroke-width": 0.9 }));
      }
    }
    const txt = (r.value === 0 || r.value == null) ? "\u2014" : r.value.toExponential(2);
    svg.appendChild(el("text", { x: width - padRight + 6, y: y + 11, "font-size": 10,
                                 fill: "#333" }, txt + (r.note ? "  " + r.note : "")));
  });

  const axisY = padTop + rows.length * rowH + 4;
  svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: axisY, y2: axisY,
                               stroke: "#999", "stroke-width": 0.8 }));
  for (const frac of [0, 0.25, 0.5, 0.75, 1]) {
    const x = padLeft + frac * plotW;
    svg.appendChild(el("line", { x1: x, x2: x, y1: axisY, y2: axisY + 4, stroke: "#999" }));
    const tickVal = logScale
      ? Math.pow(10, logFloor + frac * (Math.log10(max) - logFloor))
      : frac * max;
    svg.appendChild(el("text", { x, y: axisY + 15, "text-anchor": "middle", "font-size": 10,
                                 fill: "#5c6670" },
                       (!logScale && frac === 0) ? "0" : tickVal.toExponential(1)));
  }
  svg.appendChild(el("text", { x: padLeft + plotW / 2, y: axisY + 29, "text-anchor": "middle",
                               "font-size": 10, fill: "#5c6670" },
                     options.xLabel || "mutational supply per site per replication"));
  container.appendChild(svg);
}

// ---------------------------------------------------------------------------------
/**
 * Rate along the genome. series: [{label, colour, points: [[position, rate], ...]}].
 * Log10 value axis. Optional shaded spans for gene boundaries.
 */
export function genomeTrack(container, series, options = {}) {
  container.innerHTML = "";
  const width = options.width || LAYOUT.width, height = options.height || 260;
  const padLeft = LAYOUT.padLeft, padRight = LAYOUT.padRight, padBottom = 40;
  const x0 = options.xMin || 1, x1 = options.xMax || 8213;

  // An explicit gene band (labelled boxes above the plot) is opt-in via `geneBand`, so the
  // older faint `spans` shading used elsewhere keeps its appearance. Overlapping reading
  // frames -- gag/pol, and tat/rev/vpu/env -- are packed onto as many rows as they need,
  // because drawing them on one row would either hide genes or imply they do not overlap.
  const band = (options.geneBand || []).slice().sort((a, b) => a.start - b.start);
  const bandRows = [];
  for (const g of band) {
    let r = 0;
    while (bandRows[r] && bandRows[r].some((o) => !(g.start > o.end || g.end < o.start))) r++;
    (bandRows[r] = bandRows[r] || []).push(g);
    g._row = r;
  }
  const BOX_H = 12, ROW_GAP = 2;
  const bandH = bandRows.length ? bandRows.length * (BOX_H + ROW_GAP) + 6 : 0;
  const padTop = 24 + bandH;
  const plotW = width - padLeft - padRight, plotH = height - padTop - padBottom;

  // The ribbon's extent counts toward the y-domain: an envelope drawn outside the computed
  // range would be silently clipped at the plot edge and read as a narrower spread.
  const ribbon = options.ribbon || null;
  const ribbonVals = ribbon ? ribbon.flatMap((r) => [r[1], r[2]]).filter((v) => v > 0) : [];
  const all = series.flatMap((s) => s.points.map((p) => p[1])).filter((v) => v > 0)
    .concat(ribbonVals);
  if (!all.length) {
    container.appendChild(Object.assign(document.createElement("p"),
      { className: "caveat", textContent: "No rates to plot for this selection." }));
    return;
  }
  const lo = Math.log10(Math.min(...all)), hi = Math.log10(Math.max(...all));
  const yLo = Math.floor(lo), yHi = Math.ceil(hi);
  const sx = (p) => padLeft + ((p - x0) / (x1 - x0)) * plotW;
  const sy = (v) => padTop + plotH - ((Math.log10(v) - yLo) / (yHi - yLo)) * plotH;

  const svg = svgRoot(width, height, options.title);
  svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title || ""));

  for (const span of options.spans || []) {
    svg.appendChild(el("rect", { x: sx(span.start), y: padTop,
                                 width: Math.max(sx(span.end) - sx(span.start), 1), height: plotH,
                                 fill: span.colour || "#0b5c8a", "fill-opacity": 0.055 }));
    svg.appendChild(el("text", { x: (sx(span.start) + sx(span.end)) / 2, y: padTop - 3,
                                 "text-anchor": "middle", "font-size": 9, fill: "#5c6670" },
                       span.label || ""));
  }

  // Labelled gene band. The gene containing the marker is filled solid so the eye lands on
  // it; the rest are tinted. A box too narrow for an inside label gets the label above it.
  for (const g of band) {
    const gx = sx(g.start), gw = Math.max(sx(g.end) - sx(g.start), 2);
    const gy = 20 + g._row * (BOX_H + ROW_GAP);
    const hot = Number.isFinite(options.marker) &&
                options.marker >= g.start && options.marker <= g.end;
    svg.appendChild(el("rect", { x: gx, y: gy, width: gw, height: BOX_H, rx: 2,
                                 fill: hot ? "#2f4b7c" : "#ccd8e6",
                                 stroke: hot ? "#22395e" : "#9fb0c6", "stroke-width": 0.7 }));
    // Faint guide down the plot marks where this gene starts, so the band and the curve line up.
    svg.appendChild(el("line", { x1: gx, x2: gx, y1: padTop, y2: padTop + plotH,
                                 stroke: "#c9d4e2", "stroke-width": 0.6, "stroke-dasharray": "2 4" }));
    if (gw >= 20) {
      svg.appendChild(el("text", { x: gx + gw / 2, y: gy + BOX_H - 3, "text-anchor": "middle",
                                   "font-size": 8.5, fill: hot ? "#ffffff" : "#3c4a5c" },
                         g.label || ""));
    } else {
      svg.appendChild(el("text", { x: gx + gw / 2, y: gy - 2, "text-anchor": "middle",
                                   "font-size": 8, fill: "#3c4a5c" }, g.label || ""));
    }
  }

  // Variant ticks: where a submitted sequence differs from NL4-3. Drawn on the baseline
  // rather than over the curve, so a dense cluster reads as a band without hiding the rates.
  if (options.variants && options.variants.length) {
    const vy = padTop + plotH;
    for (const v of options.variants) {
      if (v < x0 || v > x1) continue;
      svg.appendChild(el("line", { x1: sx(v), x2: sx(v), y1: vy - 6, y2: vy,
                                   stroke: "#b5485d", "stroke-width": 0.8, "stroke-opacity": 0.75 }));
    }
    svg.appendChild(el("text", { x: padLeft + plotW + 6, y: vy, "font-size": 8.5,
                                 fill: "#b5485d" },
                       `${options.variants.length} difference${options.variants.length === 1 ? "" : "s"}`));
  }

  for (const e of niceLogTicks(yLo, yHi)) {
    const y = sy(Math.pow(10, e));
    svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: y, y2: y,
                                 stroke: "#e3e8ec", "stroke-width": 0.8 }));
    svg.appendChild(el("text", { x: padLeft - 6, y: y + 3.5, "text-anchor": "end",
                                 "font-size": 10, fill: "#5c6670" }, supLabel(e)));
  }

  // Inter-sample envelope, drawn under the series so the median line stays readable. Points
  // are [position, lo, hi]; a gap (lo or hi non-positive) breaks the polygon rather than
  // bridging it, so an uncovered stretch does not read as a narrow spread.
  if (ribbon && ribbon.length) {
    let runFwd = [], runBack = [];
    const flush = () => {
      if (runFwd.length > 1) {
        svg.appendChild(el("path", {
          d: "M" + runFwd.join(" L") + " L" + runBack.reverse().join(" L") + " Z",
          fill: options.ribbonColour || "#0b5c8a", "fill-opacity": options.ribbonOpacity || 0.18,
          stroke: "none" }));
      }
      runFwd = []; runBack = [];
    };
    for (const r of ribbon) {
      if (!(r[1] > 0 && r[2] > 0) || r[0] < x0 || r[0] > x1) { flush(); continue; }
      runFwd.push(`${sx(r[0]).toFixed(1)} ${sy(r[2]).toFixed(1)}`);
      runBack.push(`${sx(r[0]).toFixed(1)} ${sy(r[1]).toFixed(1)}`);
    }
    flush();
    if (options.ribbonLabel) {
      svg.appendChild(el("text", { x: padLeft + plotW + 6, y: padTop + 11, "font-size": 9,
                                   fill: options.ribbonColour || "#0b5c8a" }, options.ribbonLabel));
    }
  }

  for (const s of series) {
    const pts = s.points.filter((p) => p[1] > 0 && p[0] >= x0 && p[0] <= x1);
    if (!pts.length) continue;
    const d = pts.map((p, i) => (i ? "L" : "M") + sx(p[0]).toFixed(1) + " " + sy(p[1]).toFixed(1)).join(" ");
    const attrs = { d, fill: "none", stroke: s.colour || "#0b5c8a",
                    "stroke-width": s.width || 1, "stroke-opacity": s.opacity || 1 };
    // A dashed series reads as a comparator rather than a measurement, which is what the
    // NL4-3 overlay on the sequence page is.
    if (s.dash) attrs["stroke-dasharray"] = s.dash;
    svg.appendChild(el("path", attrs));
    const last = pts[pts.length - 1];
    svg.appendChild(el("text", { x: padLeft + plotW + 6, y: sy(last[1]) + 3.5, "font-size": 10,
                                 fill: s.colour || "#0b5c8a" }, s.label));
  }

  const axisY = padTop + plotH;
  svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: axisY, y2: axisY, stroke: "#999" }));
  const step = Math.max(1, Math.round((x1 - x0) / 8 / 500) * 500);
  for (let p = x0; p <= x1; p += step) {
    svg.appendChild(el("line", { x1: sx(p), x2: sx(p), y1: axisY, y2: axisY + 4, stroke: "#999" }));
    svg.appendChild(el("text", { x: sx(p), y: axisY + 15, "text-anchor": "middle",
                                 "font-size": 10, fill: "#5c6670" }, p.toLocaleString()));
  }
  svg.appendChild(el("text", { x: padLeft + plotW / 2, y: axisY + 30, "text-anchor": "middle",
                               "font-size": 10, fill: "#5c6670" },
                     options.xLabel || "NL4-3 reference position"));
  svg.appendChild(el("text", { x: 11, y: padTop + plotH / 2, "font-size": 10, fill: "#5c6670",
                               "text-anchor": "middle",
                               transform: `rotate(-90 11 ${padTop + plotH / 2})` },
                     options.yLabel || "rate per site per replication"));

  // Optional marker at a selected position, drawn last so it sits above the series.
  if (Number.isFinite(options.marker) && options.marker >= x0 && options.marker <= x1) {
    const mx = sx(options.marker);
    svg.appendChild(el("line", { x1: mx, x2: mx, y1: padTop, y2: padTop + plotH,
                                 stroke: "#b5485d", "stroke-width": 1.2,
                                 "stroke-dasharray": "4 3" }));
    svg.appendChild(el("text", { x: mx + 4, y: padTop + 10, "font-size": 9, fill: "#b5485d" },
                       options.markerLabel || String(options.marker)));
  }

  // Optional click-to-select over the plot area. The overlay is transparent and sits on
  // top, so a click anywhere in the plot maps to a reference position.
  if (typeof options.onSelect === "function") {
    const hit = el("rect", { x: padLeft, y: padTop, width: plotW, height: plotH,
                             fill: "transparent", cursor: "crosshair" });
    hit.addEventListener("click", (ev) => {
      const box = svg.getBoundingClientRect();
      // getBoundingClientRect is in CSS pixels; the SVG may be scaled by its viewBox.
      const scale = box.width / width;
      const px = (ev.clientX - box.left) / scale;
      const pos = Math.round(x0 + ((px - padLeft) / plotW) * (x1 - x0));
      options.onSelect(Math.min(Math.max(pos, x0), x1));
    });
    svg.appendChild(hit);
  }
  container.appendChild(svg);
}

// ---------------------------------------------------------------------------------
/**
 * Signed fold-change dot plot, centred at 0 because 0 is the semantic zero for a log2
 * ratio -- a change that is neither more nor less accessible.
 * rows: [{label, value, note?}]
 */
export function dotPlot(container, rows, options = {}) {
  container.innerHTML = "";
  if (!rows.length) {
    container.appendChild(Object.assign(document.createElement("p"),
      { className: "caveat", textContent: options.emptyMessage || "Nothing to plot." }));
    return;
  }
  const rowH = 16, padTop = 26, padBottom = 36, padLeft = 108, padRight = 96;
  const width = options.width || 720;
  const height = padTop + rows.length * rowH + padBottom;
  const plotW = width - padLeft - padRight;
  const span = Math.max(0.5, ...rows.map((r) => Math.abs(r.value)));
  const sx = (v) => padLeft + ((v + span) / (2 * span)) * plotW;
  const svg = svgRoot(width, height, options.title);
  svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title || ""));
  svg.appendChild(el("line", { x1: sx(0), x2: sx(0), y1: padTop - 4,
                               y2: padTop + rows.length * rowH, stroke: "#999",
                               "stroke-dasharray": "3 2" }));
  rows.forEach((r, i) => {
    const y = padTop + i * rowH + rowH / 2;
    svg.appendChild(el("text", { x: padLeft - 8, y: y + 3.5, "text-anchor": "end" }, r.label));
    svg.appendChild(el("line", { x1: sx(0), x2: sx(r.value), y1: y, y2: y,
                                 stroke: "#c8ced4", "stroke-width": 0.8 }));
    svg.appendChild(el("circle", { cx: sx(r.value), cy: y, r: 3.4,
                                   fill: r.value >= 0 ? "#0b5c8a" : "#c2561a" }));
    if (r.note) {
      svg.appendChild(el("text", { x: padLeft + plotW + 6, y: y + 3.5, "font-size": 10,
                                   fill: "#5c6670" }, r.note));
    }
  });
  const axisY = padTop + rows.length * rowH + 4;
  svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: axisY, y2: axisY, stroke: "#999" }));
  for (const v of [-span, -span / 2, 0, span / 2, span]) {
    svg.appendChild(el("line", { x1: sx(v), x2: sx(v), y1: axisY, y2: axisY + 4, stroke: "#999" }));
    svg.appendChild(el("text", { x: sx(v), y: axisY + 15, "text-anchor": "middle",
                                 "font-size": 10, fill: "#5c6670" }, v.toFixed(1)));
  }
  svg.appendChild(el("text", { x: padLeft + plotW / 2, y: axisY + 30, "text-anchor": "middle",
                               "font-size": 10, fill: "#5c6670" },
                     options.xLabel || "log2 fold change vs NL4-3  (left = less accessible)"));
  container.appendChild(svg);
}

// ---------------------------------------------------------------------------------
/**
 * Pie chart of mutational supply across a set of amino-acid changes. Each slice is one
 * target amino acid; the slice angle is its share of the total supply. A side legend
 * carries the absolute rate and percentage, so the slices stay uncluttered.
 *
 * slices: [{label, value, highlight?}]  (value = mu_aa; non-positive values dropped)
 */
export function pieChart(container, slices, options = {}) {
  container.innerHTML = "";
  const data = slices.filter((s) => s.value > 0).sort((a, b) => b.value - a.value);
  if (!data.length) {
    container.appendChild(Object.assign(document.createElement("p"),
      { className: "caveat", textContent: options.emptyMessage || "No changes to plot." }));
    return;
  }
  const total = data.reduce((t, s) => t + s.value, 0);
  const width = options.width || 680, height = options.height || 300;
  // The circle is sized from the box rather than fixed, so asking for a bigger figure makes
  // the pie bigger instead of just adding white space to the right of it.
  const legendW = options.legendWidth || 215;
  const r = options.radius ||
            Math.max(80, Math.min((height - 58) / 2, (width - legendW - 56) / 2));
  const cx = 26 + r;
  const cy = Math.max(r + 26, height / 2);
  const fsSlice = r >= 130 ? 13 : 11;
  const fsLegend = r >= 130 ? 12.5 : 11;
  const legendStep = r >= 130 ? 19 : 17;

  // A slice may carry its own colour. That is how the position pie is kept in step with the
  // supply-flow chart beside it: both are coloured by destination residue, so a reader can
  // carry one colour between the two figures instead of learning a new key for each. Without
  // a per-slice colour the palette is positional, which means the same residue changes colour
  // as soon as the ordering changes -- actively misleading when two figures sit side by side.
  const palette = ["#0b5c8a", "#c2561a", "#1f7a3f", "#6b3fa0", "#8a6d1f", "#2f7fa8",
                   "#a4262c", "#5c6670", "#c08a1e", "#3a7d5d", "#7a4ea0", "#1d6a8a"];
  const colour = (i, s) => s.colour || (s.highlight ? "#c2561a" : palette[i % palette.length]);
  // Slice selection: `options.onSelect(label)` makes the wedges and legend rows clickable,
  // and `options.selected` names the one drawn as chosen.
  const pick = options.onSelect;
  const isSel = (s) => options.selected != null && s.label === options.selected;
  const wire = (node, s) => {
    if (!pick) return node;
    node.style.cursor = "pointer";
    node.addEventListener("click", () => pick(s.label));
    node.addEventListener("keydown", (e) => {
      if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(s.label); }
    });
    node.setAttribute("tabindex", "0");
    node.setAttribute("role", "button");
    node.setAttribute("aria-label", `Select ${s.label}`);
    return node;
  };

  const svg = svgRoot(width, height, options.title);
  if (options.title) svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title));

  let angle = -Math.PI / 2;   // start at 12 o'clock
  data.forEach((s, i) => {
    const frac = s.value / total;
    const a0 = angle, a1 = angle + frac * 2 * Math.PI;
    angle = a1;
    const large = a1 - a0 > Math.PI ? 1 : 0;
    const x0 = cx + r * Math.cos(a0), y0 = cy + r * Math.sin(a0);
    const x1 = cx + r * Math.cos(a1), y1 = cy + r * Math.sin(a1);
    // a single 100% slice can't be drawn as an arc (start==end); draw a full circle
    if (data.length === 1) {
      svg.appendChild(wire(el("circle", { cx, cy, r, fill: colour(i, s),
        "fill-opacity": isSel(s) ? 1 : 0.85,
        stroke: isSel(s) ? "#1b1f23" : "none", "stroke-width": isSel(s) ? 2 : 0 }), s));
    } else {
      const d = `M ${cx} ${cy} L ${x0.toFixed(2)} ${y0.toFixed(2)} ` +
                `A ${r} ${r} 0 ${large} 1 ${x1.toFixed(2)} ${y1.toFixed(2)} Z`;
      svg.appendChild(wire(el("path", { d, fill: colour(i, s),
        "fill-opacity": isSel(s) ? 1 : 0.85,
        stroke: isSel(s) ? "#1b1f23" : "#fff", "stroke-width": isSel(s) ? 2 : 1 }), s));
    }
    // label the slice with the amino acid if it is big enough to read
    if (frac > 0.05) {
      const mid = (a0 + a1) / 2;
      const lx = cx + r * 0.62 * Math.cos(mid), ly = cy + r * 0.62 * Math.sin(mid);
      svg.appendChild(el("text", { x: lx, y: ly + 3, "text-anchor": "middle",
                                   "font-size": fsSlice,
                                   fill: "#fff", "font-weight": 600 }, s.label));
    }
  });

  // legend
  const lx = cx + r + 32;
  let ly = cy - data.length * legendStep / 2 + legendStep / 2;
  data.forEach((s, i) => {
    const g = el("g", {});
    if (isSel(s)) {
      g.appendChild(el("rect", { x: lx - 4, y: ly - legendStep + 5, width: legendW - 10,
                                 height: legendStep, rx: 2, fill: "#eef4f8" }));
    }
    g.appendChild(el("rect", { x: lx, y: ly - 8, width: 10, height: 10,
                               fill: colour(i, s), "fill-opacity": 0.85 }));
    const pct = ((s.value / total) * 100).toFixed(1);
    g.appendChild(el("text", { x: lx + 15, y: ly, "font-size": fsLegend, fill: "#1b1f23",
                               "font-weight": isSel(s) ? 600 : 400 },
      `${s.label}  ${pct}%  (${s.value.toExponential(2)})`));
    svg.appendChild(wire(g, s));
    ly += legendStep;
  });
  container.appendChild(svg);
}

