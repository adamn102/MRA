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
  const svg = svgRoot(width, height, options.title);

  svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title || ""));

  rows.forEach((r, i) => {
    const y = padTop + i * rowH;
    svg.appendChild(el("text", { x: padLeft - 6, y: y + 11, "text-anchor": "end",
                                 "font-weight": r.highlight ? 600 : 400 }, r.label));
    const frac = max > 0 ? (r.value || 0) / max : 0;
    const w = Math.max(frac * plotW, r.value > 0 ? 1.5 : 0);
    svg.appendChild(el("rect", { x: padLeft, y: y + 3, width: w, height: rowH - 7,
                                 fill: r.highlight ? "#c2561a" : "#0b5c8a",
                                 "fill-opacity": r.highlight ? 0.85 : 0.6 }));
    // Interval as a whisker on the same linear scale.
    if (r.lo != null && r.hi != null && max > 0) {
      const x1 = padLeft + (r.lo / max) * plotW, x2 = padLeft + (r.hi / max) * plotW;
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
    svg.appendChild(el("text", { x, y: axisY + 15, "text-anchor": "middle", "font-size": 10,
                                 fill: "#5c6670" },
                       frac === 0 ? "0" : (frac * max).toExponential(1)));
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
  const width = options.width || 1020, height = options.height || 260;
  const padLeft = 62, padRight = 118, padTop = 24, padBottom = 40;
  const plotW = width - padLeft - padRight, plotH = height - padTop - padBottom;
  const x0 = options.xMin || 1, x1 = options.xMax || 8213;

  const all = series.flatMap((s) => s.points.map((p) => p[1])).filter((v) => v > 0);
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

  for (const e of niceLogTicks(yLo, yHi)) {
    const y = sy(Math.pow(10, e));
    svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: y, y2: y,
                                 stroke: "#e3e8ec", "stroke-width": 0.8 }));
    svg.appendChild(el("text", { x: padLeft - 6, y: y + 3.5, "text-anchor": "end",
                                 "font-size": 10, fill: "#5c6670" }, supLabel(e)));
  }

  for (const s of series) {
    const pts = s.points.filter((p) => p[1] > 0 && p[0] >= x0 && p[0] <= x1);
    if (!pts.length) continue;
    const d = pts.map((p, i) => (i ? "L" : "M") + sx(p[0]).toFixed(1) + " " + sy(p[1]).toFixed(1)).join(" ");
    svg.appendChild(el("path", { d, fill: "none", stroke: s.colour || "#0b5c8a",
                                 "stroke-width": s.width || 1, "stroke-opacity": s.opacity || 1 }));
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
  const cx = 150, cy = height / 2, r = 110;

  // A categorical palette; the highlighted slice (if any) takes the lead colour.
  const palette = ["#0b5c8a", "#c2561a", "#1f7a3f", "#6b3fa0", "#8a6d1f", "#2f7fa8",
                   "#a4262c", "#5c6670", "#c08a1e", "#3a7d5d", "#7a4ea0", "#1d6a8a"];
  const colour = (i, s) => s.highlight ? "#c2561a" : palette[i % palette.length];

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
      svg.appendChild(el("circle", { cx, cy, r, fill: colour(i, s), "fill-opacity": 0.85 }));
    } else {
      const d = `M ${cx} ${cy} L ${x0.toFixed(2)} ${y0.toFixed(2)} ` +
                `A ${r} ${r} 0 ${large} 1 ${x1.toFixed(2)} ${y1.toFixed(2)} Z`;
      svg.appendChild(el("path", { d, fill: colour(i, s), "fill-opacity": 0.85,
                                   stroke: "#fff", "stroke-width": 1 }));
    }
    // label the slice with the amino acid if it is big enough to read
    if (frac > 0.05) {
      const mid = (a0 + a1) / 2;
      const lx = cx + r * 0.62 * Math.cos(mid), ly = cy + r * 0.62 * Math.sin(mid);
      svg.appendChild(el("text", { x: lx, y: ly + 3, "text-anchor": "middle", "font-size": 11,
                                   fill: "#fff", "font-weight": 600 }, s.label));
    }
  });

  // legend
  const lx = cx + r + 40;
  let ly = cy - data.length * 8;
  data.forEach((s, i) => {
    svg.appendChild(el("rect", { x: lx, y: ly - 8, width: 10, height: 10,
                                 fill: colour(i, s), "fill-opacity": 0.85 }));
    const pct = ((s.value / total) * 100).toFixed(1);
    svg.appendChild(el("text", { x: lx + 15, y: ly, "font-size": 11, fill: "#1b1f23" },
      `${s.label}  ${pct}%  (${s.value.toExponential(2)})`));
    ly += 17;
  });
  container.appendChild(svg);
}
