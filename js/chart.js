// SPDX-License-Identifier: MIT
//
// chart.js -- the charts this portal draws, in plain SVG. No dependencies.
//
//   genomeTrack  rate along the genome, optionally with a second series overlaid
//   pieChart     one-step amino-acid supply at one codon
//
// SVG rather than canvas so the output is inspectable, scales on a high-DPI screen,
// and can be saved by right-click. Rates span orders of magnitude, so the genome track
// defaults to a log10 axis fitted to the values in view. A caller can ask for a linear
// axis, or for the wider fixed axis (powers of ten, or zero to the largest rate).

const NS = "http://www.w3.org/2000/svg";

/**
 * Two range inputs that share one track and never cross. `onChange(lo, hi)` fires while a
 * thumb moves. `set` updates both thumbs without firing, so a gene click or a double-click
 * can move the window without looping back into itself.
 */
export function bindZoomRange(loEl, hiEl, { fill, readout, onChange } = {}) {
  const paint = () => {
    const min = Number(loEl.min);
    const span = (Number(loEl.max) - min) || 1;
    const a = (Number(loEl.value) - min) / span;
    const b = (Number(hiEl.value) - min) / span;
    if (fill) {
      fill.style.left = (a * 100) + "%";
      fill.style.width = Math.max((b - a) * 100, 0) + "%";
    }
    if (readout) {
      readout.textContent = Number(loEl.value).toLocaleString() + "\u2013" +
        Number(hiEl.value).toLocaleString();
    }
  };
  const clamp = (which) => {
    const min = Number(loEl.min), max = Number(loEl.max);
    let lo = Number(loEl.value), hi = Number(hiEl.value);
    if (lo >= hi) {
      if (which === "lo") lo = hi - 1;
      else hi = lo + 1;
    }
    lo = Math.max(min, Math.min(lo, max - 1));
    hi = Math.max(min + 1, Math.min(hi, max));
    if (lo >= hi) lo = hi - 1;
    loEl.value = String(lo);
    hiEl.value = String(hi);
    paint();
    if (onChange) onChange(lo, hi);
  };
  loEl.addEventListener("input", () => clamp("lo"));
  hiEl.addEventListener("input", () => clamp("hi"));
  loEl.addEventListener("pointerdown", () => { loEl.style.zIndex = "5"; hiEl.style.zIndex = "4"; });
  hiEl.addEventListener("pointerdown", () => { hiEl.style.zIndex = "5"; loEl.style.zIndex = "4"; });
  hiEl.style.zIndex = "4";
  paint();
  return {
    set(lo, hi) {
      loEl.value = String(lo);
      hiEl.value = String(hi);
      paint();
    },
  };
}

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

/** A 1-2-5 step for a value span. Unlike `niceStep`, this is allowed to be far below 1. */
function niceValueStep(span) {
  const raw = span / 5;
  if (!(raw > 0)) return 1;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const err = raw / pow;
  const nice = err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1;
  return nice * pow;
}

/** Tick values on a linear axis, including the ends when they land on the step. */
function linearTickValues(yMin, yMax) {
  const span = yMax - yMin;
  if (!(span > 0)) return [];
  const step = niceValueStep(span);
  const ticks = [];
  const first = Math.ceil((yMin - step * 1e-8) / step);
  const last = Math.floor((yMax + step * 1e-8) / step);
  for (let i = first; i <= last && ticks.length < 8; i++) {
    const v = i * step;
    if (v < yMin - step * 1e-6 || v > yMax + step * 1e-6) continue;
    ticks.push(v);
  }
  return ticks.length ? ticks : [(yMin + yMax) / 2];
}

/**
 * Positive tick values on a fitted log axis. Whole powers of ten when the span is wide;
 * 1-2-5 (then finer) mantissas when the data occupies less than a decade, so a tight
 * axis still has several gridlines.
 */
function fittedLogTicks(yLo, yHi) {
  const span = yHi - yLo;
  const mantissas = span >= 2 ? [1]
    : span >= 1.15 ? [1, 2, 5]
    : span >= 0.5 ? [1, 2, 3, 5, 7]
    : [1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9];
  const ticks = [];
  for (let e = Math.floor(yLo) - 1; e <= Math.ceil(yHi); e++) {
    for (const m of mantissas) {
      const lv = Math.log10(m) + e;
      if (lv <= yLo + span * 0.03 || lv >= yHi - span * 0.03) continue;
      ticks.push(m * Math.pow(10, e));
    }
  }
  return ticks;
}

/** A 1-2-5 step that puts about six labels across a nucleotide span. */
function niceStep(span) {
  const raw = Math.max(span, 1) / 6;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const err = raw / pow;
  const nice = err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1;
  return Math.max(1, nice * pow);
}

const supLabel = (e) => "10" + String(e).replace(/-/g, "\u2212").split("").map(
  (c) => ({ "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3", "4": "\u2074",
            "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078", "9": "\u2079",
            "\u2212": "\u207b" }[c] || c)).join("");

/** Compact tick text. Powers of ten stay as 10ⁿ; other values are a one-decimal mantissa. */
function formatAxisValue(v) {
  if (!(v > 0)) return "0";
  const e = Math.floor(Math.log10(v) + 1e-12);
  const m = v / Math.pow(10, e);
  const mRound = Math.round(m * 10) / 10;
  if (Math.abs(mRound - 1) < 1e-6) return supLabel(e);
  const mStr = Math.abs(mRound - Math.round(mRound)) < 1e-6
    ? String(Math.round(mRound))
    : mRound.toFixed(1);
  return mStr + "\u00d7" + supLabel(e);
}

// Survives a redraw, so the second click of a double-click is still recognisable after the
// first click has replaced the plot. A lone click waits briefly, and is cancelled if the
// second click arrives, so zooming back out does not also move the marker.
let plotPointer = { t: 0, x: 0, y: 0, timer: 0 };
function plotClickOrReset(ev, pos, onSelect, onZoomReset) {
  const now = performance.now();
  const second = plotPointer.t && (now - plotPointer.t) < 400 &&
    Math.hypot(ev.clientX - plotPointer.x, ev.clientY - plotPointer.y) < 12;
  clearTimeout(plotPointer.timer);
  plotPointer.timer = 0;
  if (second && typeof onZoomReset === "function") {
    plotPointer.t = 0;
    onZoomReset();
    return;
  }
  plotPointer.t = now;
  plotPointer.x = ev.clientX;
  plotPointer.y = ev.clientY;
  if (typeof onSelect !== "function") return;
  if (typeof onZoomReset !== "function") { onSelect(pos); return; }
  plotPointer.timer = setTimeout(() => {
    plotPointer.timer = 0;
    onSelect(pos);
  }, 280);
}

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

// ---------------------------------------------------------------------------------
/**
 * Rate along the genome. series: [{label, colour, points: [[position, rate], ...]}].
 * The value axis is log10 unless `yScale` is "linear". `autoScale` (default true) fits
 * the axis to the rates inside the current x window; false keeps a fixed axis — whole
 * powers of ten on a log scale, or zero to the largest rate on a linear scale.
 * Optional shaded spans for gene boundaries.
 */
export function genomeTrack(container, series, options = {}) {
  container.innerHTML = "";
  const width = options.width || LAYOUT.width;
  // Extra room between the gene band and the rate plot, so a position slider can sit there
  // without covering either. The plot area keeps the height it was given.
  const sliderGap = Math.max(0, options.sliderGap || 0);
  const height = (options.height || 260) + sliderGap;
  let padLeft = LAYOUT.padLeft;
  // One series is named inside the plot. Several series need a legend in the right margin,
  // and the margin has to be wide enough for the longest name or the text is clipped.
  const sideLegend = series.length > 1;
  const legendText = sideLegend
    ? Math.max(...series.map((s) => (s.label || "").length), 1) * 7.1
    : 0;
  const padRight = sideLegend
    ? Math.max(LAYOUT.padRight, Math.ceil(36 + legendText) + 10)
    : LAYOUT.padRight;
  const padBottom = 40;
  const x0 = options.xMin || 1;
  let x1 = options.xMax || 8213;
  if (!(x1 > x0)) x1 = x0 + 1;

  // An explicit gene band (labelled boxes above the plot) is opt-in via `geneBand`, so the
  // older faint `spans` shading used elsewhere keeps its appearance. Overlapping reading
  // frames -- gag/pol, and tat/rev/vpu/env -- are packed onto as many rows as they need,
  // because drawing them on one row would either hide genes or imply they do not overlap.
  const band = (options.geneBand || []).slice().sort((a, b) => a.start - b.start);
  // Rev is spliced. Its first exon is short enough to fit in the gap between vpr and env,
  // which would leave the two exons on different rows. Both belong on the bottom track.
  const pinBottom = new Set(["Rev"]);
  const bandRows = [];
  for (const g of band) {
    if (pinBottom.has(g.label)) continue;
    let r = 0;
    while (bandRows[r] && bandRows[r].some((o) => !(g.start > o.end || g.end < o.start))) r++;
    (bandRows[r] = bandRows[r] || []).push(g);
    g._row = r;
  }
  const pinned = band.filter((g) => pinBottom.has(g.label));
  if (pinned.length) {
    const r = bandRows.length;
    (bandRows[r] = bandRows[r] || []).push(...pinned);
    for (const g of pinned) g._row = r;
  }
  const BOX_H = 12, ROW_GAP = 2;
  const bandH = bandRows.length ? bandRows.length * (BOX_H + ROW_GAP) + 6 : 0;
  const padTop = 24 + bandH + sliderGap;
  const plotH = height - padTop - padBottom;

  // The ribbon's extent counts toward the y-domain: an envelope drawn outside the computed
  // range would be silently clipped at the plot edge and read as a narrower spread.
  const ribbon = options.ribbon || null;
  const yScale = options.yScale === "linear" ? "linear" : "log";
  const autoScale = options.autoScale !== false;
  // Auto-scale reads only the window on screen, which is the data being plotted. The fixed
  // axis reads every point, so zooming does not move the scale.
  const collectRates = (restrict) => {
    const vals = [];
    for (const s of series) {
      for (const p of s.points) {
        if (!(p[1] > 0)) continue;
        if (restrict && (p[0] < x0 || p[0] > x1)) continue;
        vals.push(p[1]);
      }
    }
    if (ribbon) {
      for (const r of ribbon) {
        if (restrict && (r[0] < x0 || r[0] > x1)) continue;
        if (r[1] > 0) vals.push(r[1]);
        if (r[2] > 0) vals.push(r[2]);
      }
    }
    return vals;
  };
  let basis = autoScale ? collectRates(true) : collectRates(false);
  if (!basis.length) basis = collectRates(false);
  if (!basis.length) {
    container.appendChild(Object.assign(document.createElement("p"),
      { className: "caveat", textContent: "No rates to plot for this selection." }));
    return;
  }
  const dataMin = Math.min(...basis), dataMax = Math.max(...basis);
  let sy;
  let yTicks;
  if (yScale === "linear") {
    let yMin, yMax;
    if (autoScale) {
      const span = dataMax - dataMin;
      const pad = span > 0 ? span * 0.12 : (dataMax || 1) * 0.08;
      yMin = Math.max(0, dataMin - pad);
      yMax = dataMax + pad;
      if (!(yMax > yMin)) yMax = yMin + (dataMax || 1) * 0.1;
    } else {
      const every = collectRates(false);
      const mx = Math.max(...every, dataMax);
      const step = niceValueStep(mx || 1);
      yMin = 0;
      yMax = Math.max(step, Math.ceil((mx * 1.02) / step) * step);
    }
    sy = (v) => padTop + plotH - ((v - yMin) / (yMax - yMin)) * plotH;
    yTicks = linearTickValues(yMin, yMax).map((v) => ({ v, label: formatAxisValue(v) }));
  } else {
    const lo = Math.log10(dataMin), hi = Math.log10(dataMax);
    let yLo, yHi;
    if (autoScale) {
      const span = hi - lo;
      if (span < 0.12) {
        const mid = (lo + hi) / 2;
        yLo = mid - 0.08;
        yHi = mid + 0.08;
      } else {
        const pad = span * 0.12;
        yLo = lo - pad;
        yHi = hi + pad;
      }
    } else {
      const every = collectRates(false);
      const a = Math.log10(Math.min(...every));
      const b = Math.log10(Math.max(...every));
      yLo = Math.floor(a);
      yHi = Math.ceil(b);
      if (yHi <= yLo) yHi = yLo + 1;
    }
    sy = (v) => padTop + plotH - ((Math.log10(v) - yLo) / (yHi - yLo)) * plotH;
    const powers = autoScale ? fittedLogTicks(yLo, yHi) : niceLogTicks(yLo, yHi).map((e) => Math.pow(10, e));
    yTicks = powers.map((v) => ({ v, label: formatAxisValue(v) }));
  }
  // A mantissa label such as 2.5×10⁻⁴ is wider than 10⁻⁴. Grow the left margin so it
  // stays clear of the rotated axis title, and so the position slider follows the plot.
  const tickChars = yTicks.reduce((m, t) => Math.max(m, t.label.length), 2);
  const needLeft = 24 + Math.ceil(tickChars * 6.5);
  if (needLeft > padLeft) padLeft = needLeft;
  const plotW = width - padLeft - padRight;
  const sx = (p) => padLeft + ((p - x0) / (x1 - x0)) * plotW;

  const svg = svgRoot(width, height, options.title);
  svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title || ""));

  for (const span of options.spans || []) {
    if (span.end < x0 || span.start > x1) continue;
    const ss = Math.max(span.start, x0), se = Math.min(span.end, x1);
    svg.appendChild(el("rect", { x: sx(ss), y: padTop,
                                 width: Math.max(sx(se) - sx(ss), 1), height: plotH,
                                 fill: span.colour || "#0b5c8a", "fill-opacity": 0.055 }));
    svg.appendChild(el("text", { x: (sx(ss) + sx(se)) / 2, y: padTop - 3,
                                 "text-anchor": "middle", "font-size": 9, fill: "#5c6670" },
                       span.label || ""));
  }

  // Labelled gene band. The gene containing the marker is filled solid so the eye lands on
  // it; the rest are tinted. A box too narrow for an inside label gets the label above it.
  for (const g of band) {
    if (g.end < x0 || g.start > x1) continue;
    const gs = Math.max(g.start, x0), ge = Math.min(g.end, x1);
    const gx = sx(gs), gw = Math.max(sx(ge) - sx(gs), 2);
    const gy = 20 + g._row * (BOX_H + ROW_GAP);
    const hot = Number.isFinite(options.marker) &&
                options.marker >= g.start && options.marker <= g.end;
    const box = el("rect", { x: gx, y: gy, width: gw, height: BOX_H, rx: 2,
                             fill: hot ? "#2f4b7c" : "#ccd8e6",
                             stroke: hot ? "#22395e" : "#9fb0c6", "stroke-width": 0.7 });
    // Faint guide down the plot marks where this gene starts, so the band and the curve line up.
    // A gene clipped by the zoom window starts off-screen; its true boundary is not at the edge.
    // The guide stays outside the click target, so a click on the plot still selects a position.
    if (g.start >= x0) {
      svg.appendChild(el("line", { x1: gx, x2: gx, y1: padTop, y2: padTop + plotH,
                                   stroke: "#c9d4e2", "stroke-width": 0.6, "stroke-dasharray": "2 4" }));
    }
    let labelNode = null;
    if (gw >= 20) {
      labelNode = el("text", { x: gx + gw / 2, y: gy + BOX_H - 3, "text-anchor": "middle",
                               "font-size": 8.5, fill: hot ? "#ffffff" : "#3c4a5c" },
                     g.label || "");
    } else if (g._row === 0) {
      labelNode = el("text", { x: gx + gw / 2, y: gy - 2, "text-anchor": "middle",
                               "font-size": 8, fill: "#3c4a5c" }, g.label || "");
    } else {
      // A label above this box would sit in the track above it. Keep it in its own row.
      labelNode = el("text", { x: gx + gw + 3, y: gy + BOX_H - 3, "text-anchor": "start",
                               "font-size": 8, fill: "#3c4a5c" }, g.label || "");
    }
    // The range is the gene's own nucleotides, not the clipped sliver currently on screen.
    if (typeof options.onZoom === "function") {
      const grp = el("g", { role: "button" });
      grp.style.cursor = "pointer";
      grp.setAttribute("aria-label", `Zoom to ${g.label || "gene"}`);
      grp.addEventListener("click", (ev) => {
        ev.stopPropagation();
        ev.preventDefault();
        options.onZoom(g.start, g.end);
      });
      grp.append(box, labelNode);
      svg.appendChild(grp);
    } else {
      svg.append(box, labelNode);
    }
  }

  // Variant ticks: where a submitted sequence differs from NL4-3. Drawn on the baseline
  // rather than over the curve, so a dense cluster reads as a band without hiding the rates.
  if (options.variants && options.variants.length) {
    const visible = options.variants.filter((v) => v >= x0 && v <= x1);
    const vy = padTop + plotH;
    for (const v of visible) {
      svg.appendChild(el("line", { x1: sx(v), x2: sx(v), y1: vy - 6, y2: vy,
                                   stroke: "#b5485d", "stroke-width": 0.8, "stroke-opacity": 0.75 }));
    }
    if (visible.length && !sideLegend) {
      svg.appendChild(el("text", { x: padLeft + plotW + 6, y: vy, "font-size": 8.5,
                                   fill: "#b5485d" },
                         `${visible.length} difference${visible.length === 1 ? "" : "s"}`));
    }
  }

  for (const tick of yTicks) {
    const y = sy(tick.v);
    if (y < padTop - 1 || y > padTop + plotH + 1) continue;
    svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: y, y2: y,
                                 stroke: "#e3e8ec", "stroke-width": 0.8 }));
    svg.appendChild(el("text", { x: padLeft - 6, y: y + 3.5, "text-anchor": "end",
                                 "font-size": 10, fill: "#5c6670" }, tick.label));
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

  const drawn = [];
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
    drawn.push({ s, last: pts[pts.length - 1] });
  }
  if (drawn.length === 1) {
    // Bottom-right corner of the plot, flush with the right edge and just above the axis.
    const { s } = drawn[0];
    svg.appendChild(el("text", {
      x: padLeft + plotW, y: padTop + plotH - 6, "text-anchor": "end", "font-size": 12,
      fill: s.colour || "#0b5c8a", stroke: "#ffffff", "stroke-width": 3.5,
      "paint-order": "stroke", "stroke-linejoin": "round",
    }, s.label || ""));
  } else if (drawn.length > 1) {
    const lx = padLeft + plotW + 12;
    const rowH = 16;
    let ly = padTop + 14;
    for (const { s } of drawn) {
      const colour = s.colour || "#0b5c8a";
      const swatch = el("line", {
        x1: lx, x2: lx + 16, y1: ly - 4, y2: ly - 4,
        stroke: colour, "stroke-width": Math.max(s.width || 2, 2), "stroke-linecap": "round",
      });
      if (s.dash) swatch.setAttribute("stroke-dasharray", s.dash);
      svg.appendChild(swatch);
      svg.appendChild(el("text", { x: lx + 22, y: ly, "font-size": 12, fill: colour },
                         s.label || ""));
      ly += rowH;
    }
  }

  const axisY = padTop + plotH;
  svg.appendChild(el("line", { x1: padLeft, x2: padLeft + plotW, y1: axisY, y2: axisY, stroke: "#999" }));
  const step = niceStep(x1 - x0);
  const ticks = [];
  for (let p = Math.ceil(x0 / step) * step; p <= x1 + 1e-6; p += step) ticks.push(p);
  const near = (p) => ticks.some((t) => Math.abs(t - p) < step * 0.45);
  if (!near(x0)) ticks.unshift(x0);
  if (!near(x1)) ticks.push(x1);
  for (const p of ticks) {
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

  // Hover readout of the smoothed rate on every line at one nucleotide. The guide and the
  // dots are omitted from a PNG; the box itself sits outside the SVG.
  const hoverSeries = [];
  for (const s of series) {
    const rates = new Map();
    for (const p of s.points) if (p[1] > 0) rates.set(p[0], p[1]);
    if (rates.size) hoverSeries.push({
      label: s.name || s.label || "",
      colour: s.colour || "#0b5c8a",
      rates,
    });
  }
  const guide = el("line", {
    y1: padTop, y2: padTop + plotH, stroke: "#5c6670", "stroke-width": 1,
    "stroke-dasharray": "2 2", visibility: "hidden", "pointer-events": "none",
  });
  guide.setAttribute("data-export-skip", "");
  const dots = hoverSeries.map((s) => {
    const c = el("circle", {
      r: 3.2, fill: s.colour, stroke: "#ffffff", "stroke-width": 1,
      visibility: "hidden", "pointer-events": "none",
    });
    c.setAttribute("data-export-skip", "");
    return c;
  });
  svg.appendChild(guide);
  for (const c of dots) svg.appendChild(c);
  if (getComputedStyle(container).position === "static") container.style.position = "relative";
  const tip = document.createElement("div");
  tip.className = "track-tip";
  tip.setAttribute("data-export-skip", "");
  tip.hidden = true;
  container.appendChild(tip);
  const hideHover = () => {
    guide.setAttribute("visibility", "hidden");
    for (const c of dots) c.setAttribute("visibility", "hidden");
    tip.hidden = true;
  };
  const showHover = (pos, ev) => {
    const x = sx(pos);
    guide.setAttribute("x1", x);
    guide.setAttribute("x2", x);
    guide.setAttribute("visibility", "visible");
    tip.replaceChildren();
    const head = document.createElement("div");
    head.className = "track-tip-pos";
    head.textContent = "NL4-3 " + pos.toLocaleString();
    const unit = document.createElement("div");
    unit.className = "track-tip-unit";
    unit.textContent = "supply per replication";
    tip.append(head, unit);
    const table = document.createElement("table");
    hoverSeries.forEach((s, i) => {
      const rate = s.rates.get(pos);
      const tr = document.createElement("tr");
      const name = document.createElement("td");
      const sw = document.createElement("i");
      sw.className = "track-tip-swatch";
      sw.style.background = s.colour;
      name.append(sw, document.createTextNode(s.label));
      name.title = s.label;
      const val = document.createElement("td");
      val.className = "num";
      val.textContent = rate > 0 ? rate.toExponential(2) : "\u2014";
      tr.append(name, val);
      table.appendChild(tr);
      const dot = dots[i];
      if (rate > 0) {
        dot.setAttribute("cx", x);
        dot.setAttribute("cy", sy(rate));
        dot.setAttribute("visibility", "visible");
      } else dot.setAttribute("visibility", "hidden");
    });
    tip.appendChild(table);
    tip.hidden = false;
    const host = container.getBoundingClientRect();
    let left = ev.clientX - host.left + 14;
    let top = ev.clientY - host.top + 14;
    tip.style.left = left + "px";
    tip.style.top = top + "px";
    if (left + tip.offsetWidth > host.width - 8) left = ev.clientX - host.left - tip.offsetWidth - 14;
    if (top + tip.offsetHeight > host.height - 8) top = ev.clientY - host.top - tip.offsetHeight - 10;
    tip.style.left = Math.max(4, left) + "px";
    tip.style.top = Math.max(4, top) + "px";
  };

  // Optional click-to-select over the plot area. The overlay is transparent and sits on
  // top, so a click anywhere in the plot maps to a reference position. When onZoom is
  // set, a drag across the curve zooms the x-axis to that span and a click still selects.
  // Moving the pointer reads the smoothed rate of every line at that nucleotide.
  if (hoverSeries.length || typeof options.onSelect === "function" || typeof options.onZoom === "function") {
    const hit = el("rect", { x: padLeft, y: padTop, width: plotW, height: plotH,
                             fill: "transparent", cursor: "crosshair" });
    const posAt = (ev) => {
      const box = svg.getBoundingClientRect();
      // getBoundingClientRect is in CSS pixels; the SVG may be scaled by its viewBox.
      const scale = box.width / width;
      const px = (ev.clientX - box.left) / scale;
      const pos = Math.round(x0 + ((px - padLeft) / plotW) * (x1 - x0));
      return Math.min(Math.max(pos, x0), x1);
    };
    let brushing = false;
    hit.addEventListener("pointermove", (ev) => {
      if (brushing) { hideHover(); return; }
      if (hoverSeries.length) showHover(posAt(ev), ev);
    });
    hit.addEventListener("pointerleave", hideHover);
    if (typeof options.onZoom === "function") {
      const brush = el("rect", {
        fill: "#2f4b7c", "fill-opacity": 0.16, stroke: "#2f4b7c", "stroke-width": 1,
        visibility: "hidden", "pointer-events": "none",
      });
      let startX = null, startPos = null;
      hit.addEventListener("pointerdown", (ev) => {
        if (ev.button !== 0) return;
        clearTimeout(plotPointer.timer);
        plotPointer.timer = 0;
        startX = ev.clientX;
        startPos = posAt(ev);
        brushing = false;
        try { hit.setPointerCapture(ev.pointerId); } catch { /* pointer is not active */ }
        ev.preventDefault();
      });
      hit.addEventListener("pointermove", (ev) => {
        if (startX === null || Math.abs(ev.clientX - startX) < 5) return;
        brushing = true;
        hideHover();
        const pos = posAt(ev);
        const xA = sx(Math.min(startPos, pos));
        const xB = sx(Math.max(startPos, pos));
        brush.setAttribute("x", xA);
        brush.setAttribute("y", padTop);
        brush.setAttribute("width", Math.max(xB - xA, 1));
        brush.setAttribute("height", plotH);
        brush.setAttribute("visibility", "visible");
      });
      const finish = (ev) => {
        if (startPos === null) return;
        const origin = startPos;
        const pos = posAt(ev);
        const dragged = brushing;
        startX = null;
        startPos = null;
        brushing = false;
        brush.setAttribute("visibility", "hidden");
        if (dragged && Math.abs(pos - origin) >= 2) {
          clearTimeout(plotPointer.timer);
          plotPointer.timer = 0;
          plotPointer.t = 0;
          options.onZoom(Math.min(origin, pos), Math.max(origin, pos));
        } else {
          plotClickOrReset(ev, pos, options.onSelect, options.onZoomReset);
        }
      };
      hit.addEventListener("pointerup", finish);
      hit.addEventListener("pointercancel", () => {
        startX = null;
        startPos = null;
        brushing = false;
        brush.setAttribute("visibility", "hidden");
      });
      svg.appendChild(brush);
      svg.appendChild(hit);
    } else if (typeof options.onSelect === "function") {
      hit.addEventListener("click", (ev) => options.onSelect(posAt(ev)));
      svg.appendChild(hit);
    } else {
      svg.appendChild(hit);
    }
  }
  if (sliderGap) {
    const bandBottom = bandRows.length
      ? 20 + (bandRows.length - 1) * (BOX_H + ROW_GAP) + BOX_H
      : 16;
    // The open strip between the gene boxes and the rate plot. A page that asked for the
    // gap positions its slider from these viewBox units.
    svg.dataset.sliderTop = String(bandBottom + 10);
    svg.dataset.sliderLeft = String(padLeft);
    svg.dataset.sliderWidth = String(plotW);
    svg.dataset.sliderSvgWidth = String(width);
    svg.dataset.sliderSvgHeight = String(height);
  }
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
  // the pie bigger instead of just adding white space to the right of it. `legend: false`
  // drops the side key and centres the circle; Browse puts that key in a table under the pie.
  const showLegend = options.legend !== false;
  const legendW = showLegend ? (options.legendWidth || 215) : 0;
  const r = options.radius ||
            Math.max(80, Math.min((height - 58) / 2,
                                   (width - legendW - (showLegend ? 56 : 48)) / 2));
  const cx = showLegend ? 26 + r : width / 2;
  const cy = Math.max(r + 26, height / 2);
  const fsSlice = r >= 130 ? 13 : r >= 80 ? 11 : 9;
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
  // The dashboard puts the title in the card heading, where it lines up with the chart
  // beside it. Drawing it again inside the SVG left a second, smaller title and a gap.
  if (options.title && options.drawTitle !== false) {
    svg.appendChild(el("text", { x: 0, y: 13, "font-size": 12 }, options.title));
  }

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

  if (showLegend) {
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
      // The rate itself is dropped where the caller asks for a compact legend: it is already
      // on the selector under the pie, and carrying it here costs about 80 px of legend
      // width that the circle can have instead.
      g.appendChild(el("text", { x: lx + 15, y: ly, "font-size": fsLegend, fill: "#1b1f23",
                                 "font-weight": isSel(s) ? 600 : 400 },
        options.compactLegend ? `${s.label}  ${pct}%`
                              : `${s.label}  ${pct}%  (${s.value.toExponential(2)})`));
      svg.appendChild(wire(g, s));
      ly += legendStep;
    });
  }
  container.appendChild(svg);
}

