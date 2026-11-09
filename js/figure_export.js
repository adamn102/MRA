// SPDX-License-Identifier: MIT
//
// figure_export.js -- one PNG download button per figure.
//
// Until now each page had a single "Download as PNG" button that captured one fixed region:
// on the batch tab it took the two site charts together, on the browse tab the whole position
// dashboard. The genome-wide track -- the figure most likely to be wanted on its own -- could
// not be saved at all. Anyone who wanted just one panel had to screenshot it by hand.
//
// This module instead treats <figure class="chart"> as the unit of export. Every such figure
// on the page gets its own button, so the set of downloadable plots is exactly the set of
// plots, and a figure added later is covered without touching this file.
//
// WHAT GOES INTO THE IMAGE. A figure's caption carries the settings it was drawn under
// (smoothing window, supply class, which position), so the caption is captured with the plot
// and the PNG is self-describing. Interactive controls are not: a range slider rendered into
// a still image is noise, and html2canvas draws it at the wrong size anyway. Anything marked
// data-export-skip is left out, which is also how the button keeps itself out of its own
// image.
//
// PROVENANCE. A CSV export gets a comment header naming the atlas version and the date; a PNG
// has nowhere to put one, so the same information is drawn into the image as a footer line.
// A figure that ends up in a talk or a manuscript then still says which atlas it came from.

/** Characters that are awkward in a filename on at least one platform. */
const UNSAFE = /[^a-z0-9]+/gi;

/**
 * Build the download filename for a figure.
 *
 * Pure, so the naming rule is testable without a DOM: the caller passes the figure's declared
 * export name (data-export-name, falling back to its id) and the date.
 */
export function exportFileName(name, dateISO) {
  const slug = String(name || "figure").trim().replace(UNSAFE, "_")
    .replace(/^_+|_+$/g, "").toLowerCase() || "figure";
  const day = /^\d{4}-\d{2}-\d{2}$/.test(dateISO || "")
    ? dateISO : new Date().toISOString().slice(0, 10);
  return `mra_${slug}_${day}.png`;
}

/** The one-line provenance footer drawn into the image. */
export function figureStamp(release, dateISO) {
  const day = /^\d{4}-\d{2}-\d{2}$/.test(dateISO || "")
    ? dateISO : new Date().toISOString().slice(0, 10);
  const parts = ["HIV-1 Mutation Rate Atlas"];
  if (release && release.version) parts.push("atlas " + release.version);
  if (release && release.doi) parts.push(release.doi);
  parts.push("retrieved " + day);
  return parts.join(" \u00b7 ");
}

/** Does this figure hold anything that would actually show up in an image? */
export function figureIsDrawn(fig) {
  return !!(fig && fig.querySelector("svg, canvas, img"));
}

/**
 * Capture one <figure> as a PNG and hand it to the browser as a download.
 *
 * Resolves to the filename written. Throws if the figure is empty or html2canvas is missing,
 * so the caller decides how to report it.
 */
export async function downloadFigure(fig, options = {}) {
  if (!figureIsDrawn(fig)) {
    throw new Error("this figure has not been drawn yet");
  }
  if (typeof html2canvas === "undefined") {
    throw new Error("the PNG library did not load; reload the page and try again");
  }
  const name = fig.dataset.exportName || fig.id || "figure";
  const filename = exportFileName(name, options.dateISO);

  // The stamp is part of the image, not part of the page: add it, capture, take it away again
  // in a finally so a failed capture cannot leave it behind.
  const stamp = fig.ownerDocument.createElement("div");
  stamp.className = "figure-stamp";
  stamp.textContent = figureStamp(options.release, options.dateISO);
  fig.appendChild(stamp);
  try {
    const canvas = await html2canvas(fig, {
      backgroundColor: "#ffffff",
      scale: options.scale || 2,
      logging: false,
      useCORS: true,
      ignoreElements: (node) => !!(node.hasAttribute && node.hasAttribute("data-export-skip")),
    });
    const blob = await new Promise((res) => canvas.toBlob(res, "image/png"));
    if (!blob) throw new Error("the canvas could not be converted to an image");
    const url = URL.createObjectURL(blob);
    const a = fig.ownerDocument.createElement("a");
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
    return filename;
  } finally {
    stamp.remove();
  }
}

/**
 * Give one figure its own download button, idempotently.
 *
 * Returns the button, or null if the figure already had one. The button row is marked
 * data-export-skip so it never appears in the image it produces.
 */
export function attachFigureDownload(fig, options = {}) {
  if (!fig || fig.querySelector(":scope > .figure-tools")) return null;
  const doc = fig.ownerDocument;
  const row = doc.createElement("div");
  row.className = "figure-tools";
  row.setAttribute("data-export-skip", "");

  const btn = doc.createElement("button");
  btn.type = "button";
  btn.className = "secondary figure-download";
  btn.textContent = "Download this figure (PNG)";
  const label = btn.textContent;

  btn.addEventListener("click", async () => {
    btn.disabled = true;
    btn.textContent = "Saving\u2026";
    try {
      await downloadFigure(fig, options);
      btn.textContent = "Saved";
    } catch (err) {
      btn.textContent = "Could not save";
      if (typeof options.onError === "function") options.onError(err, fig);
      else console.error("figure export failed:", err);
    } finally {
      // Leave the outcome visible long enough to read, then go back to the offer.
      setTimeout(() => { btn.disabled = false; btn.textContent = label; }, 2000);
    }
  });

  row.appendChild(btn);
  fig.appendChild(row);
  return btn;
}

/**
 * Bring every <figure class="chart"> on the page into line with what it currently holds: a
 * figure offers a download exactly when there is something to download. Call it after a
 * render; it is idempotent, so calling it again costs nothing.
 *
 * Syncing rather than only adding matters because a figure can be emptied again -- selecting
 * a position outside a modelled codon clears the two site charts -- and a download button on
 * a blank figure is an offer the page cannot keep.
 *
 * Returns {attached, removed, drawn} for the caller and for tests.
 */
export function attachFigureDownloads(options = {}) {
  const root = options.root || document;
  let attached = 0, removed = 0, drawn = 0;
  for (const fig of root.querySelectorAll("figure.chart")) {
    const tools = fig.querySelector(":scope > .figure-tools");
    if (figureIsDrawn(fig)) {
      drawn++;
      if (!tools && attachFigureDownload(fig, options)) attached++;
    } else if (tools) {
      tools.remove();
      removed++;
    }
  }
  return { attached, removed, drawn };
}
