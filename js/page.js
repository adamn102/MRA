// SPDX-License-Identifier: MIT
//
// page.js -- the shell every page shares: header, nav, version banner and selector,
// footer with licence and citation, and the bootstrap that resolves which atlas
// version this page is showing.
//
// Written once here so that the per-page files contain only what is specific to that
// view. Nothing in this file fetches anything except same-origin assets.

import { loadReleases, requestedVersionFromUrl, resolveVersion, renderVersionBanner,
         renderVersionSelector, citationText, provenanceHeader, shareableUrl } from "./version.js";

const PAGES = [
  { href: "index.html", label: "Overview" },
  { href: "browse.html", label: "Browse the atlas" },
  { href: "sequence.html", label: "Your sequence" },
  { href: "drm.html", label: "Resistance & escape" },
  { href: "explore.html", label: "Explore the genome" },
  { href: "compare.html", label: "Compare & genome track" },
  { href: "download.html", label: "Download" },
];

export const PRIVACY_TEXT =
  "Sequences you paste or upload never leave your computer. This is a static site: " +
  "there is no server, no database, no analytics, no telemetry and no error reporter. " +
  "Alignment and all rate calculations run inside your browser. The one piece of " +
  "third-party code is html2canvas, which turns a chart into a PNG for download; it is " +
  "served from this site rather than a CDN, and makes no network requests of its own. " +
  "You can confirm all of this by opening your browser's network panel \u2014 after the page " +
  "and its data files load, submitting a sequence produces no further requests \u2014 or by " +
  "disconnecting from the network and using the tool offline.";

/** Build the header, nav and version controls; returns the resolved version. */
export async function bootstrapPage(activeHref) {
  const header = document.createElement("header");
  header.className = "site";
  const inner = document.createElement("div");
  inner.className = "inner";
  const brand = document.createElement("span");
  brand.className = "brand";
  brand.textContent = "HIV-1 Mutation Rate Atlas";
  const nav = document.createElement("nav");
  nav.className = "site";
  const slot = document.createElement("div");
  slot.className = "version-slot";
  inner.appendChild(brand); inner.appendChild(nav); inner.appendChild(slot);
  header.appendChild(inner);
  document.body.insertBefore(header, document.body.firstChild);

  const bannerSlot = document.createElement("div");
  bannerSlot.id = "version-banner-slot";
  document.body.insertBefore(bannerSlot, header.nextSibling);

  let releases, resolved;
  try {
    const requested = requestedVersionFromUrl();
    releases = await loadReleases();
    resolved = resolveVersion(releases, requested);
  } catch (err) {
    showFatal(err.message);
    throw err;
  }

  // Every link the UI generates carries an explicit version, never the implicit default.
  const v = resolved.version.replace(/^v/, "");
  for (const p of PAGES) {
    const a = document.createElement("a");
    a.href = p.href + "?v=" + v;
    a.textContent = p.label;
    if (p.href === activeHref) a.setAttribute("aria-current", "page");
    nav.appendChild(a);
  }
  renderVersionSelector(slot, resolved, releases);
  renderVersionBanner(bannerSlot, resolved, releases);
  buildFooter(resolved.release);
  registerOfflineSupport();
  return { ...resolved, releases };
}

// Offline support. Service workers need a secure context, so this is a no-op when the
// folder is opened straight off the filesystem; everything else still works.
function registerOfflineSupport() {
  if (!("serviceWorker" in navigator)) return;
  // Development escape hatch: ?nosw=1 skips registration and tears down any worker
  // already installed for this origin. Without it, an edit to a module can be masked
  // by yesterday's cached copy and look like a code bug.
  if (new URLSearchParams(location.search).has("nosw")) {
    navigator.serviceWorker.getRegistrations()
      .then((rs) => rs.forEach((r) => r.unregister()))
      .catch(() => {});
    return;
  }
  if (location.protocol !== "https:" && location.hostname !== "localhost" &&
      location.hostname !== "127.0.0.1") return;
  navigator.serviceWorker.register("sw.js").catch((err) => {
    console.warn("offline support unavailable:", err.message);
  });
}

function buildFooter(release) {
  const footer = document.createElement("footer");
  footer.className = "site";
  const inner = document.createElement("div");
  inner.className = "inner";
  inner.innerHTML =
    "<p><strong>Licence.</strong> Portal code: MIT. Atlas data and every file this site " +
    "exports: CC-BY-4.0 \u2014 attribution is required. See <code>LICENSE</code> and " +
    "<code>LICENSE-DATA</code> in the repository.</p>" +
    "<p><strong>Cite this atlas.</strong> Every page and every export carries this block:</p>";
  const pre = document.createElement("pre");
  pre.className = "citation";
  pre.textContent = citationText(release);
  inner.appendChild(pre);
  footer.appendChild(inner);
  document.body.appendChild(footer);
}

/** A privacy panel. Call on any page that accepts a sequence. */
export function privacyPanel() {
  const div = document.createElement("div");
  div.className = "privacy";
  div.innerHTML = "<strong>Your sequence stays on your machine.</strong> " + PRIVACY_TEXT;
  return div;
}

export function note(html) {
  const div = document.createElement("div");
  div.className = "note";
  div.innerHTML = html;
  return div;
}

export function caveat(html) {
  const p = document.createElement("p");
  p.className = "caveat";
  p.innerHTML = html;
  return p;
}

export function showError(container, message) {
  const div = document.createElement("div");
  div.className = "error-box";
  div.textContent = message;
  container.innerHTML = "";
  container.appendChild(div);
}

export function showFatal(message) {
  const main = document.querySelector("main") || document.body;
  const div = document.createElement("div");
  div.className = "error-box";
  div.innerHTML = "<strong>The portal could not start.</strong> " + message;
  main.insertBefore(div, main.firstChild);
}

/** Provenance header for an exported file: version, model, retrieval date. */
export function exportProvenance(release, extraLines = []) {
  let text = provenanceHeader(release, "# ");
  for (const line of extraLines) text += "# " + line + "\n";
  return text;
}

/** A "copy link to this view" button. The link always names its version explicitly. */
export function shareButton(params) {
  const btn = document.createElement("button");
  btn.className = "secondary";
  btn.textContent = "Copy link to this view";
  btn.addEventListener("click", async () => {
    const url = shareableUrl(params.version, params.extra || {});
    try {
      await navigator.clipboard.writeText(url);
      btn.textContent = "Link copied";
    } catch {
      // Clipboard access can be refused; show the URL so the user can copy it by hand.
      btn.textContent = url;
    }
    setTimeout(() => { btn.textContent = "Copy link to this view"; }, 2500);
  });
  return btn;
}

export const NEF_NOTE =
  "<strong>nef is only partly covered.</strong> The model covers NL4-3 reference " +
  "positions 1\u20138,213, which is 72 of roughly 206 nef codons. A query for a nef " +
  "position above 72 returns nothing because the model has no coverage there, not " +
  "because of an error.";

export const CUSTOM_CAVEAT =
  "<strong>This is not a measurement of your virus.</strong> It is the MRA model, " +
  "fitted on NL4-3, applied to your sequence. What genuinely differs between sequences " +
  "is <em>codon identity</em> \u2014 two viruses with the same amino acid at a position but " +
  "different synonymous codons reach different amino acids at different rates \u2014 and, " +
  "second, local sequence context. Position-specific effects were measured on NL4-3 and " +
  "transfer only where your sequence aligns to it. Every row below records which of " +
  "those three components it rests on.";
