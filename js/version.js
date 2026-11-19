// SPDX-License-Identifier: MIT
//
// version.js -- the version spine of the MRA portal.
//
// THE GOVERNING RULE
//   A link printed in a paper must resolve to the exact atlas version that paper
//   used, forever. A reader following that citation in 2030 must not silently
//   receive a newer model with different numbers.
//
// How that rule is enforced here:
//   1. Assets live under a version directory: data/v0.1/..., data/v1.0/...
//      Nothing is ever overwritten in place. A new release is a new directory.
//   2. data/releases.json (stable path, never versioned) lists every release.
//   3. The version is URL-addressable as ?v=0.1 and defaults to the latest
//      release when absent.
//   4. The model bundle carries its own version string inside the JSON. If the
//      embedded version disagrees with the directory it was fetched from, the
//      portal REFUSES TO START. That mismatch means assets were copied between
//      release directories, which would silently serve the wrong model.
//   5. Every shareable link the UI generates carries an explicit version.
//   6. Viewing anything other than the latest release shows a banner.
//
// No framework. Plain ES module. Read it top to bottom.

export const RELEASES_URL = "data/releases.json";

// ---------------------------------------------------------------------------
// Version strings
// ---------------------------------------------------------------------------
// Canonical form is "v0.1". A user typing ?v=0.1 means the same thing, so we
// normalise. Anything else is rejected rather than guessed at -- a typo in a
// printed URL must fail loudly, not fall back to latest.

export function normaliseVersion(raw) {
  if (raw === null || raw === undefined) return null;
  const s = String(raw).trim().toLowerCase();
  if (s === "") return null;
  const withV = s.startsWith("v") ? s : "v" + s;
  if (!/^v[0-9]+(\.[0-9]+){1,2}$/.test(withV)) {
    throw new Error(
      "Not a valid atlas version string: " + JSON.stringify(raw) +
      '. Expected something like "v0.1".'
    );
  }
  return withV;
}

// Read ?v= from the current URL. Returns null when absent (meaning "latest").
export function requestedVersionFromUrl(search) {
  const qs = new URLSearchParams(search === undefined ? window.location.search : search);
  return normaliseVersion(qs.get("v"));
}

// Build a shareable link that always names its version explicitly.
export function shareableUrl(version, extraParams) {
  const url = new URL(window.location.href);
  url.searchParams.set("v", version.replace(/^v/, ""));
  if (extraParams) {
    for (const [k, val] of Object.entries(extraParams)) url.searchParams.set(k, val);
  }
  return url.toString();
}

// ---------------------------------------------------------------------------
// The release manifest
// ---------------------------------------------------------------------------

export async function loadReleases() {
  const resp = await fetch(RELEASES_URL, { cache: "no-cache" });
  if (!resp.ok) {
    throw new Error("Could not load " + RELEASES_URL + " (HTTP " + resp.status + ")");
  }
  const releases = await resp.json();
  if (!releases.latest || !Array.isArray(releases.releases)) {
    throw new Error("releases.json is malformed: needs 'latest' and 'releases'.");
  }
  return releases;
}

export function findRelease(releases, version) {
  const hit = releases.releases.find((r) => r.version === version);
  if (!hit) {
    const known = releases.releases.map((r) => r.version).join(", ");
    throw new Error(
      "Atlas version " + version + " does not exist. Known versions: " + known + "."
    );
  }
  return hit;
}

// Resolve which version this page is showing, given the URL and the manifest.
export function resolveVersion(releases, requested) {
  if (requested === null) {
    return { version: releases.latest, release: findRelease(releases, releases.latest), isLatest: true };
  }
  const release = findRelease(releases, requested);
  return { version: requested, release: release, isLatest: requested === releases.latest };
}

export function assetUrl(version, filename) {
  return "data/" + version + "/" + filename;
}

// ---------------------------------------------------------------------------
// Loading an asset, with the version cross-check and an optional checksum check
// ---------------------------------------------------------------------------

// Fetch the model bundle and refuse it if its embedded version disagrees with
// the directory it came from. This is rule 4 above and it is a hard stop.
export async function loadModelBundle(version) {
  const url = assetUrl(version, "mra_web_model.json");
  const resp = await fetch(url);
  if (!resp.ok) throw new Error("Could not load model bundle " + url + " (HTTP " + resp.status + ")");
  const bundle = await resp.json();
  if (bundle.version !== version) {
    throw new Error(
      "VERSION MISMATCH -- refusing to start. The bundle at " + url +
      " declares version " + JSON.stringify(bundle.version) +
      " but was served from the " + version + " directory. Assets have been copied " +
      "between release directories; the portal would otherwise serve the wrong model."
    );
  }
  return bundle;
}

// SHA-256 of a fetched asset, compared against the release manifest.
// crypto.subtle is only available in secure contexts (https, or localhost), so
// this silently reports "skipped" when opened straight off the filesystem.
export async function verifyAssetChecksum(version, filename, release) {
  const expected = release.assets && release.assets[filename] && release.assets[filename].sha256;
  if (!expected) return { status: "no-checksum", filename };
  if (!(window.crypto && window.crypto.subtle)) return { status: "skipped", filename };
  const resp = await fetch(assetUrl(version, filename));
  if (!resp.ok) return { status: "fetch-failed", filename };
  const buf = await resp.arrayBuffer();
  const digest = await window.crypto.subtle.digest("SHA-256", buf);
  const hex = Array.from(new Uint8Array(digest))
    .map((b) => b.toString(16).padStart(2, "0"))
    .join("");
  return { status: hex === expected ? "ok" : "MISMATCH", filename, expected, actual: hex };
}

// ---------------------------------------------------------------------------
// The not-latest banner
// ---------------------------------------------------------------------------

export function renderVersionBanner(container, resolved, releases) {
  container.innerHTML = "";
  if (resolved.isLatest) return;
  const latest = findRelease(releases, releases.latest);
  const div = document.createElement("div");
  div.className = "version-banner";
  div.setAttribute("role", "status");
  const changed = latest.description ? " What changed: " + latest.description : "";
  div.innerHTML =
    "<strong>You are viewing atlas " + resolved.version + ", which is not the latest release.</strong> " +
    "The current release is " + releases.latest + " (" + latest.release_date + ")." + changed +
    ' <a href="?v=' + releases.latest.replace(/^v/, "") + '">Go to ' + releases.latest + "</a>. " +
    "This older version is kept permanently so that citations to it keep resolving.";
  container.appendChild(div);
}

// The atlas version, as text. Older releases stay reachable by their ?v= address and are
// archived with the repository; the header does not offer a switcher.
export function renderVersionSelector(container, resolved) {
  container.innerHTML = "";
  const span = document.createElement("span");
  span.textContent = "Atlas version " + resolved.version;
  container.appendChild(span);
}

// ---------------------------------------------------------------------------
// Citation block -- every page and every export carries one
// ---------------------------------------------------------------------------

export function citationText(release, retrievedISO) {
  const doi = release.doi ? "DOI: " + release.doi : "DOI: [not yet assigned]";
  const lines = [
    "HIV-1 Mutation Rate Atlas (MRA), atlas version " + release.version + ".",
    release.paper ? release.paper : "[Paper citation to be completed on publication]",
    doi,
    "Atlas data licensed CC-BY-4.0; portal code licensed MIT.",
    "Retrieved " + (retrievedISO || new Date().toISOString().slice(0, 10)) + " from the MRA portal.",
  ];
  return lines.join("\n");
}

// Prefix for every exported file: who made it, from what, when.
export function provenanceHeader(release, commentChar) {
  const c = commentChar === undefined ? "# " : commentChar;
  return citationText(release)
    .split("\n")
    .map((l) => c + l)
    .join("\n") + "\n";
}
