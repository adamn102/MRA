// SPDX-License-Identifier: MIT
//
// sw.js -- make the portal work offline after the first load.
//
// A clinic or a conference room is exactly where someone will want this tool and exactly
// where the network will be unreliable, so "works offline after first load" is a
// requirement rather than a nicety.
//
// CACHING STRATEGY, and why it is split
//   data/<version>/...   CACHE FIRST, forever -- for the MODEL OUTPUTS only. Those are
//                        immutable by design: nothing is overwritten in place, a new
//                        release is a new directory, so a cached copy can never be stale.
//   the catalog overlays NETWORK FIRST, HTTP cache bypassed. bnab_catalog.csv and
//                        drm_catalog.csv live under a version directory but are curated
//                        layers that get revised, so they are explicitly NOT frozen
//                        (see MUTABLE_OVERLAYS). Editing one in place and treating the
//                        directory as immutable is what served a stale table before.
//   data/releases.json   NETWORK FIRST. It is the one file that changes when a new
//                        release is cut, so a cached copy would hide the new version.
//   everything else      NETWORK FIRST, cache as fallback. Page and code updates should
//                        arrive promptly, but an offline visit must still work.
//
// This worker never caches anything a user submits: submitted sequences are not fetched
// over the network at all, so they never pass through here.

// Bumped v1 -> v2 when the bNAb/DRM catalogs became mutable overlays: the activate handler
// deletes every cache whose key is not this one, so the bump purges the old frozen copies that
// were pinned under the previous cache-first-forever rule.
// Bumped when the shell changes. v6: the tab restructure -- Browse became the landing page,
// Overview and Download merged into reference.html, Table analysis and Curated tables merged
// into tables.html. Without a bump a returning reader keeps yesterday's nav from cache.
const CACHE = "mra-portal-v6";

// Curated annotation overlays: under a version directory, but revised over time rather than
// frozen. Served network-first with the HTTP cache bypassed, so an edit is never masked by an
// old "immutable" copy. Keep this in step with the /data/v*/*.csv rules in _headers.
const MUTABLE_OVERLAYS = /\/(bnab_catalog|drm_catalog)\.csv$/;

const SHELL = [
  // The five live pages, then the redirect stubs for addresses that have moved. A stub is
  // cheap and must stay cached: an offline reader following an old bookmark should still be
  // forwarded rather than meeting a blank page.
  "./", "./index.html", "./tables.html", "./sequence.html", "./batch.html", "./reference.html",
  "./browse.html", "./compare.html", "./download.html", "./drm.html", "./explore.html",
  "./css/mra.css", "./vendor/html2canvas/html2canvas.min.js",
  "./js/version.js", "./js/page.js", "./js/data.js", "./js/table.js", "./js/chart.js", "./js/track_view.js", "./js/explore_core.js",
  "./js/mra_engine.js", "./js/align.js", "./js/align.worker.js", "./js/analysis.js",
  "./js/mutations.js", "./js/browse.js", "./js/tables.js", "./js/sequence.js",
  "./js/curated_tables.js", "./js/tables_custom.js", "./js/reference.js",
  "./js/figure_export.js",
  "./js/explore_redirect.js", "./js/compare_redirect.js", "./js/browse_redirect.js",
  "./js/drm_redirect.js", "./js/reference_redirect.js",
  "./js/batch.js", "./js/batch_core.js", "./js/batch_view.js",
  "./js/codon_sensitivity.js",
  "./data/releases.json",
];

self.addEventListener("install", (event) => {
  event.waitUntil((async () => {
    const cache = await caches.open(CACHE);
    // addAll fails the whole install if any single file 404s, which would leave the
    // portal with no offline support and no clue why. Add them individually instead.
    await Promise.all(SHELL.map(async (url) => {
      try { await cache.add(new Request(url, { cache: "reload" })); }
      catch (err) { console.warn("sw: could not pre-cache", url, err.message); }
    }));
    self.skipWaiting();
  })());
});

self.addEventListener("activate", (event) => {
  event.waitUntil((async () => {
    for (const key of await caches.keys()) {
      if (key !== CACHE) await caches.delete(key);
    }
    await self.clients.claim();
  })());
});

self.addEventListener("fetch", (event) => {
  const req = event.request;
  if (req.method !== "GET") return;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return;     // nothing off-origin, ever

  const isVersionedAsset = /\/data\/v[0-9.]+\//.test(url.pathname);
  const isReleases = url.pathname.endsWith("/data/releases.json");
  const isMutableOverlay = MUTABLE_OVERLAYS.test(url.pathname);

  // Frozen model outputs under a release directory: cache first, forever. The mutable overlays
  // are excluded here and fall through to the network-first branch below.
  if (isVersionedAsset && !isReleases && !isMutableOverlay) {
    event.respondWith((async () => {
      const cache = await caches.open(CACHE);
      const hit = await cache.match(req);
      if (hit) return hit;
      const resp = await fetch(req);
      if (resp.ok) cache.put(req, resp.clone());
      return resp;
    })());
    return;
  }

  event.respondWith((async () => {
    const cache = await caches.open(CACHE);
    try {
      // Overlays bypass the HTTP cache (cache: "reload"), so a browser that stored the old copy
      // as "immutable" before this change is not pinned to it; everything else is a normal fetch.
      const resp = isMutableOverlay
        ? await fetch(url.href, { cache: "reload" })
        : await fetch(req);
      if (resp.ok) cache.put(req, resp.clone());
      return resp;
    } catch (err) {
      const hit = await cache.match(req) || await cache.match(new URL("./index.html", self.location).href);
      if (hit) return hit;
      throw err;
    }
  })());
});
