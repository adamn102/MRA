// SPDX-License-Identifier: MIT
//
// sw.js -- make the portal work offline after the first load.
//
// A clinic or a conference room is exactly where someone will want this tool and exactly
// where the network will be unreliable, so "works offline after first load" is a
// requirement rather than a nicety.
//
// CACHING STRATEGY, and why it is split
//   data/<version>/...   CACHE FIRST, forever. A release directory is immutable by
//                        design: nothing is ever overwritten in place, a new release is
//                        a new directory. So a cached asset can never be stale, and
//                        serving it from cache costs no correctness.
//   data/releases.json   NETWORK FIRST. It is the one file that changes when a new
//                        release is cut, so a cached copy would hide the new version.
//   everything else      NETWORK FIRST, cache as fallback. Page and code updates should
//                        arrive promptly, but an offline visit must still work.
//
// This worker never caches anything a user submits: submitted sequences are not fetched
// over the network at all, so they never pass through here.

const CACHE = "mra-portal-v1";

const SHELL = [
  "./", "./index.html", "./browse.html", "./sequence.html", "./compare.html",
  "./download.html", "./drm.html", "./explore.html", "./css/mra.css", "./vendor/html2canvas/html2canvas.min.js",
  "./js/version.js", "./js/page.js", "./js/data.js", "./js/table.js", "./js/chart.js",
  "./js/mra_engine.js", "./js/align.js", "./js/align.worker.js", "./js/analysis.js",
  "./js/mutations.js", "./js/index.js", "./js/browse.js", "./js/explore.js", "./js/sequence.js",
  "./js/compare.js", "./js/download.js", "./js/drm.js",
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

  if (isVersionedAsset && !isReleases) {
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
      const resp = await fetch(req);
      if (resp.ok) cache.put(req, resp.clone());
      return resp;
    } catch (err) {
      const hit = await cache.match(req) || await cache.match(new URL("./index.html", self.location).href);
      if (hit) return hit;
      throw err;
    }
  })());
});
