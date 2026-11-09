// SPDX-License-Identifier: MIT
//
// tables.js -- the MRA Table tab.
//
// One page carrying every compiled table: the five drug-resistance classes, the three bNAb
// epitope classes, and last a custom query for a list the reader brings themselves. This
// merges what were two tabs, "Curated tables" and "Table analysis". They were separated on
// the theory that a curated catalog and an ad-hoc lookup are different activities; in
// practice a reader arrives wanting "the rate for these mutations" and should not have to
// know in advance whether their mutations are in the catalog.
//
// This module only assembles. The two halves live in curated_tables.js and tables_custom.js
// and are exported as init functions, because both used to run at module load and claim the
// page for themselves -- putting them on one page means exactly one bootstrapPage call, here.
//
// The section order and anchors come from TABLE_SECTIONS in page.js, which the nav dropdown
// is also built from, so the menu cannot offer a section the page does not render.

import { bootstrapPage, showError, TABLE_SECTIONS } from "./page.js";
import { loadSmallAssets, loadAminoAcidLayer, loadDrmCatalog, loadBnabCatalog } from "./data.js";
import { MraEngine } from "./mra_engine.js";
import { initCuratedTables } from "./curated_tables.js";
import { initCustomQuery } from "./tables_custom.js";

const page = await bootstrapPage("tables.html");

const status = document.getElementById("tables-status");
let assets, engine, drm, bnab, aaLayer;
try {
  assets = await loadSmallAssets(page.version);
  engine = new MraEngine({ bundle: assets.bundle, reference: assets.reference,
                           intervals: assets.intervals, hxb2: assets.hxb2 });
  [drm, bnab, aaLayer] = await Promise.all([
    loadDrmCatalog(page.version),
    loadBnabCatalog(page.version),
    loadAminoAcidLayer(page.version),
  ]);
  status.textContent =
    `${new Set(drm.records.map((r) => r.mutation)).size} DRMs across ` +
    `${drm.taxonomy.classes.length} drug classes, and ` +
    `${bnab.records.length} curated bNAb escape rows, loaded.`;
} catch (err) {
  showError(status, "Could not load the catalogs: " + err.message);
  throw err;
}

const rendered = initCuratedTables({
  page, engine, drm, bnab, host: document.getElementById("table-sections"),
});

// The index mirrors the nav dropdown, but is built from what was actually rendered plus the
// custom section, so it can never advertise a table the catalogs did not supply.
const indexSlot = document.getElementById("section-index");
const ids = new Set([...rendered, "custom"]);
for (const s of TABLE_SECTIONS) {
  if (!ids.has(s.id)) continue;
  const a = document.createElement("a");
  a.href = "#" + s.id;
  a.textContent = s.label;
  indexSlot.appendChild(a);
}

initCustomQuery({ page, engine, aaLayer });

// A link from the retired drm.html arrives as ?class=NNRTI; land on that section. Done after
// both halves render, because the anchor does not exist until its section does.
const qs = new URLSearchParams(window.location.search);
const wanted = qs.get("class");
if (wanted) {
  const match = TABLE_SECTIONS.find(
    (s) => s.key === wanted || s.id === wanted.toLowerCase());
  const target = match && document.getElementById(match.id);
  if (target) target.scrollIntoView();
} else if (window.location.hash) {
  const target = document.getElementById(window.location.hash.slice(1));
  if (target) target.scrollIntoView();
}
