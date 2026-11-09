// SPDX-License-Identifier: MIT
//
// Navigation after the tab restructure. Run:
//   node test/nav_cases.mjs --version v0.1
//
// The MRA Table dropdown makes a promise -- nine named sections, each reachable -- that is
// spread across three places: TABLE_SECTIONS in page.js, the anchors the page renders, and
// the catalogs the sections are built from. Nothing stops those drifting apart, and when they
// do the symptom is a menu entry that scrolls nowhere, which no syntax check or unit test
// sees. This suite ties them together, and checks that every retired address still forwards.
//
// Needs jsdom for the rendered-anchor checks; those SKIP without it, as in page_boot_cases.

import fs from "node:fs";
import path from "node:path";

const version = process.argv.includes("--version")
  ? process.argv[process.argv.indexOf("--version") + 1] : "v0.1";
const REPO = process.cwd();
let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

const { TABLE_SECTIONS } = await import(path.join(REPO, "js/page.js"));

// --- the menu the user asked for ---------------------------------------------------------
const WANTED = [
  ["pi",             "Protease Inhibitors (PI)"],
  ["nrti",           "Nucleoside RT Inhibitors (NRTI)"],
  ["nnrti",          "Non-nucleoside RT Inhibitors (NNRTI)"],
  ["ini",            "Integrase strand-transfer Inhibitors (INI)"],
  ["cai",            "Capsid Inhibitors (CAI)"],
  ["bnab-cd4bs",     "bNAb \u2013 CD4bs"],
  ["bnab-v3-glycan", "bNAb \u2013 V3 Glycan"],
  ["bnab-v2-apex",   "bNAb \u2013 V2 Apex"],
  ["custom",         "Custom"],
];

console.log("The dropdown offers exactly the requested sections, in order");
check("nine entries", TABLE_SECTIONS.length === 9, String(TABLE_SECTIONS.length));
for (let i = 0; i < WANTED.length; i++) {
  const got = TABLE_SECTIONS[i];
  check(`entry ${i + 1} is ${WANTED[i][1]}`,
        got && got.id === WANTED[i][0] && got.label === WANTED[i][1],
        got ? `${got.id} / ${got.label}` : "missing");
}
check("section ids are unique",
      new Set(TABLE_SECTIONS.map((s) => s.id)).size === TABLE_SECTIONS.length);

// --- every section names a class the shipped catalogs actually carry ----------------------
console.log("Each section is backed by data");
const taxonomy = JSON.parse(
  fs.readFileSync(path.join(REPO, "data", version, "drm_taxonomy.json"), "utf8"));
const drmClasses = new Set(taxonomy.classes.map((c) => c.class));
const bnabCsv = fs.readFileSync(path.join(REPO, "data", version, "bnab_catalog.csv"), "utf8");
const bnabHeader = bnabCsv.split(/\r?\n/)[0].split(",");
const epitopeCol = bnabHeader.indexOf("epitope_class");
const bnabClasses = new Set(bnabCsv.split(/\r?\n/).slice(1).filter(Boolean)
  .map((l) => l.split(",")[epitopeCol]));

for (const s of TABLE_SECTIONS) {
  if (s.kind === "drm") {
    check(`${s.id}: ${s.key} is a drug class in drm_taxonomy.json`, drmClasses.has(s.key));
  } else if (s.kind === "bnab") {
    check(`${s.id}: ${s.key} is an epitope class in bnab_catalog.csv`, bnabClasses.has(s.key));
  }
}
check("every drug class in the catalog has a section",
      [...drmClasses].every((c) => TABLE_SECTIONS.some((s) => s.key === c)),
      [...drmClasses].join(", "));
check("every bNAb epitope class in the catalog has a section",
      [...bnabClasses].every((c) => TABLE_SECTIONS.some((s) => s.key === c)),
      [...bnabClasses].join(", "));

// --- the retired addresses still forward --------------------------------------------------
console.log("Retired addresses forward to pages that exist");
const STUBS = [
  ["browse.html", "js/browse_redirect.js", "index.html"],
  ["drm.html", "js/drm_redirect.js", "tables.html"],
  ["download.html", "js/reference_redirect.js", "reference.html"],
  ["compare.html", "js/compare_redirect.js", "sequence.html"],
  ["explore.html", "js/explore_redirect.js", "index.html"],
];
for (const [stub, module, target] of STUBS) {
  check(`${stub} exists`, fs.existsSync(path.join(REPO, stub)));
  const src = fs.readFileSync(path.join(REPO, module), "utf8");
  check(`${stub} forwards to ${target}`, src.includes(`"${target}"`), module);
  check(`${target} exists`, fs.existsSync(path.join(REPO, target)));
  const html = fs.readFileSync(path.join(REPO, stub), "utf8");
  check(`${stub} declares #target for the manual link`, /id="target"/.test(html));
}
check("no stub forwards to another stub",
      !STUBS.some(([, , target]) => STUBS.some(([s]) => s === target)),
      STUBS.map(([, , t]) => t).join(", "));

// --- the page really renders the anchors the menu names ------------------------------------
let JSDOM;
for (const spec of [process.env.JSDOM_PATH, "jsdom",
                    path.join(REPO, "node_modules/jsdom/lib/api.js"),
                    "/tmp/domtest/node_modules/jsdom/lib/api.js"]) {
  if (!spec) continue;
  try { ({ JSDOM } = await import(spec)); break; } catch { /* try the next location */ }
}
if (!JSDOM) {
  console.log("jsdom not found - SKIPPING rendered-anchor checks.");
} else {
  const html = fs.readFileSync(path.join(REPO, "tables.html"), "utf8");
  const dom = new JSDOM(html, { url: `https://example.org/tables.html?v=${version.replace(/^v/, "")}`,
                                pretendToBeVisual: true });
  const win = dom.window;
  globalThis.window = win;
  globalThis.document = win.document;
  for (const g of ["Node", "Element", "HTMLElement", "Event", "CustomEvent", "Blob", "URL",
                   "FileReader", "File", "DOMParser", "XMLSerializer", "SVGElement",
                   "getComputedStyle", "requestAnimationFrame", "cancelAnimationFrame"]) {
    if (win[g] !== undefined) globalThis[g] = win[g];
  }
  globalThis.sessionStorage = { getItem: () => null, setItem: () => {}, removeItem: () => {} };
  globalThis.localStorage = globalThis.sessionStorage;
  globalThis.fetch = async (u) => {
    const rel = String(u).replace(/^https?:\/\/[^/]+\//, "").split("?")[0];
    const p = path.join(REPO, decodeURIComponent(rel));
    if (!fs.existsSync(p)) {
      return { ok: false, status: 404, statusText: "not found",
               json: async () => { throw new Error("404 " + rel); },
               text: async () => { throw new Error("404 " + rel); } };
    }
    const buf = fs.readFileSync(p);
    return { ok: true, status: 200,
             json: async () => JSON.parse(buf.toString("utf8")),
             text: async () => buf.toString("utf8"),
             arrayBuffer: async () =>
               buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength) };
  };

  await import(`${REPO}/js/tables.js?nav=${Date.now()}`);
  const doc = win.document;

  console.log("Every menu entry lands on a section that exists");
  for (const s of TABLE_SECTIONS) {
    const el = doc.getElementById(s.id);
    check(`#${s.id} is rendered`, !!el, el ? el.tagName : "missing");
  }
  check("all 5 DRM sections carry a table",
        TABLE_SECTIONS.filter((s) => s.kind === "drm")
          .every((s) => doc.querySelector(`#${s.id} table`)),
        String(TABLE_SECTIONS.filter((s) => s.kind === "drm").length));
  check("all 3 bNAb sections carry a table",
        TABLE_SECTIONS.filter((s) => s.kind === "bnab")
          .every((s) => doc.querySelector(`#${s.id} table`)));
  check("the custom section comes last in document order",
        [...doc.querySelectorAll("section.table-section")].pop()?.id === "custom");
  check("each DRM section has its own stratum control, not one shared page-wide",
        TABLE_SECTIONS.filter((s) => s.kind === "drm")
          .every((s) => doc.querySelectorAll(`#${s.id} select`).length >= 2));

  console.log("The nav dropdown is built and matches the index on the page");
  const menuLinks = [...doc.querySelectorAll("nav.site .nav-menu a")];
  check("the nav renders a 9-item dropdown", menuLinks.length === 9, String(menuLinks.length));
  check("every dropdown link points at tables.html with a hash",
        menuLinks.every((a) => /tables\.html\?v=[^#]*#/.test(a.getAttribute("href"))),
        menuLinks[0]?.getAttribute("href") || "");
  check("every dropdown hash names a rendered section",
        menuLinks.every((a) => doc.getElementById(a.getAttribute("href").split("#")[1])));
  const indexLinks = [...doc.querySelectorAll("#section-index a")];
  check("the on-page index lists the same nine sections",
        indexLinks.length === 9, String(indexLinks.length));
  check("the index and the dropdown agree on order",
        indexLinks.map((a) => a.getAttribute("href")).join() ===
        TABLE_SECTIONS.map((s) => "#" + s.id).join());

  console.log("The dropdown is reachable without a mouse");
  const group = doc.querySelector("nav.site .nav-group");
  check("the trigger is still a real link to the whole page",
        !!group && group.querySelector("a")?.getAttribute("href")?.startsWith("tables.html"));
  const css = fs.readFileSync(path.join(REPO, "css/mra.css"), "utf8");
  check("the menu opens on hover", /\.nav-group:hover \.nav-menu/.test(css));
  check("the menu opens on keyboard focus", /\.nav-group:focus-within \.nav-menu/.test(css));
}

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL NAV CASES PASSED");
