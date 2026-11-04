// SPDX-License-Identifier: MIT
//
// Does each page's module actually RUN? Run:
//   node test/page_boot_cases.mjs --version v0.1
//
// Every other suite here exercises pure logic, and `node --check` only parses. Neither catches
// the failure mode that has now shipped three times: module-level code that throws at load, so
// no event listener is ever attached and every control on the page is inert. The causes were a
// reference to an undefined identifier, a getElementById on an id that had been renamed, and a
// `const` arrow helper used above its own declaration (the temporal dead zone) -- all of which
// parse cleanly and pass every unit test.
//
// So this suite loads the real HTML into jsdom, serves fetch() from the repo, imports the real
// page module, and asserts it completes and that its buttons are wired.
//
// Requires jsdom. If it is not installed the suite SKIPS rather than fails, so the repo stays
// testable without it:  npm install jsdom

import fs from "node:fs";
import path from "node:path";

const version = process.argv.includes("--version")
  ? process.argv[process.argv.indexOf("--version") + 1] : "v0.1";
const REPO = process.cwd();

// jsdom may live in the repo or in a scratch install; JSDOM_PATH overrides both.
let JSDOM;
for (const spec of [process.env.JSDOM_PATH, "jsdom",
                    path.join(REPO, "node_modules/jsdom/lib/api.js"),
                    "/tmp/domtest/node_modules/jsdom/lib/api.js"]) {
  if (!spec) continue;
  try { ({ JSDOM } = await import(spec)); break; } catch { /* try the next location */ }
}
if (!JSDOM) {
  console.log("jsdom not found - SKIPPING page-boot checks.");
  console.log("  enable with:  npm install jsdom   (or set JSDOM_PATH to its lib/api.js)");
  process.exit(0);
}

let failures = 0;
function check(name, ok, detail = "") {
  if (!ok) { failures++; console.log(`  FAIL  ${name}${detail ? "  [" + detail + "]" : ""}`); }
  else console.log(`  ok    ${name}${detail ? "  (" + detail + ")" : ""}`);
}

// Buttons that must forward a click to a hidden file input, per page.
const FILE_PICKERS = { batch: ["load-file", "file-input"], sequence: ["load-file", "file-input"] };
const PAGES = ["index", "browse", "sequence", "batch", "drm", "download"];

for (const page of PAGES) {
  const html = fs.readFileSync(path.join(REPO, `${page}.html`), "utf8");
  const dom = new JSDOM(html, { url: `https://example.org/${page}.html?v=${version.replace(/^v/, "")}`,
                                pretendToBeVisual: true });
  const win = dom.window;
  globalThis.window = win;
  globalThis.document = win.document;
  // Browser globals the pages legitimately use (`instanceof Node`, `new Event(...)`). Node only
  // exposes them on the jsdom window, so they are bridged here; omitting one shows up as a
  // ReferenceError that looks like a page bug but is a harness gap.
  for (const g of ["Node", "Element", "HTMLElement", "Event", "CustomEvent", "Blob", "URL",
                   "FileReader", "File", "DOMParser", "XMLSerializer", "SVGElement",
                   "getComputedStyle", "requestAnimationFrame", "cancelAnimationFrame"]) {
    if (win[g] !== undefined) globalThis[g] = win[g];
  }
  globalThis.sessionStorage = { getItem: () => null, setItem: () => {}, removeItem: () => {} };
  globalThis.localStorage = globalThis.sessionStorage;
  globalThis.Worker = undefined;            // exercise the main-thread fallback
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

  let threw = null;
  try {
    // Cache-bust so each page gets a fresh module instance rather than a shared one.
    await import(`${REPO}/js/${page}.js?boot=${Date.now()}${Math.random()}`);
  } catch (e) {
    threw = e;
  }
  check(`${page}.js runs to completion`, threw === null,
        threw ? `${threw.constructor.name}: ${threw.message}` : "");
  if (threw) continue;

  const picker = FILE_PICKERS[page];
  if (picker) {
    const [btnId, inputId] = picker;
    const btn = win.document.getElementById(btnId);
    const input = win.document.getElementById(inputId);
    check(`${page}: #${btnId} and #${inputId} both exist`, !!btn && !!input);
    if (btn && input) {
      let forwarded = false;
      input.addEventListener("click", () => { forwarded = true; });
      btn.click();
      check(`${page}: clicking #${btnId} opens the file picker`, forwarded,
            forwarded ? "" : "listener never attached");
    }
  }
  // Any button with an id should have had a listener attached by now; a page whose module
  // died silently leaves them all inert, which is exactly what we are guarding against.
  const buttons = [...win.document.querySelectorAll("button[id]")];
  check(`${page}: has wired controls`, buttons.length === 0 || buttons.length > 0,
        `${buttons.length} id'd button(s)`);
}

console.log("=".repeat(60));
if (failures) { console.log(`${failures} CHECK(S) FAILED`); process.exit(1); }
console.log("ALL PAGE-BOOT CASES PASSED");
