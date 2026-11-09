// SPDX-License-Identifier: MIT
//
// Forwards the retired drm.html to the MRA Table page. The curated tables are now sections of
// that page rather than a tab of their own, so a link naming a drug class is turned into the
// matching anchor: drm.html?class=NNRTI lands on tables.html#nnrti. A module rather than a
// meta refresh, so the parameters survive and the page's CSP can stay script-src 'self'.

const qs = new URLSearchParams(window.location.search);
const cls = qs.get("class");
// The anchor ids are the lowercased class codes; see TABLE_SECTIONS in page.js.
const hash = cls ? "#" + cls.toLowerCase() : window.location.hash;
const href = "tables.html" + window.location.search + hash;
document.getElementById("target").href = href;
window.location.replace(href);
