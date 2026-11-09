// SPDX-License-Identifier: MIT
//
// Forwards the retired explore.html to the landing page (index.html), carrying its query
// link lands on the position it named. Written as a module rather than a meta refresh so the
// parameters survive and the page's CSP can stay script-src 'self'.

const qs = window.location.search;
const href = "index.html" + qs;
document.getElementById("target").href = href;
window.location.replace(href);
