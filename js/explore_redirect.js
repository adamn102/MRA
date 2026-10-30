// SPDX-License-Identifier: MIT
//
// Forwards the retired explore.html to browse.html, carrying its query string so a shared
// link lands on the position it named. Written as a module rather than a meta refresh so the
// parameters survive and the page's CSP can stay script-src 'self'.

const qs = window.location.search;
const href = "browse.html" + qs;
document.getElementById("target").href = href;
window.location.replace(href);
