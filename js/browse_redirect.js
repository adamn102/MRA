// SPDX-License-Identifier: MIT
//
// Forwards the old browse.html to index.html. Browsing the atlas is now the landing page, so
// the view did not move -- only its address did. A module rather than a meta refresh, so the
// query string survives (shared links carry ?v= and a position) and the page's CSP can stay
// script-src 'self'.

const qs = window.location.search;
const hash = window.location.hash;
const href = "index.html" + qs + hash;
document.getElementById("target").href = href;
window.location.replace(href);
