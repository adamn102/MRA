// SPDX-License-Identifier: MIT
//
// Forwards the retired batch.html to sequence.html, carrying its query string so a shared
// link still lands on the same atlas version. A module rather than a meta refresh, so
// parameters survive and the page's CSP can stay script-src 'self'.

const qs = window.location.search;
const href = "sequence.html" + qs;
document.getElementById("target").href = href;
window.location.replace(href);
