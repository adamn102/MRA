// SPDX-License-Identifier: MIT
//
// Forwards the retired compare.html to sequence.html, carrying its query string so a shared
// link still lands sensibly. A module rather than a meta refresh, so parameters survive and the
// page's CSP can stay script-src 'self'. The comparison views are now part of Your sequence.

const qs = window.location.search;
const href = "sequence.html" + qs;
document.getElementById("target").href = href;
window.location.replace(href);
