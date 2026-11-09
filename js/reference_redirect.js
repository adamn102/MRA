// SPDX-License-Identifier: MIT
//
// Forwards the retired download.html and index.html-as-Overview to reference.html. The
// Overview and Download tabs are now one page, "Reference and Download". A module rather than
// a meta refresh, so the query string survives and the page's CSP can stay script-src 'self'.

const qs = window.location.search;
const hash = window.location.hash;
const href = "reference.html" + qs + hash;
document.getElementById("target").href = href;
window.location.replace(href);
