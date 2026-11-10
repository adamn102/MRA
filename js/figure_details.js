// SPDX-License-Identifier: MIT
//
// figure_details.js -- move a figure's explanation behind a "Details" button.
//
// Every chart on this site carried its explanation as a paragraph under it. Read once, that
// prose is what makes the figure interpretable; read every time, it is what made the page
// look like notes rather than a dashboard, and on the position dashboard it put several
// lines of text between two figures meant to be compared.
//
// So the text is kept -- none of it is deleted, because a reader meeting a supply-flow chart
// for the first time genuinely needs it -- and moved into a dialog opened from a small
// button on the figure. A one-line label stays visible, so a figure still says what it is at
// a glance and a downloaded PNG is not anonymous.
//
// Native <dialog>: it gives focus trapping, Escape-to-close and the backdrop for free, and it
// is in every browser this portal supports. Where showModal is missing (jsdom, and any very
// old engine) the dialog is opened non-modally instead, so the content is still reachable.

/**
 * The one button row a plot carries, created on first use. Both the PNG download and the
 * Details button go in it: two modules each making their own row left two stacked strips of
 * buttons under every figure, which is most of what made the page look unfinished.
 */
function escapeHtml(s) {
  return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

export function toolRow(host) {
  let row = host.querySelector(":scope > .plot-tools");
  if (!row) {
    row = (host.ownerDocument || document).createElement("div");
    row.className = "plot-tools";
    row.setAttribute("data-export-skip", "");   // never part of the figure's own PNG
    host.appendChild(row);
  }
  return row;
}

/**
 * Give `host` a Details button and a dialog, or update the ones it already has.
 *
 * Idempotent: the pages redraw a figure on every selection, and re-attaching must not leave
 * a row of buttons behind.
 *
 * @param {Element} host   the element the button is appended to (a plot card or a figure)
 * @param {string}  html   the explanation, as HTML
 * @param {object}  [opts] title: dialog heading; label: button text
 * @returns {Element|null} the button, or null if there is nothing to explain
 */
export function setDetails(host, html, opts = {}) {
  if (!host) return null;
  if (!html || !String(html).trim()) {            // nothing to say: remove a stale button
    const stale = host.querySelector(":scope > .plot-tools .details-button");
    if (stale) {
      const d = stale.parentNode.querySelector("dialog.details-dialog");
      if (d) d.remove();
      stale.remove();
    }
    return null;
  }
  const row = toolRow(host);

  let btn = row.querySelector(".details-button");
  let dialog = row.querySelector("dialog");
  if (!btn) {
    btn = document.createElement("button");
    btn.className = "secondary details-button";
    btn.type = "button";
    row.appendChild(btn);

    dialog = document.createElement("dialog");
    dialog.className = "details-dialog";
    const body = document.createElement("div");
    body.className = "details-body";
    const close = document.createElement("button");
    close.className = "secondary details-close";
    close.type = "button";
    close.textContent = "Close";
    close.addEventListener("click", () => dialog.close());
    dialog.append(body, close);
    row.appendChild(dialog);

    btn.addEventListener("click", () => {
      if (typeof dialog.showModal === "function") dialog.showModal();
      else dialog.setAttribute("open", "");       // jsdom and very old engines
    });
    // Clicking the backdrop closes it. The dialog fills its own box, so a click whose
    // coordinates fall outside that box landed on the backdrop.
    dialog.addEventListener("click", (e) => {
      if (e.target !== dialog) return;
      const r = dialog.getBoundingClientRect();
      const out = e.clientX < r.left || e.clientX > r.right ||
                  e.clientY < r.top || e.clientY > r.bottom;
      if (out) dialog.close();
    });
  }

  btn.textContent = opts.label || "Details";
  // The heading is inserted as HTML. Titles now include a submitted sequence name, so
  // a header that contains < or & must not be read as markup.
  const heading = opts.title ? `<h3>${escapeHtml(opts.title)}</h3>` : "";
  dialog.querySelector(".details-body").innerHTML = heading + html;
  dialog.setAttribute("aria-label", opts.title || "Figure details");
  return btn;
}
