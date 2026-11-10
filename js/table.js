// SPDX-License-Identifier: MIT
//
// table.js -- a sortable, filterable, downloadable HTML table. No dependencies.
//
// Deliberately small and obvious: build a <table>, sort on header click, filter on a
// text box, export what is currently shown as CSV with a provenance header. That is all
// any view in this portal needs from a "table library".

export function formatRate(x) {
  if (x === null || x === undefined || Number.isNaN(x)) return "\u2014";
  if (x === 0) return "0";
  return x.toExponential(2).replace("e", "e");
}

export function formatInterval(lo, hi) {
  if (lo === null || hi === null) return "\u2014";
  return formatRate(lo) + " \u2013 " + formatRate(hi);
}

/**
 * @param {HTMLElement} container
 * @param {object} spec
 *   columns: [{key, label, type: "text"|"num"|"rate", title?, render?}]
 *   rows:    array of records
 *   sort:    {key, dir} initial sort
 *   filterPlaceholder, emptyMessage
 *   downloadName, provenance   -- when both given, a Download CSV button appears
 *   csvHeader: "label" writes the visible column headings into the CSV; otherwise the
 *   field names are written
 */
export function renderTable(container, spec) {
  const state = {
    sortKey: spec.sort ? spec.sort.key : null,
    sortDir: spec.sort ? spec.sort.dir : "asc",
    filter: "",
  };

  container.innerHTML = "";
  const controls = document.createElement("div");
  controls.className = "controls";
  container.appendChild(controls);

  if (spec.filterPlaceholder !== false) {
    const wrap = document.createElement("div");
    const label = document.createElement("label");
    label.textContent = "Filter";
    const input = document.createElement("input");
    input.type = "text";
    input.placeholder = spec.filterPlaceholder || "type to filter any column";
    input.size = 28;
    input.addEventListener("input", () => { state.filter = input.value.toLowerCase(); draw(); });
    wrap.appendChild(label); wrap.appendChild(input);
    controls.appendChild(wrap);
  }

  const countEl = document.createElement("div");
  countEl.className = "caveat";
  controls.appendChild(countEl);

  if (spec.downloadName) {
    const wrap = document.createElement("div");
    const btn = document.createElement("button");
    btn.className = "secondary";
    btn.textContent = "Download CSV";
    btn.addEventListener("click", () => downloadCsv(visibleRows(), spec));
    wrap.appendChild(btn);
    controls.appendChild(wrap);
  }

  const scroll = document.createElement("div");
  scroll.className = "table-scroll";
  const table = document.createElement("table");
  table.className = "data";
  const thead = document.createElement("thead");
  const tbody = document.createElement("tbody");
  table.appendChild(thead); table.appendChild(tbody);
  scroll.appendChild(table);
  container.appendChild(scroll);

  function visibleRows() {
    let rows = spec.rows;
    if (state.filter) {
      rows = rows.filter((r) => spec.columns.some((c) => {
        const v = r[c.key];
        return v !== null && v !== undefined && String(v).toLowerCase().includes(state.filter);
      }));
    }
    if (state.sortKey) {
      const col = spec.columns.find((c) => c.key === state.sortKey);
      const dir = state.sortDir === "asc" ? 1 : -1;
      rows = rows.slice().sort((a, b) => {
        const x = a[state.sortKey], y = b[state.sortKey];
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        if (col && (col.type === "num" || col.type === "rate")) return dir * (x - y);
        return dir * String(x).localeCompare(String(y));
      });
    }
    return rows;
  }

  function draw() {
    thead.innerHTML = "";
    const tr = document.createElement("tr");
    for (const c of spec.columns) {
      const th = document.createElement("th");
      th.textContent = c.label + (state.sortKey === c.key ? (state.sortDir === "asc" ? " \u25b2" : " \u25bc") : "");
      th.setAttribute("scope", "col");
      if (c.title) th.title = c.title;
      th.addEventListener("click", () => {
        if (state.sortKey === c.key) state.sortDir = state.sortDir === "asc" ? "desc" : "asc";
        else { state.sortKey = c.key; state.sortDir = c.type === "num" || c.type === "rate" ? "desc" : "asc"; }
        draw();
      });
      tr.appendChild(th);
    }
    thead.appendChild(tr);

    const rows = visibleRows();
    countEl.textContent = rows.length === spec.rows.length
      ? `${rows.length.toLocaleString()} rows`
      : `${rows.length.toLocaleString()} of ${spec.rows.length.toLocaleString()} rows`;

    tbody.innerHTML = "";
    if (!rows.length) {
      const td = document.createElement("td");
      td.colSpan = spec.columns.length;
      td.textContent = spec.emptyMessage || "No rows match.";
      const r = document.createElement("tr"); r.appendChild(td); tbody.appendChild(r);
      return;
    }
    // Cap the DOM at 2,000 rows; the full set still goes into a download. Rendering
    // 62,518 <tr> elements locks the page up for seconds on a laptop.
    const CAP = 2000;
    for (const row of rows.slice(0, CAP)) {
      const tr2 = document.createElement("tr");
      for (const c of spec.columns) {
        const td = document.createElement("td");
        if (c.render) {
          const out = c.render(row);
          if (out instanceof Node) td.appendChild(out); else td.innerHTML = out;
        } else if (c.type === "rate") {
          td.textContent = formatRate(row[c.key]); td.className = "num";
        } else if (c.type === "num") {
          td.textContent = row[c.key] === null ? "\u2014" : row[c.key]; td.className = "num";
        } else {
          td.textContent = row[c.key] === null || row[c.key] === undefined ? "\u2014" : row[c.key];
        }
        tr2.appendChild(td);
      }
      tbody.appendChild(tr2);
    }
    if (rows.length > CAP) {
      const td = document.createElement("td");
      td.colSpan = spec.columns.length;
      td.className = "caveat";
      td.textContent = `Showing the first ${CAP.toLocaleString()} of ${rows.length.toLocaleString()} rows. ` +
                       `Download the CSV for the complete set.`;
      const r = document.createElement("tr"); r.appendChild(td); tbody.appendChild(r);
    }
  }

  draw();
  return { redraw: draw, visibleRows };
}

/** Build a CSV string with a provenance header comment and hand it to the browser. */
export function toCsv(rows, columns, provenance, opts = {}) {
  const esc = (v) => {
    if (v === null || v === undefined) return "";
    const s = String(v);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  };
  const lines = [];
  if (provenance) lines.push(provenance.trimEnd());
  // A table can ask for its visible headings in the file. The default stays the field
  // name, which is what the other downloads have always written.
  const header = opts.header === "label" ? columns.map((c) => c.label) : columns.map((c) => c.key);
  lines.push(header.map(esc).join(","));
  for (const r of rows) lines.push(columns.map((c) => esc(r[c.key])).join(","));
  return lines.join("\n") + "\n";
}

export function downloadCsv(rows, spec) {
  const text = toCsv(rows, spec.columns, spec.provenance, { header: spec.csvHeader });
  downloadText(spec.downloadName, text);
}

export function downloadText(filename, text, mime = "text/csv") {
  // A Blob URL is created and revoked locally. Nothing is uploaded anywhere.
  const blob = new Blob([text], { type: mime + ";charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}
