(function () {
  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#39;",
    }[c]));
  }

  function getUrlRoot() {
    // Sphinx provides this on built pages (including ReadTheDocs)
    if (window.DOCUMENTATION_OPTIONS && window.DOCUMENTATION_OPTIONS.URL_ROOT) {
      return window.DOCUMENTATION_OPTIONS.URL_ROOT; // e.g. "../../" on RTD
    }
    return "";
  }

  async function loadIndex() {
    const base = getUrlRoot();
    const url = base + "_static/api_index.json";
    const res = await fetch(url);
    if (!res.ok) throw new Error("Could not load api_index.json from " + url);
    return await res.json(); // array of {name, qual, url}
  }

  function pickBestExact(exactList) {
    // Prefer the canonical ProDy parsePDB if it exists
    const preferred = exactList.find(
      (e) => e.qual === "prody.proteins.pdbfile.parsePDB"
    );
    if (preferred) return preferred;

    // Otherwise: choose shortest qualified name (usually the most "core")
    return exactList.slice().sort((a, b) => (a.qual || "").length - (b.qual || "").length)[0];
  }

  function buildUI() {
    const form = document.getElementById("smartSearchForm");
    const input = document.getElementById("smartSearchInput");
    const status = document.getElementById("smartSearchStatus");

    const resultsBox = document.createElement("div");
    resultsBox.id = "smartSearchResults";
    resultsBox.style.marginTop = "12px";
    form.parentElement.appendChild(resultsBox);

    function clearResults() {
      resultsBox.innerHTML = "";
    }

    function showResults(items, query) {
      clearResults();
      if (!items.length) {
        resultsBox.innerHTML = `No matches for <b>${escapeHtml(query)}</b>.`;
        return;
      }
      const ul = document.createElement("ul");
      ul.style.paddingLeft = "18px";

      items.slice(0, 50).forEach((it) => {
        const li = document.createElement("li");
        const a = document.createElement("a");
        a.href = it.url;
        a.textContent = it.qual; // show full path
        li.appendChild(a);
        ul.appendChild(li);
      });

      resultsBox.appendChild(ul);
    }

    return { form, input, status, clearResults, showResults };
  }

  (async function main() {
    const ui = buildUI();

    let index = [];
    try {
      index = await loadIndex();
      ui.status.textContent =
        "Enter an exact function name to jump directly. Partial text shows suggestions.";
    } catch (e) {
      ui.status.textContent =
        "API index not found (_static/api_index.json). Smart Search cannot work.";
      return;
    }

    // Map short-name (lowercased) -> list of entries
    const byLowerName = new Map();
    for (const it of index) {
      const key = (it.name || "").toLowerCase();
      if (!byLowerName.has(key)) byLowerName.set(key, []);
      byLowerName.get(key).push(it);
    }

    function getPartialMatches(qLower) {
      const matches = index.filter((it) =>
        (it.name || "").toLowerCase().includes(qLower)
      );

      // startsWith first, then alphabetical
      matches.sort((a, b) => {
        const an = (a.name || "").toLowerCase();
        const bn = (b.name || "").toLowerCase();
        const aStarts = an.startsWith(qLower) ? 0 : 1;
        const bStarts = bn.startsWith(qLower) ? 0 : 1;
        if (aStarts !== bStarts) return aStarts - bStarts;
        return an.localeCompare(bn);
      });

      return matches;
    }

    // ✅ ENTER / Search button: exact match redirects, else show suggestions
    ui.form.addEventListener("submit", function (ev) {
      ev.preventDefault();
      const q = (ui.input.value || "").trim();
      ui.clearResults();
      if (!q) return;

      const qLower = q.toLowerCase();

      const exact = byLowerName.get(qLower);
      if (exact && exact.length) {
        const best = pickBestExact(exact);
        window.location.href = best.url;
        return;
      }

      const matches = getPartialMatches(qLower);
      ui.showResults(matches, q);
    });

    // ✅ Typing: ONLY show suggestions (do NOT redirect)
    ui.input.addEventListener("input", function () {
      const q = (ui.input.value || "").trim();
      ui.clearResults();
      if (!q) return;
      const qLower = q.toLowerCase();

      // If exact exists, we don't redirect while typing — we just show it first
      const exact = byLowerName.get(qLower) || [];
      const partial = getPartialMatches(qLower);

      // Put exact matches (if any) at top, then partial
      const seen = new Set();
      const merged = [];
      for (const it of exact.concat(partial)) {
        const key = it.qual + "||" + it.url;
        if (seen.has(key)) continue;
        seen.add(key);
        merged.push(it);
      }

      ui.showResults(merged, q);
    });
  })();
})();
