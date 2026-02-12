(async function () {
  const form = document.getElementById("smartSearchForm");
  const input = document.getElementById("smartSearchInput");
  const status = document.getElementById("smartSearchStatus");

  // Create a result container (list of links)
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
      resultsBox.innerHTML = `<p>No matches for <b>${escapeHtml(query)}</b>.</p>`;
      return;
    }

    const ul = document.createElement("ul");
    ul.style.paddingLeft = "18px";

    items.slice(0, 50).forEach((it) => {
      const li = document.createElement("li");
      const a = document.createElement("a");
      a.href = it.url;
      a.textContent = it.qual;
      li.appendChild(a);
      ul.appendChild(li);
    });

    resultsBox.appendChild(ul);
  }

  function escapeHtml(s) {
    return s.replace(/[&<>"']/g, (c) => ({
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#039;",
    }[c]));
  }

  async function loadIndex() {
    // Works on RTD and local builds
    const base = document.documentElement.dataset.contentRoot || "";
    const url = base + "_static/api_index.json";

    const res = await fetch(url);
    if (!res.ok) throw new Error("Could not load api_index.json");
    return await res.json(); // array of {name, qual, url}
  }

  let index = [];
  try {
    index = await loadIndex();
    status.textContent = "Exact match jumps directly. Partial searches show suggestions below.";
  } catch (e) {
    status.textContent = "API index not found. Smart Search cannot work.";
    return;
  }

  // Build a fast lookup for exact matches (case-insensitive)
  const byLowerName = new Map(); // lower -> [entries...]
  for (const it of index) {
    const key = (it.name || "").toLowerCase();
    if (!byLowerName.has(key)) byLowerName.set(key, []);
    byLowerName.get(key).push(it);
  }

  function doSearch(qRaw) {
    const q = (qRaw || "").trim();
    clearResults();
    if (!q) return;

    const qLower = q.toLowerCase();

    // ✅ 1) EXACT MATCH (case-insensitive): go directly to the best one
    const exact = byLowerName.get(qLower);
    if (exact && exact.length) {
      // Prefer the canonical parsePDB location if present (optional but nice)
      const preferred = exact.find(e => e.qual === "prody.proteins.pdbfile.parsePDB");
      window.location.href = (preferred || exact[0]).url;
      return;
    }

    // ✅ 2) Otherwise show matches here (NOT Sphinx search results)
    // Rule:
    // - If user typed "parse": show anything containing "parse" (case-insensitive)
    // - Sort: startsWith first, then contains
    const matches = index.filter(it =>
      (it.name || "").toLowerCase().includes(qLower)
    );

    matches.sort((a, b) => {
      const an = a.name.toLowerCase(), bn = b.name.toLowerCase();
      const aStarts = an.startsWith(qLower) ? 0 : 1;
      const bStarts = bn.startsWith(qLower) ? 0 : 1;
      if (aStarts !== bStarts) return aStarts - bStarts;
      return an.localeCompare(bn);
    });

    showResults(matches, q);
  }

  // Enter key works because we handle submit
  form.addEventListener("submit", function (ev) {
    ev.preventDefault();
    doSearch(input.value);
  });

  // Optional: live suggestions while typing
  input.addEventListener("input", function () {
    const q = input.value.trim();
    if (!q) {
      clearResults();
      return;
    }
    // only show suggestions for partial queries, not for exact (exact would redirect on Enter)
    doSearch(q);
  });
})();
