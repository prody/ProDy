document.addEventListener("DOMContentLoaded", () => {
  const form = document.getElementById("smartSearchForm");
  const input = document.getElementById("smartSearchInput");
  const status = document.getElementById("smartSearchStatus");

  if (!form || !input) return;

  const URL_ROOT =
    (window.DOCUMENTATION_OPTIONS && window.DOCUMENTATION_OPTIONS.URL_ROOT) || "./";

  const apiIndexUrl = URL_ROOT + "_static/api_index.json";

  async function loadIndex() {
    const r = await fetch(apiIndexUrl, { cache: "no-store" });
    if (!r.ok) throw new Error(`Failed to load api_index.json (${r.status})`);
    return await r.json();
  }

  function goToSearch(query) {
    window.location.href = URL_ROOT + "search.html?q=" + encodeURIComponent(query);
  }

  form.addEventListener("submit", async (e) => {
    e.preventDefault();

    const q = (input.value || "").trim();
    if (!q) return;

    if (status) status.textContent = "Searching...";

    try {
      const idx = await loadIndex();

      const key = q.toLowerCase();

      // 1) Exact match by simple name (parsepdb -> parsePDB)
      if (idx.exact_simple && idx.exact_simple[key]) {
        window.location.href = URL_ROOT + idx.exact_simple[key];
        return;
      }

      // 2) Exact match by full name if user pasted it
      if (idx.exact_full && idx.exact_full[key]) {
        window.location.href = URL_ROOT + idx.exact_full[key];
        return;
      }

      // 3) No exact match -> normal Sphinx search (prefix search)
      goToSearch(q);
    } catch (err) {
      console.error(err);
      if (status) status.textContent = "Index not loaded; using normal search.";
      goToSearch(q);
    }
  });

  // Optional: make Enter always submit even if browser is weird
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      form.requestSubmit();
    }
  });
});
