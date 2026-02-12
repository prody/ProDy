async function smartSearch() {
  const q = document.getElementById("smart-search-input").value.trim();
  if (!q) return;

  // where RTD stores search index
  const base = window.location.pathname.split("/").slice(0,3).join("/");
  const indexUrl = base + "/_static/api_index.json";

  try {
    const res = await fetch(indexUrl);
    const data = await res.json();

    // EXACT match only
    if (data[q] && data[q].length > 0) {
      const target = data[q][0];
      const parts = target.split(".");
      const mod = parts.slice(1,-1).join("/");
      window.location.href =
        base + "/reference/" + mod + ".html#" + target;
      return;
    }
  } catch (e) {
    console.warn("Smart search fallback:", e);
  }

  // fallback to normal RTD search
  window.location.href = base + "/search.html?q=" + encodeURIComponent(q);
}
document.addEventListener("DOMContentLoaded", () => {
  const input = document.getElementById("smart-search-input");
  if (!input) return;

  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      smartSearch();
    }
  });
});
