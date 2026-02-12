function runSmartSearch(query) {
  const q = query.trim();
  if (!q) return;

  const base = window.location.pathname.split("/").slice(0, 3).join("/");

  // Redirect to normal RTD search
  window.location.href = base + "/search.html?q=" + encodeURIComponent(q);
}

document.addEventListener("DOMContentLoaded", function () {
  const form = document.getElementById("smartSearchForm");
  const input = document.getElementById("smartSearchInput");

  if (!form || !input) return;

  form.addEventListener("submit", function (e) {
    e.preventDefault();
    runSmartSearch(input.value);
  });

  // THIS makes Enter work
  input.addEventListener("keydown", function (e) {
    if (e.key === "Enter") {
      e.preventDefault();
      runSmartSearch(input.value);
    }
  });
});
