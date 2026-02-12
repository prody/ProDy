(async function () {
  async function loadIndex() {
    const r = await fetch("_static/api_index.json", { cache: "no-store" });
    if (!r.ok) throw new Error("Cannot load _static/api_index.json");
    return await r.json();
  }

  function normalize(s) {
    return (s || "").trim().toLowerCase();
  }

  function renderResults(container, hits) {
    if (!hits.length) {
      container.innerHTML = "<p>No matches.</p>";
      return;
    }
    container.innerHTML = hits
      .map(
        (h) =>
          `<div style="margin:6px 0;">
             <a href="${h.url}"><b>${h.name}</b></a>
             <div style="opacity:.75;font-size:13px;">${h.full}</div>
           </div>`
      )
      .join("");
  }

  let index = [];
  try {
    index = await loadIndex();
  } catch (e) {
    const out = document.getElementById("smartSearchResults");
    if (out) out.innerHTML = `<p style="color:#b00;">${e.message}</p>`;
    return;
  }

  const form = document.getElementById("smartSearchForm");
  const input = document.getElementById("smartSearchInput");
  const out = document.getElementById("smartSearchResults");

  if (!form || !input || !out) return;

  form.addEventListener("submit", (ev) => {
    ev.preventDefault();

    const q = normalize(input.value);
    if (!q) {
      out.innerHTML = "<p>Type something.</p>";
      return;
    }

    // 1) exact match on function name (case-insensitive)
    const exact = index.find((x) => normalize(x.name) === q);
    if (exact) {
      // go directly to that function doc
      window.location.href = exact.url;
      return;
    }

    // 2) otherwise show partial matches
    const hits = index.filter((x) => {
      const name = normalize(x.name);
      const full = normalize(x.full);
      return name.includes(q) || full.includes(q);
    });

    renderResults(out, hits);
  });
})();
