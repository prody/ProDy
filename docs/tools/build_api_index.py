import json
import os
import pkgutil
import importlib
import inspect

def iter_modules(pkg, prefix):
    # walk packages but skip tests (they cause side effects / missing files on RTD)
    for m in pkgutil.walk_packages(pkg.__path__, prefix):
        name = m.name
        if name.startswith("prody.tests"):
            continue
        yield name

def main():
    import prody

    out_path = os.path.join(os.path.dirname(__file__), "..", "_static", "api_index.json")
    out_path = os.path.abspath(out_path)

    entries = []
    seen = set()

    for modname in iter_modules(prody, "prody."):
        try:
            mod = importlib.import_module(modname)
        except Exception:
            # skip modules that fail to import on RTD
            continue

        for attr_name, obj in vars(mod).items():
            if attr_name.startswith("_"):
                continue
            if not inspect.isfunction(obj):
                continue
            if getattr(obj, "__module__", None) != modname:
                continue

            full = f"{modname}.{attr_name}"
            if full in seen:
                continue
            seen.add(full)

            # Guess the page path used by your docs:
            # prody.proteins.pdbfile.parsePDB -> reference/proteins/pdbfile.html#prody.proteins.pdbfile.parsePDB
            parts = modname.split(".")
            if len(parts) >= 3:
                section = parts[1]
                page = parts[2]
                url = f"reference/{section}/{page}.html#{full}"
            else:
                url = f"reference/index.html#{full}"

            entries.append({
                "name": attr_name,
                "full": full,
                "url": url
            })

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(entries, f, indent=2)

    print(f"Wrote {len(entries)} entries -> {out_path}")

if __name__ == "__main__":
    main()
