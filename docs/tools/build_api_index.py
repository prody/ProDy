import json
import os
import pkgutil
import importlib
import inspect

OUTFILE = os.path.join(os.path.dirname(__file__), "..", "_static", "api_index.json")

def is_tests_module(modname: str) -> bool:
    return modname.startswith("prody.tests")

def module_to_html(modname: str) -> str:
    # prody.proteins.pdbfile -> reference/proteins/pdbfile.html
    parts = modname.split(".")
    if parts[0] == "prody":
        parts = parts[1:]
    return "reference/" + "/".join(parts) + ".html"

def main():
    import prody  # must be importable during doc build

    # We store a list so partial search can show multiple matches
    entries = []  # {name, qual, url}

    for m in pkgutil.walk_packages(prody.__path__, prefix="prody."):
        modname = m.name
        if is_tests_module(modname):
            continue

        try:
            mod = importlib.import_module(modname)
        except Exception:
            continue

        html = module_to_html(modname)

        for name, obj in inspect.getmembers(mod):
            if name.startswith("_"):
                continue

            try:
                if getattr(obj, "__module__", None) != modname:
                    continue
            except Exception:
                continue

            if inspect.isfunction(obj) or inspect.isclass(obj):
                qual = f"{modname}.{name}"
                url = f"{html}#{qual}"
                entries.append({"name": name, "qual": qual, "url": url})

    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    with open(OUTFILE, "w", encoding="utf-8") as f:
        json.dump(entries, f, indent=2)

    print(f"Wrote {len(entries)} entries to {OUTFILE}")

if __name__ == "__main__":
    main()
