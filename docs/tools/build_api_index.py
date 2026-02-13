# -*- coding: utf-8 -*-
import json
import os
import sys
import importlib
from pathlib import Path


def iter_module_names(pkg_name: str, pkg_root: Path):
    """
    Walk the package directory on disk WITHOUT importing subpackages.
    This avoids importing prody.tests.* which breaks on RTD.
    """
    for py in pkg_root.rglob("*.py"):
        # Skip caches / hidden
        if "__pycache__" in py.parts:
            continue

        rel = py.relative_to(pkg_root)

        # Skip tests entirely (this fixes your NRAS_BRAF_GDP_model_0.cif crash)
        if "tests" in rel.parts:
            continue

        # Build module name
        parts = list(rel.parts)
        if parts[-1] == "__init__.py":
            parts = parts[:-1]
        else:
            parts[-1] = parts[-1].replace(".py", "")

        if not parts:
            yield pkg_name
        else:
            yield pkg_name + "." + ".".join(parts)


def main():
    import prody

    pkg_name = "prody"
    pkg_root = Path(prody.__file__).resolve().parent

    here = Path(__file__).resolve().parent.parent  # docs/
    out_dir = here / "_static"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "api_index.json"

    results = []

    for modname in sorted(set(iter_module_names(pkg_name, pkg_root))):
        try:
            mod = importlib.import_module(modname)
        except Exception:
            # Some modules may require optional deps or side effects.
            # We skip them so docs build doesn't fail.
            continue

        for name, obj in vars(mod).items():
            if callable(obj) and getattr(obj, "__module__", "") == modname:
                qual = f"{modname}.{name}"
                results.append(
                    {
                        "name": name,          # parsePDB
                        "qualname": qual,      # prody.proteins.pdbfile.parsePDB
                    }
                )

    with out_file.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"[smart_search] wrote {len(results)} entries to {out_file}")


if __name__ == "__main__":
    main()
