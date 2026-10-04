# -*- coding: utf-8 -*-
import os
import sys
import sphinx_rtd_theme

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here.
# Point to the root of the repo to find the 'prody' package
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------
project = 'ProDy'
copyright = '2010-2026, Bahar Lab'
author = 'Bahar Lab'

# The short X.Y version (optional: you can import prody to get this dynamically)
# import prody
# version = prody.__version__
# release = prody.__version__
version = '2.5'
release = '2.5.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings.
extensions = [
    # Core Sphinx extensions
    'sphinx.ext.autodoc',       # Generates docs from docstrings
    'sphinx.ext.autosummary',   # Generates summary tables
    'sphinx.ext.doctest',       # Tests code snippets in docs
    'sphinx.ext.mathjax',       # Renders LaTeX math
    'sphinx.ext.viewcode',      # Adds links to source code
    'sphinx.ext.napoleon',      # Parses NumPy-style docstrings (CRITICAL for ProDy)
    'sphinx.ext.intersphinx',   # Links to other libraries' docs

    # Third-party extensions
    'sphinx_rtd_theme',         # The theme itself
    'sphinxcontrib.jquery',     # FIXED: Restores jQuery to make Search work
    'nbsphinx',                 # Renders Jupyter Notebooks
    'sphinx_copybutton',        # Adds "copy" button to code blocks
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_rtd_theme'

# 1. PATH TO YOUR LOGO (Relative to the 'docs' folder)
html_logo = "_static/logo.png"

# 2. THEME OPTIONS
html_theme_options = {
    # 'logo_only': True,        # Set to True if you want to hide the text "ProDy" entirely
    'logo_only': False,       # Set to False if you want Logo + Text
    'display_version': True,  # Show the version number (e.g., v2.5) below the logo
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
}

# Ensure this is set so Sphinx looks in the _static folder
html_static_path = ['_static']
html_js_files = [
    "smart_search.js",
]

# -- Extension Configuration -------------------------------------------------

# 1. Napoleon Settings (for NumPy docstrings)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# 2. nbsphinx Settings (Jupyter Notebooks)
# Prevent notebooks from running during the build (prevents timeouts)
nbsphinx_execute = 'never'

# 3. Intersphinx mapping (Links to external docs)
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
}
import os
import sys
import json
import os
from sphinx.util.inventory import InventoryFile

def build_api_index(app, exception):
    # If build failed, don't do anything
    if exception is not None:
        return

    inv_path = os.path.join(app.outdir, "objects.inv")
    if not os.path.isfile(inv_path):
        print("WARNING: objects.inv not found, skipping api_index.json generation")
        return

    with open(inv_path, "rb") as f:
        inv = InventoryFile.load(f, "", lambda base, uri: uri)

    # Map lowercase simple name -> URL
    # Also keep full names in case you want them later
    exact_simple = {}
    exact_full = {}

    # Prefer functions over methods if there's a collision
    role_priority = ["py:function", "py:method"]

    for role in role_priority:
        if role not in inv:
            continue
        for fullname, data in inv[role].items():
            # data = (project, version, uri, dispname)
            uri = data[2]
            # Sphinx inventory sometimes uses "$" as placeholder for fullname
            uri = uri.replace("$", fullname)

            simple = fullname.split(".")[-1].lower()
            full_lower = fullname.lower()

            # fill full map always
            exact_full[full_lower] = uri

            # fill simple map only if empty (function wins because we iterate functions first)
            if simple not in exact_simple:
                exact_simple[simple] = uri

    out_path = os.path.join(app.outdir, "_static", "api_index.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump({"exact_simple": exact_simple, "exact_full": exact_full}, f)

    print(f"Wrote {out_path} with {len(exact_simple)} simple entries")

def setup(app):
    app.connect("build-finished", build_api_index)
