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
