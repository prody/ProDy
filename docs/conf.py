import sphinx_rtd_theme

install_requires=[
    'numpy',
    'scipy',
    'matplotlib',
    'pyparsing',
    'biopython',
    'requests', # if you use it
],

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

extensions = [
    'sphinx.ext.autodoc',      # Reads your Python code
    'sphinx.ext.autosummary',  # Generates summary tables
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',      # Renders math equations
    'sphinx.ext.viewcode',     # Adds links to source code
    'sphinx.ext.napoleon',     # <--- CRITICAL: Parses NumPy/Scientific docstrings
    # 'sphinx.ext.intersphinx', # Optional: links to python docs
    'sphinx_rtd_theme',    # Ensure this is in the list too!
    'sphinxcontrib.jquery', # <--- ADD THIS

]

# Napoleon settings (to handle ProDy's scientific docstrings)
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
