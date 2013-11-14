# -*- coding: utf-8 -*-
import os
import sys


sys.path.append(os.path.abspath('sphinxext'))
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

extensions = ['sphinx.ext.todo',
              'sphinx.ext.autodoc',
              'sphinx.ext.doctest',
              'sphinx.ext.coverage',
              'sphinx.ext.extlinks',
              'sphinx.ext.graphviz',
              'sphinx.ext.ifconfig',
              'sphinx.ext.viewcode',
              'sphinx.ext.intersphinx',
              'sphinx.ext.inheritance_diagram',
              'matplotlib.sphinxext.mathmpl',
              'matplotlib.sphinxext.only_directives',
              'IPython.sphinxext.ipython_console_highlighting',
              'IPython.sphinxext.ipython_directive',
              'googleanalytics',]

templates_path = ['_theme']
source_suffix = '.rst'
master_doc = 'index'

project = u'ProDy'
copyright = u'2010-2013, Ahmet Bakan'


#import prody
version = release = '1.5-dev' #prody.__version__

exclude_patterns = ['_build']


add_module_names = False

show_authors = True

pygments_style = 'sphinx'

modindex_common_prefix = ['prody.']
doctest_global_setup = "from prody import *"

# -- Options for HTML output ---------------------------------------------------
html_theme = '_theme'
html_theme_options = {}
html_theme_path = ['.']

html_title = "ProDy"
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'

html_index = 'index.html'

generic_sidebars = ['toolbox.html', 'releasenotes.html', 'howtocite.html']
html_sidebars = {
    'index': generic_sidebars,
    'genindex': generic_sidebars,
    'py-modindex': generic_sidebars,
    'search': generic_sidebars,
    'tutorial': generic_sidebars,
    'bibliography': generic_sidebars,
    'changes': generic_sidebars,
    'credits': generic_sidebars,
    'features': generic_sidebars,
    'getprody': generic_sidebars,
    'license': generic_sidebars,
    'examples/index': generic_sidebars,
    'reference/index': generic_sidebars,
    'reports/index': generic_sidebars,
    'scripts/index': generic_sidebars,
    'todo': generic_sidebars,
    'plugins/index': generic_sidebars,
    'plugins/getnmwiz': generic_sidebars,
    '**': ['toolbox.html', 'releasenotes.html', 'howtocite.html']}

html_copy_source = False
html_show_sourcelink = False

html_show_sphinx = True
html_show_copyright = True

extlinks = {
    'issue': ('https://github.com/abakan/ProDy/issues/%s', 'issue '),
    'bbissue': ('https://bitbucket.org/abakan/prody/issue/%s', 'issue '),
    'pdb': ('http://www.pdb.org/pdb/explore/explore.do?structureId=%s', ''),
    'wiki': ('http://en.wikipedia.org/wiki/%s', ''),
    'pfam': ('http://pfam.sanger.ac.uk/family/%s', ''),
    'pfamprotein': ('http://pfam.sanger.ac.uk/protein/%s', ''),
    'uniprot': ('http://www.uniprot.org/uniprot/%s', ''),
    'pdbhet': ('http://www.pdb.org/pdb/ligand/ligandsummary.do?hetId=%s', ''),
}


# ipython directive configuration
ipython_savefig_dir = os.path.join('_static', 'figures')


# -- Options for LaTeX output --------------------------------------------------
latex_paper_size = 'letter'
latex_font_size = '10pt'

latex_documents = [
    ('index', 'ProDy.tex',
     u'ProDy Documentation',
     u'Ahmet Bakan', 'manual')
]

latex_logo = '_static/logo.png'

latex_show_pagerefs = True
latex_show_urls = 'footnote'

latex_preamble = ''
latex_appendices = []

latex_domain_indices = True

latex_elements = {
    'classoptions': ',openany,oneside',
    'fontpkg': '\\usepackage{palatino}',
    'babel': '\\usepackage[english]{babel}'
}


autodoc_member_order = 'groupwise'
autodoc_default_flags = []
autoclass_content = 'both'
todo_include_todos = True

googleanalytics_enabled = True
googleanalytics_id = 'UA-19801227-1'

intersphinx_mapping = {
    'python': ('http://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
    'matplotlib': ('http://matplotlib.sourceforge.net/', None),
    'prodywebsite': ('http://prody.csb.pitt.edu/', None),
}

rst_epilog = u"""

.. _ProDy: http://prody.csb.pitt.edu
.. _Tutorials: http://prody.csb.pitt.edu/tutorials
.. _NMWiz: http://csb.pitt.edu/NMWiz
.. _VMD: http://www.ks.uiuc.edu/Research/vmd
.. _PDB: http://www.pdb.org

.. _MDAnalysis: http://code.google.com/p/mdanalysis
.. _pyparsing: http://pyparsing.wikispaces.com
.. _Matplotlib: http://matplotlib.org
.. _Biopython: http://biopython.org

.. _PyPI: http://pypi.python.org/pypi/ProDy
.. _GitHub: http://github.com/prody/ProDy

.. _IPython: http://ipython.org
.. _Python: http://www.python.org
.. _NumPy: http://www.numpy.org
.. _Scipy: http://www.scipy.org
.. _pip: http://www.pip-installer.org


.. |A2| replace:: Å\ :sup:`2`

.. |questions| replace:: To receive new release announcements, join our
   ProDy-News Google Group: http://groups.google.com/group/prody-news

.. |suggestions| replace:: Suggestions, feature requests, or
   problems? Submit them to the Bitbucket issue tracker:
   https://bitbucket.org/abakan/prody/issues
"""
