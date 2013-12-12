# -*- coding: utf-8 -*-
import os
import sys


sys.path.append(os.path.abspath('sphinxext'))
RTD = os.environ.get('READTHEDOCS', None) == 'True'
if RTD:
    tags.add('rtd')

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

source_suffix = '.rst'
master_doc = 'index'

project = u'ProDy'
copyright = u'2010-2014, University of Pittsburgh'


__version__ = ''
for _ in ['../prody/__init__.py', # building docs
          'ProDy/prody/__init__.py', # building complete website
          '../../ProDy/prody/__init__.py']: # building tutorial PDFs
    if os.path.isfile(_):
        with open(_) as inp:
            statement = ''.join([line for line in inp
                                 if line[:3] in ('__v', '__r')])
            exec(statement)
version, release = __version__, __release__


if release.endswith('dev'):
    import shlex
    from subprocess import Popen, PIPE
    _ = '../.git'
    if not os.path.isdir(_):
        _ = 'ProDy/.git'
    args = shlex.split('git --git-dir {} describe --tags --abbrev=0'.format(_))
    tag = Popen(args, stdout=PIPE, stderr=PIPE)
    rst_prolog = """

.. only:: html

    .. note::

        This documentation is for a development version of ProDy.
        There may be significant differences from the latest stable
        release ({}).

    """.format(tag.communicate()[0].strip())

exclude_patterns = ['_build', '_workdir']


add_module_names = False

show_authors = True

pygments_style = 'sphinx'

modindex_common_prefix = ['prody.']
doctest_global_setup = "from prody import *"

# -- Options for HTML output ---------------------------------------------------
if RTD:
    html_theme = 'default'
else:
    templates_path = ['_theme']
    html_theme = '_theme'
    html_theme_path = ['.']
    html_static_path = ['_static']
html_theme_options = {}

html_title = "ProDy"
html_favicon = '_static/favicon.ico'
html_last_updated_fmt = '%b %d, %Y'

html_index = 'index.html'

generic_sidebars = ['toolbox.html', 'releasenotes.html', 'howtocite.html']
html_sidebars = {'**': generic_sidebars}

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
.. _nose: http://nose.readthedocs.org


.. |A2| replace:: Å\ :sup:`2`
"""
