# -*- coding: utf-8 -*-
import os
import sys
from time import time

sys.path.append(os.path.abspath('_sphinxext'))
sys.path.append(os.path.abspath('_inventory'))

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
              'matplotlib.sphinxext.plot_directive',
              'matplotlib.sphinxext.only_directives',
              'sphinxcontrib.googleanalytics',
              'sphinxcontrib.googlechart',
              'sphinxcontrib.youtube',
              'ipython_console_highlighting',
              'ipython_directive']
               #, 'sphinxcontrib.spelling']
               #'sphinx.ext.pngmath',

templates_path = ['_theme']
source_suffix = '.rst'
master_doc = 'contents'

project = u'ProDy'
copyright = u'2010-2013, Ahmet Bakan'


def getRevisionNumber():
    from subprocess import PIPE, Popen
    pipe = Popen('git log --summary'.split(), stdout=PIPE, stderr=PIPE)
    return str(pipe.stdout.read().count('Author:'))

import prody
version = release = prody.__version__

exclude_patterns = ['_build', 'examples', 'tutorials/template',
                    '*acknowledgments', 'reports', 'random']


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
    'contents': generic_sidebars,
    'credits': generic_sidebars,
    'features': generic_sidebars,
    'getprody': generic_sidebars,
    'license': generic_sidebars,
    'publications': generic_sidebars,
    'examples/index': generic_sidebars,
    'reference/index': generic_sidebars,
    'reports/index': generic_sidebars,
    'scripts/index': generic_sidebars,
    'todo': generic_sidebars,
    'plugins/index': generic_sidebars,
    'plugins/getnmwiz': generic_sidebars,
    '**': ['toolbox.html', 'releasenotes.html', 'howtocite.html']}
html_additional_pages = {'index': 'index.html'}

html_copy_source = False
html_show_sourcelink = False

html_show_sphinx = True
html_show_copyright = True

extlinks = {
    'issue': ('https://bitbucket.org/abakan/prody/issue/%s', 'issue '),
    'pdb': ('http://www.pdb.org/pdb/explore/explore.do?structureId=%s', ''),
}


# ipython directive configuration
ipython_savefig_dir = os.path.join('_static', 'figures')


# -- Options for LaTeX output --------------------------------------------------
latex_paper_size = 'letter'
latex_font_size = '10pt'

latex_documents = [
    ('contents', 'ProDy.tex', u'ProDy Documentation', u'Ahmet Bakan', 'manual')
]

latex_logo = '_static/logo.png'

latex_show_pagerefs = True
latex_show_urls = 'no'

latex_preamble = ''
latex_appendices = []

latex_domain_indices = True

latex_elements = {
    'classoptions': ',openany,oneside',
    'fontpkg': '\\usepackage{palatino}',
    'babel': '\\usepackage[english]{babel}'
}


autodoc_member_order = 'groupwise'
autodoc_default_flags = []  # ['members', 'undoc-members', 'show-inheritance']
autoclass_content = 'both'
todo_include_todos = True

googleanalytics_enabled = True
googleanalytics_id = 'UA-19801227-1'

intersphinx_mapping = {
    'numpy': ('http://docs.scipy.org/doc/numpy/',
              '_inventory/numpy.inv'),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference/',
              '_inventory/scipy.inv'),
    'python': ('http://docs.python.org/',
               '_inventory/python.inv'),
    'matplotlib': ('http://matplotlib.sourceforge.net/',
                   '_inventory/matplotlib.inv'),
}

week = 7 * 24 * 3600
for pkg, (url, inv) in intersphinx_mapping.items():

    if (not os.path.isfile(inv) or time() - os.path.getmtime(inv) > week):
        import urllib2
        sys.stderr.write('Downloading {} inventory from {}\n'
                         .format(pkg, inv))
        _ = urllib2.urlopen(url + 'objects.inv')
        _ = _.read()
        with open(inv, 'w') as out:
            out.write(_)


rst_epilog = u"""

.. |nmwiz| replace:: http://www.csb.pitt.edu/NMWiz/

.. |vmd| replace:: http://www.ks.uiuc.edu/Research/vmd/

.. |pdb| replace:: http://www.pdb.org/

.. |mdanalysis| replace:: http://code.google.com/p/mdanalysis/

.. |pyparsing| replace:: http://pyparsing.wikispaces.com/

.. |matplotlib| replace:: http://matplotlib.sourceforge.net

.. |biopython| replace:: http://biopython.org/

.. |pypi| replace:: http://pypi.python.org/pypi/ProDy

.. |anm| replace:: http://ignmtest.ccbb.pitt.edu/cgi-bin/anm/anm1.cgi

.. |A2| replace:: Å\ :sup:`2`

.. |tutorials| replace:: http://csb.pitt.edu/ProDy/tutorials/

.. |questions| replace:: To receive new release announcements, join our
   ProDy-News Google Group: http://groups.google.com/group/prody-news

.. |suggestions| replace:: Suggestions, feature requests, or
   problems? Submit them to the Bitbucket issue tracker:
   https://bitbucket.org/abakan/prody/issues
"""
