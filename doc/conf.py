# -*- coding: utf-8 -*-
#
# ProDy documentation build configuration file, created by
# sphinx-quickstart on Wed Sep 15 21:31:37 2010.

import sys, os, os.path, time

sys.path.append('/home/abakan/Code/ipython/docs')

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
              'sphinxext.ipython_console_highlighting',
              'sphinxext.ipython_directive']
               #, 'sphinxcontrib.spelling']
               #'sphinx.ext.pngmath',

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'contents'

project = u'ProDy'
copyright = u'2010-2013, Ahmet Bakan'

def getRevisionNumber():
    from subprocess import PIPE, Popen
    pipe = Popen('git log --summary'.split(), stdout=PIPE, stderr=PIPE)
    return str(pipe.stdout.read().count('Author:'))
version = '1.4.2'
release =  version


exclude_patterns = ['_build', 'tutorials/*', 'tutorials/template/*']


add_module_names = False

show_authors = True

pygments_style = 'sphinx'

modindex_common_prefix = ['prody.']
doctest_global_setup = "from prody import *"

# -- Options for HTML output ---------------------------------------------------
html_theme = 'sphinxdoc'
html_theme_options = {}
html_theme_path = []

html_title = "ProDy"
html_favicon = '_static/favicon.ico'
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'

html_index = 'index.html'

generic_sidebars = ['docversion.html', 'howtocite.html', 'localtoc.html',
                    'relations.html', 'searchbox.html']
html_sidebars = {
    'index': [],#['slideshow.html', 'docversion.html', 'howtocite.html',
              #'getprody.html', 'credits.html', 'getintouch.html',
              #'searchbox.html',],
    'genindex': ['searchbox.html'],
    'py-modindex': ['searchbox.html'],
    'search': [],
    'tutorial': ['docversion.html', 'howtocite.html', 'localtoc.html',
                 'codesnippets.html', 'credits.html', 'searchbox.html'],
    'bibliography': generic_sidebars,
    'changes': generic_sidebars,
    'contents': generic_sidebars,
    'credits': generic_sidebars,
    'features': generic_sidebars,
    'getprody': ['howtocite.html', 'localtoc.html', 'relations.html',
                 'credits.html', 'searchbox.html'],
    'license': generic_sidebars,
    'publications': generic_sidebars,
    'examples/index': generic_sidebars,
    'reference/index': generic_sidebars,
    'reports/index': generic_sidebars,
    'scripts/index': generic_sidebars,
    'todo': generic_sidebars,
    'plugins/index': generic_sidebars,
    'plugins/getnmwiz': generic_sidebars,
    '**': ['docversion.html', 'howtocite.html', 'localtoc.html',
           'relations.html', 'codesnippets.html', 'searchbox.html']}
#html_sidebars = {'index': ['indexsidebar.html', 'searchbox.html']}

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
ipython_savefig_dir = '_static'

# Plot directive configuration
plot_basedir = os.path.join('_build', 'plot_directive')
plot_working_directory = os.path.join(os.getcwd(), '_doctest')

plot_formats = [('png', 80), ('pdf', 80)]

plot_pre_code = """import numpy as np
from prody import *
from matplotlib import pyplot as plt
"""

plot_template = """
{{ source_code }}

{{ only_html }}


   {% for img in images %}
   .. figure:: {{ build_dir }}/{{ img.basename }}.png
      {%- for option in options %}
      {{ option }}
      {% endfor %}

      {{ caption }}
   {% endfor %}

{{ only_latex }}

   {% for img in images %}
   .. image:: {{ build_dir }}/{{ img.basename }}.pdf
   {% endfor %}

"""
plot_rcparams = {'font.size': 10,
                 'xtick.labelsize': 'small',
                 'ytick.labelsize': 'small',
                 'figure.figsize': [5., 4.],}
plot_apply_rcparams = True


# -- Options for LaTeX output --------------------------------------------------
latex_paper_size = 'letter'
latex_font_size = '10pt'

latex_documents = [
  ('contents', 'ProDy.tex', u'ProDy Documentation',
   u'Ahmet Bakan', 'manual'),
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
autodoc_default_flags = []# ['members', 'undoc-members', 'show-inheritance']
autoclass_content = 'both'
todo_include_todos = True

googleanalytics_enabled = True
googleanalytics_id = 'UA-19801227-1'

intersphinx_mapping = {
    'matplotlib': ('http://matplotlib.sourceforge.net/', None),
    'python': ('http://docs.python.org/', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None)
}

rst_epilog = u"""

.. |pdf| image:: /_static/pdf.png

.. |more| image:: /_static/more.png

.. |example| image:: /_static/example.png

.. |bulb| image:: /_static/bulb.png

.. |new| image:: /_static/new.png

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