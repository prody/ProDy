try:
    execfile('../../conf.py')
except NameError:
    exec(open('../../conf.py').read())

sys.path[-1] = os.path.abspath('../../_sphinxext')

def is_installed(name):

    try:
        __import__(name)
    except ImportError:
        return False
    else:
        return True

extensions = filter(is_installed, extensions)

prody_inv = '_inventory/prody.inv'
if not os.path.isfile(prody_inv):
    from inventory import trim_labels
    trim_labels('../../_build/html/objects.inv', prody_inv)
intersphinx_mapping['prody'] = ('http://csb.pitt.edu/ProDy', prody_inv)

master_doc = 'index'

version = release = tutorial_version or version

latex_documents = [
    ('index',
     os.path.basename(os.getcwd()) + '.tex',
     tutorial_title,
     tutorial_author,
     'manual'),
]

latex_logo = tutorial_logo or '_static/logo.png'
latex_show_urls = 'footnote'

html_additional_pages = {}
html_domain_indices = False
html_use_index = False

latex_domain_indices = True

latex_appendices = ['acknowledgments']
