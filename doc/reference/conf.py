try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())

sys.path[-2] = os.path.abspath('../_sphinxext')
sys.path[-1] = os.path.abspath('../_inventory')

prody_inv = '_inventory/prody_noapi.inv'
if not os.path.isfile(prody_inv):
    from inventory import remove_api
    remove_api('../_build/html/objects.inv', prody_inv)
intersphinx_mapping['prody'] = ('http://csb.pitt.edu/ProDy', prody_inv)

master_doc = 'index'

latex_documents = [
    ('index',
     'ProDy_Manual.tex',
     'ProDy Reference Manual',
     'Ahmet Bakan',
     'manual'),
]

latex_show_urls = 'footnote'
latex_domain_indices = True

html_additional_pages = {}