try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())

sys.path[-1] = os.path.abspath('../_sphinxext')

prody_inv = '_inventory/prody_noapi.inv'
if not os.path.isfile(prody_inv):
    from inventory import remove_api
    remove_api('../_build/html/objects.inv', prody_inv)
intersphinx_mapping['prody'] = ('http://csb.pitt.edu/ProDy', prody_inv)

master_doc = 'index'

latex_documents = [
  ('index',
   'ProDy.tex',
   'ProDy Reference Manual',
   'Ahmet Bakan',
   'manual'),
]
latex_show_urls = 'footnote'

html_additional_pages = {}

latex_domain_indices = True

plot_working_directory = os.path.join(os.getcwd(), '..', '_doctest')