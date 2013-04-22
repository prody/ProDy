try:
    execfile('../../conf.py')
except NameError:
    exec(open('../../conf.py').read())
    
intersphinx_mapping['prody'] = ('http://csb.pitt.edu/ProDy/', None)

master_doc = 'index'

latex_documents = [
  ('index', 
   os.path.basename(os.getcwd()) + '.tex', 
   tutorial_title,
   tutorial_author, 
   'manual'),
]

html_additional_pages = {}

plot_working_directory = os.path.join(os.getcwd(), '..', '..', '_doctest')
