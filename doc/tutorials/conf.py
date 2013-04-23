try:
    execfile('../../conf.py')
except NameError:
    exec(open('../../conf.py').read())
    
version = release = tutorial_prody_version or version   
intersphinx_mapping['prody'] = ('http://csb.pitt.edu/ProDy/', None)

master_doc = 'index'

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

plot_working_directory = os.path.join(os.getcwd(), '..', '..', '_doctest')
