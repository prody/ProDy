try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())

sys.path[-1] = os.path.abspath('../_sphinxext')

prody_inv = '_pkginv/prody.inv'

if not os.path.isfile(prody_inv):

    def trim():
        inp = open('../../_build/html/objects.inv', 'rb')
        header = [inp.readline().decode('utf-8') for i in range(4)]
        if 'zlib' not in header[-1]:
            raise ValueError
        import zlib
        trimmed = []

        for line in zlib.decompress(inp.read()).decode('utf-8').splitlines():
            items = line.split(' py:')

            if len(items) > 1:
                if (items[1].startswith('method') or
                    items[1].startswith('attribute')):
                    items[0] = '.'.join(items[0].split('.')[-2:])
                elif items[1].startswith('module'):
                    _ = items[0].split('.')
                    if len(_) > 1:
                        if _[-1] == _[-2]:
                            items[0] = '.'.join(_[-2:] or _)
                        else:
                            items[0] = _[-1]
                else:
                    items[0] = items[0].split('.')[-1]
                trimmed.append(items[0] + ' py:' + items[1])
            else:
                trimmed.append(line)
        inp.close()
        out = open(prody_inv, 'wb')
        for line in header:
            out.write(line.encode('utf-8'))
        compressor = zlib.compressobj(9)

        for line in trimmed:

            out.write(compressor.compress((line + '\n').encode('utf-8')))
        out.write(compressor.flush())
        out.close()
        return

    trim()

    def get_inv():
        from sphinx.ext.intersphinx import read_inventory_v2
        f = open(prody_inv)
        f.readline()
        return read_inventory_v2(f, '.', os.path.join)

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

plot_working_directory = os.path.join(os.getcwd(), '..', '..', '_doctest')

latex_appendices = ['acknowledgments']