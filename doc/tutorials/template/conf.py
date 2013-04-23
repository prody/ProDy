tutorial_title = u'A ProDy Tutorial'
tutorial_author = u'Author One, Author Two'
tutorial_logo = u''                 # default is ProDy logo
tutorial_prody_version = u''        # default is latest ProDy version
# keep the following part as is

try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())
