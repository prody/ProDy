tutorial_title = u'TITLE'
tutorial_author = u'AUTHOR'
tutorial_logo = u''           # default is ProDy logo
tutorial_version = u''        # default is latest ProDy version

# keep the following part as is
try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())