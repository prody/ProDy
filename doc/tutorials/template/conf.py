tutorial_title = u'Tutorial Title'
tutorial_author = u'Author'

tutorial_style = u'howto' # 'manual' or 'howto' (more compact PDF output)

try:
    execfile('../conf.py')
except NameError:
    exec(open('../conf.py').read())
