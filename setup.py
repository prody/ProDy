import glob
import os
from os.path import isfile, join
import sys

from distutils.core import setup
from distutils.extension import Extension

readme = open('README.txt')
long_description = ''.join(readme.readlines())
readme.close()

__version__ = ''
for line in open('prody/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())
        break
    
def isInstalled(module_name):
    """Check if a required package is installed, by trying to import it."""
    try:
        return __import__(module_name)
    except ImportError:
        return False
    else:
        return True

if not isInstalled('numpy'):    
    print("""NumPy is not installed. This package is required for main ProDy 
features and needs to be installed before you can use ProDy.  
You can find NumPy at: http://numpy.scipy.org""")
    
PACKAGES = ['prody']

EXTENSIONS = []

if os.name != 'java' and sys.version_info[0] == 2:
    pairwise2 = ['cpairwise2module.c', 'pairwise2.py']
    if all([isfile(join('prody', fn)) for fn in pairwise2]):  
        EXTENSIONS.append(
            Extension('prody.cpairwise2',
                      ['prody/cpairwise2module.c'],
                      include_dirs=["prody"]
                      ))
    if isInstalled('numpy'):
        import numpy
        kdtree = ['__init__.py', 'KDTree.c', 'KDTree.h', 'KDTree.py', 
                  'KDTreemodule.c', 'Neighbor.h']
        if all([isfile(join('prody/KDTree', fn)) for fn in kdtree]):
            EXTENSIONS.append(
                Extension('prody.KDTree._CKDTree',
                          ['prody/KDTree/KDTree.c',
                           'prody/KDTree/KDTreemodule.c'],
                          include_dirs=[numpy.get_include()],
                          ))
        PACKAGES.append('prody.KDTree')

SCRIPTS = glob.glob('scripts/*py')

setup(
    name='ProDy',
    version=__version__,
    author='Ahmet Bakan',
    author_email='ahb12 at pitt dot edu',
    description='A Python Package for Protein Dynamics Analysis',
    long_description=long_description,
    url='http://www.csb.pitt.edu/ProDy',
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    license='GPLv3',
    keywords=('protein, dynamics, elastic network model, gaussian network model, '
              'anisotropic network model, essential dynamics analysis, '
              'principal component analysis, ProteinDataBank, PDB, GNM, ANM, PCA'),
    classifiers=[
                 'Development Status :: 4 - Beta',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License (GPL)',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Chemistry',
                ],
    scripts=SCRIPTS,
    requires=['NumPy', ],
    provides=['ProDy({0:s})'.format(__version__)]
    )
