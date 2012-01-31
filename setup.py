import glob
import os
import os.path
import sys
import shutil
import cPickle
from types import StringType, UnicodeType

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install

PY3K = sys.version_info[0] > 2

with open('README.rst') as inp:
    long_description = inp.read()

__version__ = ''
inp = open('prody/__init__.py')
for line in inp:
    if (line.startswith('__version__')):
        exec(line.strip())
        break
inp.close()

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
PACKAGE_DATA = {}
if sys.version_info[:2] > (2,6):
    PACKAGES.append('prody.tests')
    PACKAGE_DATA['prody.tests'] = ['data/pdb*.pdb', 'data/*.dat', 
                                   'data/*.coo', 'data/dcd*.dcd']

EXTENSIONS = []

if os.name != 'java' and sys.version_info[0] == 2:
    pairwise2 = ['cpairwise2.c', 'pairwise2.py']
    if all([os.path.isfile(os.path.join('prody', fn)) for fn in pairwise2]):  
        EXTENSIONS.append(
            Extension('prody.cpairwise2', ['prody/cpairwise2.c'],
                      include_dirs=["prody"]
                      ))
    if isInstalled('numpy'):
        import numpy
        kdtree = ['__init__.py', 'KDTree.c', 'KDTree.h', 'KDTree.py', 
                  'KDTreemodule.c', 'Neighbor.h']
        if all([os.path.isfile(os.path.join('prody/KDTree', fn)) 
                for fn in kdtree]):
            EXTENSIONS.append(
                Extension('prody.KDTree._CKDTree',
                          ['prody/KDTree/KDTree.c',
                           'prody/KDTree/KDTreemodule.c'],
                          include_dirs=[numpy.get_include()],
                          ))
        PACKAGES.append('prody.KDTree')

SCRIPTS = glob.glob('scripts/*py')
SCRIPTS.append('scripts/prody')

setup(
    name='ProDy',
    version=__version__,
    author='Ahmet Bakan',
    author_email='ahb12 at pitt dot edu',
    description='A Python Package for Protein Dynamics Analysis',
    long_description=long_description,
    url='http://www.csb.pitt.edu/ProDy',
    packages=PACKAGES,
    package_data=PACKAGE_DATA,
    ext_modules=EXTENSIONS,
    license='GPLv3',
    keywords=('protein, dynamics, elastic network model, '
              'Gaussian network model, anisotropic network model, '
              'essential dynamics analysis, principal component analysis, '
              'Protein Data Bank, PDB, GNM, ANM, PCA'),
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
