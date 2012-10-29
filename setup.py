import os
import sys
from os import sep as dirsep
from os.path import isfile, join
from shutil import copy

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install

PY3K = sys.version_info[0] > 2

with open('README.rst') as inp:
    long_description = inp.read()

__version__ = ''
inp = open('lib/prody/__init__.py')
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

PACKAGES = ['prody', 'prody.atomic', 'prody.dynamics', 'prody.ensemble',
            'prody.measure', 'prody.proteins', 'prody.trajectory',
            'prody.routines', 'prody.utilities']
PACKAGE_DATA = {}
if sys.version_info[:2] > (2, 6):
    PACKAGES.extend(['prody.tests',
                     'prody.tests.test_atomic',
                     'prody.tests.test_datafiles',
                     'prody.tests.test_dynamics',
                     'prody.tests.test_ensemble', 
                     'prody.tests.test_kdtree', 
                     'prody.tests.test_measure',
                     'prody.tests.test_proteins',
                     'prody.tests.test_trajectory',
                     'prody.tests.test_utilities',])
    PACKAGE_DATA['prody.tests'] = ['test_datafiles/pdb*.pdb', 
                                   'test_datafiles/*.dat', 
                                   'test_datafiles/*.coo', 
                                   'test_datafiles/dcd*.dcd',
                                   'test_datafiles/xml*.xml']
PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join('lib', *pkg.split('.'))

EXTENSIONS = []

if os.name != 'java' and sys.version_info[0] == 2:
    pairwise2 = ['cpairwise2.c', 'pairwise2.py']
    if all([isfile(join('lib', 'prody', 'proteins', fn)) 
                           for fn in pairwise2]):  
        EXTENSIONS.append(
            Extension('prody.proteins.cpairwise2', 
                      [join('lib', 'prody', 'proteins', 'cpairwise2.c')],
                      ))
    else:
        raise Exception('one or more pairwise2 module files are missing')
    if isInstalled('numpy'):
        import numpy
        kdtree_files = ['__init__.py', 'KDTree.c', 'KDTree.h',
                        'KDTreemodule.c', 'Neighbor.h', 'kdtree.py']
        if all([isfile(join('lib', 'prody', 'kdtree', fn)) 
                for fn in kdtree_files]):
            EXTENSIONS.append(
                Extension('prody.kdtree._CKDTree',
                          [join('lib', 'prody', 'kdtree', 'KDTree.c'),
                           join('lib', 'prody', 'kdtree', 'KDTreemodule.c')],
                          include_dirs=[numpy.get_include()],
                          ))
        else:
            raise Exception('one or more kdtree module files are missing')
        PACKAGES.append('prody.kdtree')
    elif isInstalled('numpy'):
        raise ImportError('numpy is not installed')

SCRIPTS = ['scripts/prody']
import platform
if (platform.system() == 'Windows' or 
    len(sys.argv) > 1 and sys.argv[1] not in ('build', 'install')):
    SCRIPTS.append('scripts/prody.bat')


if len(sys.argv) > 1 and sys.argv[1] == 'copy':
    from glob import glob
    ext = '.pyd' if platform == 'Windows' else '.so'
    base = join('build', 
                'lib*' + '.'.join([str(i) for i in sys.version_info[:2]]))
    for src in glob(join(base, 'prody', '*', '*' + ext)):
        dst = join(*src.split(dirsep)[2:])
        dst = join('lib', dst)
        sys.stderr.write('cp ' + src + ' ' + dst + '\n')
        copy(src, dst)
        
else:
    setup(
        name='ProDy',
        version=__version__,
        author='Ahmet Bakan',
        author_email='ahb12 at pitt dot edu',
        description='A Python Package for Protein Dynamics Analysis',
        long_description=long_description,
        url='http://www.csb.pitt.edu/ProDy',
        packages=PACKAGES,
        package_dir=PACKAGE_DIR,
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
