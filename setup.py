import os
import sys
import platform
from os import sep as dirsep
from os.path import isfile, join

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install

if sys.version_info[:2] < (2, 6):
    sys.stderr.write('Python 2.5 and older is not supported\n')
    sys.exit()

if os.name == 'java':
    sys.stderr.write('JavaOS is not supported\n')
    sys.exit()

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 4]:
    sys.stderr.write('numpy v1.4 or later is required, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()


__version__ = ''
with open('prody/__init__.py') as inp:
  for line in inp:
      if line.startswith('__version__'):
          exec(line.strip())
          break

with open('README.rst') as inp:
    long_description = inp.read()


PACKAGES = ['prody',
            'prody.atomic',
            'prody.database',
            'prody.dynamics',
            'prody.ensemble',
            'prody.kdtree',
            'prody.measure',
            'prody.proteins',
            'prody.sequence',
            'prody.trajectory',
            'prody.utilities',
            'prody.apps',
            'prody.apps.prody_apps',
            'prody.apps.evol_apps',
            'prody.tests',
            'prody.tests.apps',
            'prody.tests.atomic',
            'prody.tests.datafiles',
            'prody.tests.dynamics',
            'prody.tests.ensemble',
            'prody.tests.kdtree',
            'prody.tests.measure',
            'prody.tests.proteins',
            'prody.tests.sequence',
            'prody.tests.trajectory',
            'prody.tests.utilities',]
PACKAGE_DATA = {
    'prody.tests': ['datafiles/pdb*.pdb',
                    'datafiles/*.dat',
                    'datafiles/*.coo',
                    'datafiles/dcd*.dcd',
                    'datafiles/xml*.xml',
                    'datafiles/msa*',]
}

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join(*pkg.split('.'))
from glob import glob
EXTENSIONS = [
    Extension('prody.dynamics.rtbtools',
              glob(join('prody', 'dynamics', 'rtbtools.c')),
              include_dirs=[numpy.get_include()]),
    Extension('prody.dynamics.smtools',
              glob(join('prody', 'dynamics', 'smtools.c')),
              include_dirs=[numpy.get_include()]),
    Extension('prody.dynamics.saxstools',
              glob(join('prody', 'dynamics', 'saxstools.c')),
              include_dirs=[numpy.get_include()]),
    Extension('prody.sequence.msatools',
              [join('prody', 'sequence', 'msatools.c'),],
              include_dirs=[numpy.get_include()]),
    Extension('prody.sequence.msaio',
              [join('prody', 'sequence', 'msaio.c'),],
              include_dirs=[numpy.get_include()]),
    Extension('prody.sequence.seqtools',
              [join('prody', 'sequence', 'seqtools.c'),],
              include_dirs=[numpy.get_include()]),
]

CONTRIBUTED = [
    Extension('prody.proteins.cpairwise2',
              [join('prody', 'proteins', 'cpairwise2.c')]),
    Extension('prody.kdtree._CKDTree',
              [join('prody', 'kdtree', 'KDTree.c'),
               join('prody', 'kdtree', 'KDTreemodule.c')],
              include_dirs=[numpy.get_include()]),
]

for ext in CONTRIBUTED:
    if all([isfile(src) for src in ext.sources]):
        EXTENSIONS.append(ext)

SCRIPTS = ['scripts/prody', 'scripts/evol']
if (platform.system() == 'Windows' or
    len(sys.argv) > 1 and sys.argv[1] not in ('build', 'install')):
    for script in list(SCRIPTS):
        SCRIPTS.append(script + '.bat')

setup(
    name='ProDy',
    version=__version__,
    author='Cihan Kaya',
    author_email='cihank@pitt.edu',
    description='A Python Package for Protein Dynamics Analysis',
    long_description=long_description,
    url='http://www.csb.pitt.edu/ProDy',
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    package_data=PACKAGE_DATA,
    ext_modules=EXTENSIONS,
    license='MIT License',
    keywords=('protein, dynamics, elastic network model, '
              'Gaussian network model, anisotropic network model, '
              'essential dynamics analysis, principal component analysis, '
              'Protein Data Bank, PDB, GNM, ANM, SM, PCA'),
    classifiers=[
                 'Development Status :: 5 - Production/Stable',
                 'Intended Audience :: Education',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: MIT License',
                 'Operating System :: MacOS',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 2',
                 'Programming Language :: Python :: 3',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Chemistry',
                ],
    scripts=SCRIPTS,
    requires=['NumPy (>=1.7)', ],
    provides=['ProDy ({0:s})'.format(__version__)]
)
