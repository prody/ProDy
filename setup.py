import os
import sys
import platform
from os import sep as dirsep
from os.path import isfile, join

from setuptools import setup
from setuptools import Extension

import shutil

if sys.version_info[:2] < (2, 7):
    sys.stderr.write('Python 2.6 and older is not supported\n')
    sys.exit()

if sys.version_info[:2] == (2, 7) or sys.version_info[:2] <= (3, 5):
    INSTALL_REQUIRES=['numpy>=1.10', 'biopython<=1.76', 'pyparsing', 'scipy']
else:
    INSTALL_REQUIRES=['numpy>=1.10,<2', 'biopython', 'pyparsing<=3.1.1', 'scipy', 'setuptools']

if sys.version_info[0] == 3 and sys.version_info[1] < 6:
    sys.stderr.write('Python 3.5 and older is not supported\n')
    sys.exit()

if os.name == 'java':
    sys.stderr.write('JavaOS is not supported\n')
    sys.exit()

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: '
                     'http://www.numpy.org/\n')
    sys.exit()

if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 10]:
    sys.stderr.write('numpy v1.10 or later is required, you can find it at: '
                     'http://www.numpy.org/\n')
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
            'prody.chromatin',
            'prody.compounds',
            'prody.domain_decomposition',
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
    'prody.utilities': ['datafiles/*.dat'],
    'prody.tests': ['datafiles/pdb*.pdb',
                    'datafiles/*.dat',
                    'datafiles/*.coo',
                    'datafiles/dcd*.dcd',
                    'datafiles/xml*.xml',
                    'datafiles/msa*',
                    'datafiles/mmcif*cif',],
    'prody.proteins': ['tabulated_energies.txt'],
}

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join(*pkg.split('.'))
    
from glob import glob
tntDir = join('prody', 'utilities', 'tnt')
hpbSoDir = join('prody', 'proteins', 'hpbmodule',
                'hpb_Python{0}.{1}'.format(sys.version_info[0],
                                           sys.version_info[1]))
proteinsDir = join('prody', 'proteins')

try:
    shutil.copy(hpbSoDir + "/hpb.so", proteinsDir)
except FileNotFoundError:
    pass

EXTENSIONS = [
    Extension('prody.dynamics.rtbtools',
              glob(join('prody', 'dynamics', 'rtbtools.c')),
              include_dirs=[numpy.get_include()]),
    Extension('prody.dynamics.smtools',
              glob(join('prody', 'dynamics', 'smtools.c')),
              include_dirs=[numpy.get_include()]),
#    Extension('prody.dynamics.saxstools',
#              glob(join('prody', 'dynamics', 'saxstools.c')),
#              include_dirs=[numpy.get_include()]),
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

# extra arguments for compiling C++ extensions on MacOSX
if platform.system() == 'Darwin':
    os_ver = platform.mac_ver()[0]
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = os_ver
    os.environ['CC'] = 'clang'
    os.environ['CXX'] = 'clang++'
    #extra_compile_args.append('-stdlib=libc++')


# extra compilation of reg_tet.f (hpb):
# import subprocess
# subprocess.call(['gfortran', '-O3', '-fPIC', '-c',
#                  join('prody', 'proteins', 'hpbmodule', 'reg_tet.f'),
#                  '-o', join('prody', 'proteins', 'hpbmodule', 'reg_tet.o')])

CONTRIBUTED = [
    Extension('prody.kdtree._CKDTree',
              [join('prody', 'kdtree', 'KDTree.c'),
               join('prody', 'kdtree', 'KDTreemodule.c')],
              include_dirs=[numpy.get_include()]),
    Extension('prody.proteins.ccealign', 
              [join('prody', 'proteins', 'ccealign', 'ccealignmodule.cpp')], 
              include_dirs=[tntDir], language='c++'),
    #Extension('prody.proteins.hpb',
    #          [join('prody', 'proteins', 'hpbmodule', 'reg_tet.c')],
    #          include_dirs=[hpbDir], language='c++',
    #          extra_compile_args=['-O3', '-fPIC'],
    #          extra_objects=[join(hpbDir, 'libf2c', 'libf2c.a')]
    #          )
]

for ext in CONTRIBUTED:
    if all([isfile(src) for src in ext.sources]):
        EXTENSIONS.append(ext)

# SCRIPTS = ['scripts/prody', 'scripts/evol']
# if (platform.system() == 'Windows' or
#     len(sys.argv) > 1 and sys.argv[1] not in ('build', 'install')):
#     for script in list(SCRIPTS):
#         SCRIPTS.append(script + '.bat')


SCRIPTS = ['prody=prody.apps:prody_main', 'evol=prody.apps:evol_main']

setup(
    name='ProDy',
    version=__version__,
    author='James Krieger, Karolina Mikulska-Ruminska, She Zhang, Hongchun Li, Cihan Kaya, Ahmet Bakan, and others',
    author_email='jamesmkrieger@gmail.com',
    description='A Python Package for Protein Dynamics Analysis',
    long_description=long_description,
    url='http://www.csb.pitt.edu/ProDy',
    packages=PACKAGES,
    #package_dir=PACKAGE_DIR,
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
    #scripts=SCRIPTS,
    entry_points = {
        'console_scripts': SCRIPTS,
    },
    install_requires=INSTALL_REQUIRES,
    #provides=['ProDy ({0:s})'.format(__version__)]
)
