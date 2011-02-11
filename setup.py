from distutils.core import setup

readme = open('README.txt')
long_description = ''.join(readme.readlines())
readme.close()

# We now define the ProDy version number in prody/__init__.py
# Here we can't use "import prody" then "prody.__version__" as that would
# tell us the version of Prody already installed (if any).
__version__ = "Undefined"
for line in open('prody/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(
    name='ProDy',
    version=__version__,
    author='Ahmet Bakan',
    author_email='ahb12 at pitt dot edu',
    description='A Python Package for Protein Dynamics Analysis',
    long_description=long_description,
    url='http://www.csb.pitt.edu/ProDy',
    packages=['prody'],
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
    scripts=['scripts/anm.py', 'scripts/gnm.py', 'scripts/pdbselect.py'],
    requires=['NumPy', ],
    provides=['ProDy({0:s})'.format(__version__)]
    )
