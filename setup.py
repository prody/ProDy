readme = open('README.txt')
long_description = ''.join(readme.readlines())
readme.close()

from distutils.core import setup
setup(name='ProDy',
      version='0.2.0a1',
      author='Ahmet Bakan',
      author_email='ahb12@pitt.edu',
      description='A Python Package for Protein Structural Dynamics Analysis',
      long_description=long_description,
      url='http://www.pitt.edu/~ahb12/software/ProDy',
      packages=['prody', 'prody.dynamics', 'prody.proteins'],
      license='GPLv3',
      keywords=('protein, dynamics, elastic network model, gaussian network model, '
                'anisotropic network model, essential dynamics analysis, '
                'principal component analysis, ProteinDataBank, PDB, GNM, ANM, PCA'),
      classifiers=[
               'Development Status :: 3 - Alpha',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License (GPL)',
               'Operating System :: MacOS',
               'Operating System :: Microsoft :: Windows',
               'Operating System :: POSIX',
               'Programming Language :: Python',
               'Programming Language :: Python :: 2.6',
               'Topic :: Scientific/Engineering :: Bio-Informatics',
               'Topic :: Scientific/Engineering :: Chemistry',
               ],
      requires=['NumPy', ],
      provides=['ProDy(0.2.0a1)']
      )
