from distutils.core import setup
setup(name='ProDy',
      version='0.5',
      description='Protein Dynamics Analysis',
      long_description='A Python Package for Protein Structural Dynamics Analysis',
      author='Ahmet Bakan',
      author_email='ahb12@pitt.edu',
      url='http://www.pitt.edu/~ahb12/python/ProDy',
      packages=['prody', 'prody.dynamics', 'prody.proteins'],
      license='GPLv3',
      keywords=('protein, dynamics, elastic network model, gaussian network model, '
                'anisotropic network model, essential dynamics analysis, '
                'principal component analysis, ProteinDataBank, PDB, GNM, ANM, PCA')
      )
