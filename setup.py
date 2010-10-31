from distutils.core import setup
setup(name='prody',
      version='0.0.2',
      description='Protein Structure and Dynamics',
      long_description='A Python Package for Protein Structural Dynamics Analysis',
      author='Ahmet Bakan',
      author_email='ahb12@pitt.edu',
      url='http://www.pitt.edu/~ahb12/python/prody',
      packages=['prody', 'prody.dynamics', 'prody.proteins'],
      license='GPL',
      keywords=('protein, dynamics, elastic network model, gaussian network model'
                'anisotropic network model, essential dynamics analysis, '
                'principal component analysis, ProteinDataBank, PDB', 
                'GNM', 'ANM', 'PCA')
      )
