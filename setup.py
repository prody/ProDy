from distutils.core import setup
setup(name='ProDy',
      version='0.1.1',
      author='Ahmet Bakan',
      author_email='ahb12@pitt.edu',
      description='A Python Package for Protein Structural Dynamics Analysis',
      long_description=('ProDy is a Python package for protein structural dynamics analysis. '
                'ProDy features a fast PDB parser and a powerful atom seletion engine. '
                'ProDy can be used for elastic network model and principal component analysis. '
                'ProDy comes with many handy functions for comparative analysis and plotting.'),
      url='http://www.pitt.edu/~ahb12/python/ProDy',
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
      provides=['ProDy(0.1.1)']
      )
