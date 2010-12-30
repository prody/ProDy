readme = open('README.txt')
long_description = ''.join(readme.readlines())
readme.close()

from distutils.core import setup
setup(name='ProDy',
      version='0.5.1',
      author='Ahmet Bakan',
      author_email='ahb12 at pitt dot edu',
      description='A Python Package for Protein Dynamics Analysis',
      long_description=long_description,
      url='http://www.csb.pitt.edu/People/abakan/software/ProDy',
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
      requires=['NumPy', ],
      provides=['ProDy(0.5.1)']
      )
