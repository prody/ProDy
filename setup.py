# setup.py
# This file is retained for building C/C++ extensions.
# All project metadata has been moved to pyproject.toml.

import os
import platform
import shutil
import sys
from glob import glob
from os.path import isfile, join

from setuptools import Extension, setup

# The build system (pip, build) ensures numpy is available because it's
# listed in `build-system.requires` in pyproject.toml.
import numpy

# --- Logic to copy pre-compiled .so file ---
# Note: Distributing pre-compiled binaries like this is not standard practice.
# It would be better to compile this from source as part of the build.
# However, this preserves the original behavior of your script.
hpbSoDir = join('prody', 'proteins', 'hpbmodule',
                f'hpb_Python{sys.version_info[0]}.{sys.version_info[1]}')
proteinsDir = join('prody', 'proteins')

if not os.path.exists(join(proteinsDir, 'hpb.so')):
    try:
        shutil.copy(join(hpbSoDir, "hpb.so"), proteinsDir)
    except FileNotFoundError:
        # It's okay if the precompiled file doesn't exist for this platform.
        pass

# --- Define C extensions ---
EXTENSIONS = [
    Extension('prody.dynamics.rtbtools',
              glob(join('prody', 'dynamics', 'rtbtools.c')),
              include_dirs=[numpy.get_include()]),
    Extension('prody.dynamics.smtools',
              glob(join('prody', 'dynamics', 'smtools.c')),
              include_dirs=[numpy.get_include()]),
    Extension('prody.sequence.msatools',
              [join('prody', 'sequence', 'msatools.c')],
              include_dirs=[numpy.get_include()]),
    Extension('prody.sequence.msaio',
              [join('prody', 'sequence', 'msaio.c')],
              include_dirs=[numpy.get_include()]),
    Extension('prody.sequence.seqtools',
              [join('prody', 'sequence', 'seqtools.c')],
              include_dirs=[numpy.get_include()]),
]

# --- Define C++ extensions and platform-specific build args ---
tntDir = join('prody', 'utilities', 'tnt')
extra_compile_args = []

if platform.system() == 'Darwin':
    # Set environment variables for macOS clang compiler
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = platform.mac_ver()[0]
    os.environ['CC'] = 'clang'
    os.environ['CXX'] = 'clang++'
    # extra_compile_args.append('-stdlib=libc++') # May be needed for some setups

CONTRIBUTED = [
    Extension('prody.kdtree._CKDTree',
              [join('prody', 'kdtree', 'KDTree.c'),
               join('prody', 'kdtree', 'KDTreemodule.c')],
              include_dirs=[numpy.get_include()]),
    Extension('prody.proteins.ccealign',
              [join('prody', 'proteins', 'ccealign', 'ccealignmodule.cpp')],
              include_dirs=[tntDir],
              language='c++',
              extra_compile_args=extra_compile_args),
]

# Conditionally add contributed extensions if their source files exist
for ext in CONTRIBUTED:
    if all(isfile(src) for src in ext.sources):
        EXTENSIONS.append(ext)

# --- Setup call ---
# This call is now minimal. It only provides the extension modules to setuptools.
# All other configuration is automatically picked up from pyproject.toml.
if __name__ == "__main__":
    setup(
        ext_modules=EXTENSIONS,
    )
