#!/bin/sh

# A script for testing ProDy repository clone
# This script can be used to discover files missing in the repository

# Set paths and create a temporary folder
REPOPATH=`pwd`
REPONAME=${PWD##*/}
TMPDIR=`mktemp -d`

# Clone ProDy
cd $TMPDIR
git clone $REPOPATH

# Build and copy
cd $REPONAME
python setup.py build_ext --inplace --force

# Run tests
PYPATH=$PYTHONPATH
export PYTHONPATH=$TMPDIR/$REPONAME/lib/:$PYPATH
echo $PYTHONPATH
python scripts/prody test

# Restore PYTHONPATH and remove temporary files
export PYTHONPATH=$PYPATH
rm -rf $TMPDIR
