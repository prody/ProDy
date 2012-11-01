# A script for testing ProDy
# Use this script to discover missing files in
# the repository

# Set paths and create a temporary folder
REPOPATH=`pwd`
REPONAME=${PWD##*/}
TMPDIR=`mktemp -d`

# Clone ProDy
cd $TMPDIR
git clone $REPOPATH

# Build and copy
cd $REPONAME
python setup.py build
python setup.py copy

# Run tests
PYPATH=$PYTHONPATH
export PYTHONPATH=$TMPDIR/$REPONAME/lib/
echo $PYTHONPATH
python scripts/prody test

# Restore PYTHONPATH and remove temporary files
export PYTHONPATH=$PYPATH
rm -rf $TMPDIR

