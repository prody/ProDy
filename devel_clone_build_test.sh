# A script for testing ProDy
# Use this script to discover missing files in
# the repository

# Set paths and create a temporary folder
REPOPATH=`pwd`
REPONAME=${PWD##*/}
TEMPPATH=`mktemp -d`

# Clone ProDy
cd $TEMPPATH
git clone $REPOPATH

# Build and copy
cd $REPONAME
python setup.py build
python setup.py copy

# Run tests
PYTHONPATHBACKUP=$PYTHONPATH
export PYTHONPATH=$TEMPPATH/$REPONAME/lib/
echo $PYTHONPATH
python scripts/prody test

# Restore PYTHONPATH and remove temporary files
export PYTHONPATH=$PYTHONPATHBACKUP
rm -rf $TEMPPATH

