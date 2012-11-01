# A script for testing ProDy
# Run this script before committing your changes

# Create a temporary folder
TMPDIR=`mktemp -d`
rm -rf build

# Build and install into temporary folder
python setup.py build
python setup.py install --install-lib=$TMPDIR --install-scripts=$TMPDIR/scripts

# Run tests
PYPATH=$PYTHONPATH
export PYTHONPATH=$TMPDIR
echo $PYTHONPATH
python scripts/prody test

# Restore PYTHONPATH and remove temporary files
export PYTHONPATH=$PYPATH
rm -rf $TMPDIR
