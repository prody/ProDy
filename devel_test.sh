#!/bin/sh

# rebuild extensions in place and run fast tests
rm -f lib/prody/*/*.so
rm -f lib/prody/*/*.pyd
python setup.py build_ext --inplace
python scripts/prody test
