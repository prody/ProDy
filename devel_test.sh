#!/bin/sh

# rebuild extensions in place and run fast tests
python setup.py build_ext --inplace --force
python scripts/prody test
