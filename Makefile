# Makefile for some developer commands

REPOPATH = `pwd`

.PHONY: help build build3 remove test

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  build		to build extensions in place"
	@echo "  remove		to remove contributed modules"
	@echo "  test	    to clone, build and test"

build:
	python setup.py build_ext --inplace --force

build3:
	python3 setup.py build_ext --inplace --force

remove:
	rm -f lib/prody/kdtree/*c
	rm -f lib/prody/atomic/pyparsing*py
	rm -f lib/prody/apps/argparse.py
	rm -f lib/prody/proteins/*pairwise2.*

test:
	TMPDIR=`mktemp -d`; REPOPATH=`pwd`; echo $$TMPDIR; cd $$TMPDIR; \
	git clone $$REPOPATH; cd ProDy; \
	python setup.py build_ext --inplace --force; \
	export PYTHONPATH=$$TMPDIR/ProDy/lib/; \
	python scripts/prody test; \
	rm -rf $$TMPDIR

