# Makefile for ProDy Documentation

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILDDIR      = _build
DOCTESTDIR    = _doctest

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS)

.PHONY: help clean html latex pdf doctest copy stats tutorial

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files of entire documentaion"
	@echo "  pdf	    to make PDF files of reference manual and tutorials"
	@echo "  copy	    to copy license, install, and changes to package directory"
	@echo "  stats	    to update download statistics"
	@echo "  doctest    to run all doctests embedded in the documentation"

clean:
	-rm -rf $(BUILDDIR)/*

workdir:
	mkdir -p $(BUILDDIR)/workdir
	for tut in `find tutorials/ -type d -name "*_files"` ; do \
		if [ -d "$$tut" ]; then \
		 	ln -fs ../../$$tut $(BUILDDIR)/workdir/;\
		fi; \
	done

html: workdir
	rm -f _inventory/prody*inv
	cd $(BUILDDIR)/workdir; $(SPHINXBUILD) -b html -d $(BUILDDIR)/doctrees ../../ ../html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf:
	cd reference; make pdf
	cd tutorials; \
	for tut in `find . -type d` ; do \
		if [ -d "$$tut/$$tut"_files ]; then \
		 	cd $$tut; make archive html pdf; cd ..;\
		fi; \
	done

copy:
	cp -f changes.rst ../CHANGES.rst
	cp -f getprody.rst ../INSTALL.rst
	cp -f license.rst ../LICENSE.rst
	@echo
	@echo "Files copied."

stats:
	pyps update -s ProDy_stats.pkl ProDy; \
	pyps latest -o _static/ProDy_downloads.csv -d , --rst --no-md5 ProDy; \
	pyps monthly -p _static/ProDy_downloads.png ProDy_stats.pkl
	@echo
	@echo "Download statistics updated."

doctest:
	cd $(DOCTESTDIR); \
	$(SPHINXBUILD) -b doctest -d ../$(BUILDDIR)/doctrees/ .. .
	@echo
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(DOCTESTDIR)."
