# Makefile for ProDy Documentation

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILDDIR      = _build
WORKDIR       = _workdir
PDFNAME       = ProDy.pdf

.PHONY: help clean remove html pdf

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files of entire documentaion"
	@echo "  pdf	    to make PDF files of reference manual and tutorials"
	@echo "  copy	    to copy license, install, and changes to package directory"
	@echo "  stats	    to update download statistics"

clean:
	-rm -rf $(BUILDDIR)/*

remove:
	-rm -rf $(BUILDDIR)/html*


workdir:
	mkdir -p $(WORKDIR)
	for tut in `find tutorials/ -type d -name "*_files"` ; do \
		if [ -d "$$tut" ]; then \
		 	ln -fs ../$$tut $(WORKDIR)/;\
		fi; \
	done

html: workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf: workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/latex
	make -C $(BUILDDIR)/latex all-pdf
	mv -f _build/latex/$(PDFNAME) .
	@echo
	@echo "PDF file is ready, see $(PDFNAME)."