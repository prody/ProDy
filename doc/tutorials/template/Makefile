# Makefile for ProDy Tutorial

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILDDIR      = _build
DOCTESTDIR    = ../../_doctest
BASENAME      = $${PWD\#\#*/}
DOCDIR        = ../../../Documentation/tutorials/$(BASENAME)

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean pdf files test 

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  pdf        to make LaTeX files and run them through pdflatex"
	@echo "  test       to run all doctests embedded in the documentation"

clean:
	-rm -rf $(BUILDDIR)/*

pdf:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	make -C $(BUILDDIR)/latex all-pdf
	mv -f $(BUILDDIR)/latex/$(BASENAME).pdf $(DOCDIR)
	@echo
	@echo "pdflatex finished; the PDF files are in $(DOCDIR)"

files:
	tar -czf ../../../Documentation/tutorials/$(BASENAME)/$(BASENAME).tgz files/*  
	zip ../../../Documentation/tutorials/$(BASENAME)/$(BASENAME).zip files/*
	@echo 
	@echo "$(BASENAME).zip and $(BASENAME).tgz archives are are copied."


test:
	cp -f files/* $(DOCTESTDIR); cd $(DOCTESTDIR); \
	$(SPHINXBUILD) -b doctest -d $(PWD)/$(BUILDDIR)/doctrees/ $(PWD) .
	mv -f $(DOCTESTDIR)/output.txt test_output.txt
	@echo
	@echo "Testing of doctests in the sources finished, look at the" \
	      "results in test_output.txt."
