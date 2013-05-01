# Makefile for tutorials
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILDDIR      = _build
DOCTESTDIR    = ../../_doctest
BASENAME      = $${PWD\#\#*/}

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
	rm -f files.txt;
	if [ -z `ls enm_analysis_files/` ]; then echo "No required files." > files.txt; else \
	ls -hgo --time-style=+ $(BASENAME)_files/ | tail -n +2 | cut -b 14-128 > files.txt; fi
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	make -C $(BUILDDIR)/latex all-pdf
	mv -f $(BUILDDIR)/latex/$(BASENAME).pdf .
	@echo
	@echo "PDF file is ready, see $(BASENAME).pdf."

files:
	rm -f $(BASENAME)_files.zip $(BASENAME)_files.tgz
	tar -czf $(BASENAME)_files.tgz $(BASENAME)_files/*  
	zip $(BASENAME)_files.zip $(BASENAME)_files/*
	@echo 
	@echo "$(BASENAME)_files.zip and $(BASENAME)_files.tgz are created."


test:
	cp -f files/* $(DOCTESTDIR); cd $(DOCTESTDIR); \
	$(SPHINXBUILD) -b doctest -d $(PWD)/$(BUILDDIR)/doctrees/ $(PWD) .
	mv -f $(DOCTESTDIR)/output.txt .
	@echo
	@echo "Testing of tutorial finished, see results in output.txt."
