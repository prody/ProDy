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
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html pdf test copy stats tutorial

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  pdf	    to make LaTeX files and run them through pdflatex"
	@echo "  test	    to run all doctests embedded in the documentation"
	@echo "  copy	    to copy license, install, and changes to package directory"
	@echo "  stats	    to update download statistics"
	@echo "  tutorial   to prepare tutorial PDF and achieve files"

clean:
	-rm -rf $(BUILDDIR)/*

html:
	$(SPHINXBUILD) -b html -d $(BUILDDIR)/doctrees . $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) $(BUILDDIR)/latex
	@echo "Running LaTeX files through pdflatex..."
	make -C $(BUILDDIR)/latex all-pdf
	#mv -f $(BUILDDIR)/latex/ProDy.pdf ../Documentation
	@echo
	@echo "pdflatex finished; the PDF files are in $(BUILDDIR)/latex."

test:
	cd $(DOCTESTDIR); \
	$(SPHINXBUILD) -b doctest -d ../$(BUILDDIR)/doctrees/ .. . 
	@echo
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(DOCTESTDIR)."

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

tutorial:
	cd tutorials; \
	for tut in `find . -type d` ; do \
		if [ -d "$$tut/$$tut"_files ]; then \
		 	cd $$tut; \
		 	make pdf; \
		 	make files; \
		 	cd ..; \
		fi; \
	done
