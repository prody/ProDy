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
ALLSPHINXOPTS   = -d ../doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS)

.PHONY: help html clean latex pdf archive workdir

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  pdf        to make LaTeX files and run them through pdflatex"

clean:
	-rm -rf $(BUILDDIR)/* $(BUILDDIR)_files.zip $(BUILDDIR)_files.tgz $(BUILDDIR).pdf

workdir:
	mkdir -p $(BUILDDIR)/workdir
	ln -fs ../../$(BASENAME)_files $(BUILDDIR)/workdir/
	ln -fs ../../_static
	ln -fs ../../_templates
	ln -fs ../../_inventory

html: workdir lsfiles
	cd $(BUILDDIR)/workdir; $(SPHINXBUILD) -b html -d ../doctrees ../../ ../html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

latex: workdir lsfiles
	rm -f files.txt;
	cd $(BUILDDIR)/workdir; $(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) ../../ ../latex
	@echo
	@echo "LaTeX file is ready, see $(BUILDDIR)/latex."

pdf: latex
	@echo "Running LaTeX files through pdflatex..."
	make -C $(BUILDDIR)/latex all-pdf
	mv -f $(BUILDDIR)/latex/$(BASENAME).pdf .
	@echo
	@echo "PDF file is ready, see $(BASENAME).pdf."

lsfiles:
	files=`ls $(BASENAME)_files`; \
	if [ -z "$$files" ]; then echo "There are no required files." > files.txt; else \
	ls -hgo --time-style=+ $(BASENAME)_files/ | tail -n +2 | cut -b 14-128 > files.txt; fi

archive:
	rm -f $(BASENAME)_files.zip $(BASENAME)_files.tgz
	files=`ls $(BASENAME)_files`; \
	if [ -n "$$files" ]; then \
	tar -czf $(BASENAME)_files.tgz $(BASENAME)_files/*; \
	zip $(BASENAME)_files.zip $(BASENAME)_files/*; \
	echo "$(BASENAME)_files.zip and $(BASENAME)_files.tgz are created."; \
	else \
	echo "There are no tutorial files."; \
	fi