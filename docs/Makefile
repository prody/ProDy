# Makefile for ProDy Documentation

SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILD         = _build
HTML          = $(BUILD)/html
LATEX         = $(BUILD)/latex
WORKDIR       = _workdir
PDFNAME       = ProDy.pdf
PDFPATH       = $(LATEX)/$(PDFNAME)

.PHONY: help clean remove html pdf

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make HTML files of reference manual "
	@echo "  pdf	    to make PDF files of reference manual"
	@echo "  build	    to make HTML and PDF files, and copy PDF into HTML folder"

clean:
	-rm -rf $(BUILD)/*

remove:
	-rm -rf $(HTML)/* $(LATEX)/*


workdir:
	mkdir -p $(WORKDIR)

html: workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILD)/doctrees ../ ../$(HTML)
	@echo
	@echo "Build finished. The HTML pages are in $(BUILD)/html."

pdf: workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) -d ../$(BUILD)/doctrees ../ ../$(LATEX)
	make -C $(BUILD)/latex all-pdf
	@echo
	@echo "PDF file is ready, see $(PDFPATH)."

build: html pdf
	mv -f $(PDFPATH) $(HTML)/$(PDFNAME)