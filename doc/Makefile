# Makefile for ProDy Documentation

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
PAPER         = letter
BUILDDIR      = _build
WORKDIR       = _workdir

.PHONY: help clean html pdf copy stats

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files of entire documentaion"
	@echo "  pdf	    to make PDF files of reference manual and tutorials"
	@echo "  copy	    to copy license, install, and changes to package directory"
	@echo "  stats	    to update download statistics"

clean:
	-rm -rf $(BUILDDIR)/*

workdir:
	mkdir -p $(WORKDIR)
	for tut in `find tutorials/ -type d -name "*_files"` ; do \
		if [ -d "$$tut" ]; then \
		 	ln -fs ../../$$tut $(WORKDIR)/;\
		fi; \
	done

html: workdir
	cd $(WORKDIR); $(SPHINXBUILD) -b html -d ../$(BUILDDIR)/doctrees ../ ../$(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

pdf:
	rm -f reference/*pdf tutorials/*/*pdf tutorials/*/*files.zip tutorials/*/*files.tgz _inventory/prody*inv
	cd reference; make pdf
	cd tutorials; \
	for tut in `find . -type d`; do \
		if [ -d "$$tut/$$tut"_files ]; then \
		 	cd $$tut; make archive pdf; /bin/cp -rf $(WORKDIR)/* ../../$(WORKDIR); cd ..; \
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
