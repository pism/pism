ALL : installation cheatsheet userman

update: svn_update revclean
ifeq ($(shell touch revision.tex; svnversion | diff revision.tex -),)
	@echo "doc/ directory is up to date."
else
	@echo "Rebuilding manuals ..."
	@svnversion > revision.tex
	$(MAKE) ALL
endif

revision:
	@svnversion > revision.tex

revclean:
	@rm -f revision.tex


# default targets: build PDF manuals

userman:
	$(MAKE) -C userman

installation : installation.tex
	pdflatex installation.tex
	pdflatex installation.tex

cheatsheet : cheatsheet.tex
	pdflatex cheatsheet.tex
	pdflatex cheatsheet.tex


# non-default targets: browser and related

doxybib.txt: doxybib.py ice_bib.bib
	./doxybib.py

config_doc.txt: config_doc.py ../lib/pism_config.nc
	./config_doc.py > config_doc.txt

browser : doxyfile_browser doxybib.txt config_doc.txt std_names.html \
		browser_pages.txt browser_bombproof_enth.txt revision.tex
	doxygen doxyfile_browser

# rebuild browser and then open it with firefox:
browsenow : browser
	firefox doxy/html/index.html &


# utility targets:

.PHONY:  clean revision revclean svn_update userman installation browsenow browser vfigs

clean : revclean
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl
	@rm -f config_doc.txt doxybib.txt

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

