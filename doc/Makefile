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

browser :
	$(MAKE) -C browser browser

browser.tgz :
	$(MAKE) -C browser browser.tgz

# rebuild browser and then open it with firefox:
browsenow : browser
	firefox browser/doxy/html/index.html &

# utility targets:

.PHONY:  clean revision revclean svn_update userman installation browsenow browser vfigs

clean : revclean
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl *.rel *.synctex.gz
	@$(MAKE) -C userman clean
	@$(MAKE) -C browser clean

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

