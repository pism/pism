ALL : installation userman browser

update: svn_update revclean
ifeq ($(shell touch revision.tex; svnversion | diff revision.tex -),)
	@echo "doc/ directory is up to date."
else
	@echo "Rebuilding manuals ..."
	@svnversion > revision.tex
	$(MAKE) ALL
endif

userman:
	$(MAKE) -C userman

revision:
	@svnversion > revision.tex

revclean:
	@rm -f revision.tex

# default manuals and browser:

installation : installation.tex
	pdflatex installation.tex
	pdflatex installation.tex

browser : doxyfile_browser doxybib.txt config_doc.txt std_names.html \
		browser_pages.txt browser_bombproof_enth.txt revision
	doxygen doxyfile_browser

doxybib.txt: doxybib.py ice_bib.bib
	./doxybib.py

config_doc.txt: config_doc.py ../lib/pism_config.nc
	./config_doc.py > config_doc.txt

# non-default targets; sometimes useful:
#   browsenow      rebuild browser and then open it with firefox
#   vfigs          create PDF figures from a test/verifynow.py output

browsenow : browser
	firefox doxy/html/index.html &


PISMREV_VFIGS ?= 589 # revision stamp for verification results in User's Manual figures;
                     # generally not the current trunk revision
VFNOW_OUT=vfnow_r589_L5_5mar09.txt
# list old ones here; keep them, for comparison:
#VFNOW_OUT=vfnow_r327_L5_23jun08.txt
#VFNOW_OUT=vfnow_r193_L5_10oct07.txt
vfigs : verif/$(VFNOW_OUT)
	# producing figures thickerrs_G.pdf, temperrs_G.df, etc. 
	# in subdirectory figs/ for use in manual
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t G -e 'geometry' -x 'prcntVOLrelmaxETAcenterH' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/thickerrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t G -e 'temp' -x 'basecenterT' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/temperrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t G -e 'surf vels' -x 'maxUvecmaxW' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/surfvelerrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t I -e 'velocity' -x 'maxvectorprcntavvecavvector' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/velerrs
	../test/vnreport.py -f verif/$(VFNOW_OUT) -t K -e 'temp' \
	   -r ${PISMREV_VFIGS} -g pdf -o figs/temperrs

.PHONY:  clean revision revclean svn_update userman installation browsenow browser vfigs

clean : revclean
	@rm -f *.aux *.blg *.idx *.ind *.ilg *.log *.out *.toc texput.* *.bbl
	@rm -f config_doc.txt doxybib.txt

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

