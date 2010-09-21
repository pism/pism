# non-default targets: browser and related

revision.tex:
	svnversion > revision.tex
.PHONY: revision.tex

doxybib.txt: doxybib.py ../ice_bib.bib
	./doxybib.py

pism_config.nc:
	ncgen -o pism_config.nc ../../src/pism_config.cdl

config_doc.txt: config_doc.py pism_config.nc
	./config_doc.py > config_doc.txt

browser : doxyfile_browser doxybib.txt config_doc.txt std_names.html \
		browser_pages.txt browser_bombproof_enth.txt revision.tex \
                initialization-sequence.dot
	doxygen doxyfile_browser

browser.tgz : doxy/html/index.html
	cd doxy/html && tar czvf ../../browser.tgz *

clean:
	@rm -f config_doc.txt doxybib.txt pism_config.nc revision.tex texput.*