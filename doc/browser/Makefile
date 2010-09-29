# non-default targets: browser and related

ALL: browser

ice_bib.bib:
	ln -s ../ice_bib.bib ice_bib.bib

references.txt: doxybib.py ice_bib.bib
	./doxybib.py

pism_config.nc:
	ncgen -o pism_config.nc ../../src/pism_config.cdl

pism_config.txt: config_doc.py pism_config.nc
	./config_doc.py > $@

browser : my_doxyfile references.txt pism_config.txt std_names.txt \
		pages.txt bombproof_enth.txt \
                initialization-sequence.dot
	doxygen my_doxyfile

browser.tgz : html/index.html
	tar czvf browser.tgz html

my_doxyfile: doxyfile sed_script.txt
	sed -f sed_script.txt < doxyfile > $@

clean:
	@rm -f references.txt pism_config.txt pism_config.nc texput.* my_doxyfile ice_bib.bib