all: blatter_mms.cc blatter_mms.hh

blatter_mms.cc: manufactured_solutions.py
	python $^ > $@

blatter_mms.hh: manufactured_solutions.py
	python $^ --header > $@
