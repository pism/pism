all: manufactured_solutions.cc manufactured_solutions.hh

%.cc: %.py
	python $^ > $@

%.hh: %.py
	python $^ --header > $@
