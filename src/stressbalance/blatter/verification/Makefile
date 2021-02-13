all: manufactured_solutions.cc manufactured_solutions.hh

%.cc: %.py
	python3 $^ > $@

%.hh: %.py
	python3 $^ --header > $@

.PHONY: clean
clean:
	rm -f manufactured_solutions.cc manufactured_solutions.hh
