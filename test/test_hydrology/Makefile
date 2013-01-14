all:

.PHONY: clean

clean:
	@rm -f *~ *.so exactP.c *.pyc PISMNC.py runP*.txt
	@rm -rf src/
	@rm -rf build/

exactP:
	python setup.py build_ext --inplace
