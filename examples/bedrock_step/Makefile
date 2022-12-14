input.nc: bedrock_step.py
	python3 $^ $@

bedrock_step.nc ex.nc ts.nc: input.nc run.sh
	bash -x ./run.sh

figure.png: plot_results.py bedrock_step.nc
	python3 $^ $@

.PHONY: clean
clean:
	rm -f *.nc *.log *~
