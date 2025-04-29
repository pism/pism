all: result extra extra_split extra_append snapshots snapshots_split timeseries timeseries_append checkpointing compression grid_mapping beddef_lc

result: test_result.sh
	bash $^

extra: test_extra.sh
	bash $^

extra_split: test_extra.sh
	filename=ex split=yes bash $^

extra_append: test_extra_append.sh
	bash $^

snapshots: test_snapshots.sh
	bash $^

snapshots_split: test_snapshots.sh
	filename=snapshot split=yes bash $^

timeseries: test_timeseries.sh
	bash $^

timeseries_append: test_timeseries_append.sh
	bash $^

checkpointing: test_checkpoints.sh
	bash $^

# miscellaneous:

compression: test_compression.sh
	bash $^
# use "head -100" to stop ncdump from reading doing more than necessary
	ncdump -h -s result_compressed.nc | grep -v pism_config | head -100 | grep _DeflateLevel

grid_mapping: test_grid_mapping.sh
	bash $^

beddef_lc: test_beddef_lc.sh
	bash $^

clean:
	rm -f *~ *.nc
