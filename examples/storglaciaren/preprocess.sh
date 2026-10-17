# Fetch data from AWS

grid_file=grid_RGI2000-v7.0-C-08-00658.nc
wget -nc https://pism-cloud-data.s3.amazonaws.com/pism/examples/storglaciaren/$grid_file

python create_storglaciaren_res.py --grid $grid_file --dem-source arcticdem

