# Fetch data from AWS

grid_file=grid_RGI2000-v7.0-C-08-00658.nc
climate_file=carra2_RGI2000-v7.0-C-08-00658.nc
climate_mean_file=carra2_monthly_mean_RGI2000-v7.0-C-08-00658.nc
rgi_file=rgi_RGI2000-v7.0-C-08-00658-G.gpkg

for f in "$grid_file" "$climate_file" "$climate_mean_file" "$rgi_file"; do
    wget -nc "https://pism-cloud-data.s3.amazonaws.com/pism/examples/storglaciaren/$f"
done

python create_storglaciaren_res.py --grid "$grid_file" --dem-source arcticdem

