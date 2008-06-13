/* commonNetCDF.c
 *
 * Nathan Shemonski
 * 08/01/07
 */

void write_to_net_cdf(Projection p, float **b, float **H, float **up, float **a, char *file_name) {
  int stat, ncid;
  int x_dim, y_dim, t_dim;
  int x_dims[1], y_dims[1], t_dims[1];
  int lat_dims[2], lon_dims[2], ac_dims[2], bed_dims[2], H_dims[2], up_dims[2];
  int x_id, y_id, t_id/*, p_id*/;
  int lat_id, lon_id, ac_id, bed_id, H_id, up_id;
  int i;

  /* write to netCDF */
  printf("Writing file %s...", file_name);
  stat = nc_create(file_name, NC_CLOBBER|NC_WRITE, &ncid);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_def_dim(ncid, "x", p.n, &x_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "y", p.m, &y_dim);
  check_err(stat, __LINE__, __FILE__);
  stat = nc_def_dim(ncid, "t", NC_UNLIMITED, &t_dim);
  check_err(stat, __LINE__, __FILE__);

  x_dims[0] = x_dim;
  stat = nc_def_var(ncid, "x", NC_FLOAT, 1, x_dims, &x_id);
  stat = nc_put_att_text(ncid, x_id, "long_name", strlen("x-coordinate in Cartesian system"), "x-coordinate in Cartesian system");
  stat = nc_put_att_text(ncid, x_id, "standard_name", strlen("projection_x_coordinate"), "projection_x_coordinate");
  stat = nc_put_att_text(ncid, x_id, "units", strlen("m"), "m");
  check_err(stat,__LINE__,__FILE__);

  y_dims[0] = y_dim;
  stat = nc_def_var(ncid, "y", NC_FLOAT, 1, y_dims, &y_id);
  stat = nc_put_att_text(ncid, y_id, "long_name", strlen("y-coordinate in Cartesian system"), "y-coordinate in Cartesian system");
  stat = nc_put_att_text(ncid, y_id, "standard_name", strlen("projection_y_coordinate"), "projection_y_coordinate");
  stat = nc_put_att_text(ncid, y_id, "units", strlen("m"), "m");
  check_err(stat,__LINE__,__FILE__);
  
  t_dims[0] = t_dim;
  stat = nc_def_var(ncid, "t", NC_FLOAT, 1, t_dims, &t_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text(ncid, t_id, "units", strlen("seconds since 2007-01-01 00:00:00"), "seconds since 2007-01-01 00:00:00");

  lon_dims[0] = x_dim;
  lon_dims[1] = y_dim;
  stat = nc_def_var(ncid, "lon", NC_FLOAT, 2, lon_dims, &lon_id);
  stat = nc_put_att_text(ncid, lon_id, "long_name", strlen("longitude"), "longitude");
  stat = nc_put_att_text(ncid, lon_id, "standard_name", strlen("longitude"), "longitude");
  stat = nc_put_att_text(ncid, lon_id, "units", strlen("degrees_east"), "degrees_east");
  check_err(stat,__LINE__,__FILE__);

  lat_dims[0] = x_dim;
  lat_dims[1] = y_dim;
  stat = nc_def_var(ncid, "lat", NC_FLOAT, 2, lat_dims, &lat_id);
  stat = nc_put_att_text(ncid, lat_id, "long_name", strlen("latitude"), "latitude");
  stat = nc_put_att_text(ncid, lat_id, "standard_name", strlen("latitude"), "latitude");
  stat = nc_put_att_text(ncid, lat_id, "units", strlen("degrees_north"), "degrees_north");
  check_err(stat,__LINE__,__FILE__);

  if (a != NULL) {
    ac_dims[0] = x_dim;
    ac_dims[1] = y_dim;
    stat = nc_def_var(ncid, "accum", NC_FLOAT, 2, ac_dims, &ac_id);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_float(ncid, ac_id, "missing_value", NC_FLOAT, 1, &missing_value);
    stat = nc_put_att_text(ncid, ac_id, "long_name", strlen("mean ice equivalent accumulation rate"), "mean ice equivalent accumulation rate");
    stat = nc_put_att_text(ncid, ac_id, "standard_name", strlen("land_ice_surface_specific_mass_balance"), "land_ice_surface_specific_mass_balance");
    stat = nc_put_att_text(ncid, ac_id, "units", strlen("m s-1"), "m s-1");
    check_err(stat,__LINE__,__FILE__);
  }

  if (b != NULL) {
    bed_dims[0] = x_dim;
    bed_dims[1] = y_dim;
    stat = nc_def_var(ncid, "b", NC_FLOAT, 2, bed_dims, &bed_id);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_float(ncid, bed_id, "missing_value", NC_FLOAT, 1, &missing_value);
    stat = nc_put_att_text(ncid, bed_id, "long_name", strlen("bedrock_altitude"), "bedrock_altitude");
    stat = nc_put_att_text(ncid, bed_id, "standard_name", strlen("bedrock_altitude"), "bedrock_altitude");
    stat = nc_put_att_text(ncid, bed_id, "units", strlen("m"), "m");
    check_err(stat,__LINE__,__FILE__);
  }

  if (H != NULL) {
    H_dims[0] = x_dim;
    H_dims[1] = y_dim;
    stat = nc_def_var(ncid, "H", NC_FLOAT, 2, H_dims, &H_id);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_float(ncid, H_id, "missing_value", NC_FLOAT, 1, &missing_value);
    stat = nc_put_att_text(ncid, H_id, "long_name", strlen("land_ice_thickness"), "land_ice_thickness");
    stat = nc_put_att_text(ncid, H_id, "standard_name", strlen("land_ice_thickness"), "land_ice_thickness");
    stat = nc_put_att_text(ncid, H_id, "units", strlen("m"), "m");
    check_err(stat,__LINE__,__FILE__);
  }

  if (up != NULL) {
    up_dims[0] = x_dim;
    up_dims[1] = y_dim;
    stat = nc_def_var(ncid, "dbdt", NC_FLOAT, 2, up_dims, &up_id);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_float(ncid, up_id, "missing_value", NC_FLOAT, 1, &missing_value);
    stat = nc_put_att_text(ncid, up_id, "long_name", strlen("uplift_rate"), "uplift_rate");
    stat = nc_put_att_text(ncid, up_id, "standard_name", strlen("tendency_of_bedrock_altitude"), "tendency_of_bedrock_altitude");
    stat = nc_put_att_text(ncid, up_id, "units", strlen("m s-1"), "m s-1");
    check_err(stat,__LINE__,__FILE__);
  }

  nc_enddef(ncid);

  float x_val[p.n];
  for (i=0; i<p.n; ++i) {
    x_val[i] = (i-p.n/2)*5000;
  }
  stat = nc_put_var_float(ncid, x_id, x_val);
  check_err(stat,__LINE__,__FILE__);

  float y_val[p.m];
  for (i=0; i<p.m; ++i) {
    y_val[i] = (i-p.m/2)*5000;
  }
  stat = nc_put_var_float(ncid, y_id, y_val);
  check_err(stat,__LINE__,__FILE__);

  float lat_f[p.m*p.n];
  float lon_f[p.m*p.n];
  for (i=0; i<p.m*p.n; ++i) {
    lat_f[i] = p.lat[i/p.m][i%p.m];
    lon_f[i] = p.lon[i/p.m][i%p.m];
  }
  stat = nc_put_var_float(ncid, lat_id, lat_f);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_var_float(ncid, lon_id, lon_f);
  check_err(stat,__LINE__,__FILE__);

  if (a != NULL) {
    float a_f[p.m*p.n];
    for (i=0; i<p.m*p.n; ++i) {
      a_f[i] = a[i/p.m][i%p.m];
    }
    stat = nc_put_var_float(ncid, ac_id, a_f);
    check_err(stat,__LINE__,__FILE__);
  }

  if (b != NULL) {
    float b_f[p.m*p.n];
    for (i=0; i<p.m*p.n; ++i) {
      b_f[i] = b[i/p.m][i%p.m];
    }
    stat = nc_put_var_float(ncid, bed_id, b_f);
    check_err(stat,__LINE__,__FILE__);
  }

  if (H != NULL) {
    float H_f[p.m*p.n];
    for (i=0; i<p.m*p.n; ++i) {
      H_f[i] = H[i/p.m][i%p.m];
    }
    stat = nc_put_var_float(ncid, H_id, H_f);
    check_err(stat,__LINE__,__FILE__);
  }

  if (up != NULL) {
    float up_f[p.m*p.n];
    for (i=0; i<p.m*p.n; ++i) {
      up_f[i] = up[i/p.m][i%p.m];
    }
    stat = nc_put_var_float(ncid, up_id, up_f);
    check_err(stat,__LINE__,__FILE__);
  }

  nc_close(ncid);
  printf("done.\n");
}

void write_NSIDC_to_net_cdf(Projection NSIDC_proj, float **bed_raw, float **H_raw) {
  int stat, p_id, ncid; 
  float tmp;

  write_to_net_cdf(NSIDC_proj, bed_raw, H_raw, NULL, NULL, "NSIDC_unproj.nc");
  stat = nc_open("NSIDC_unproj.nc", NC_WRITE, &ncid);
  check_err(stat,__LINE__,__FILE__);

  nc_redef(ncid);

  stat = nc_def_var(ncid, "polar_stereographic", NC_INT, 0, NULL, &p_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text(ncid, p_id, "grid_mapping_name", strlen("polar_stereographic"), "polar_stereographic");
  check_err(stat,__LINE__,__FILE__);
  tmp = -39.;
  stat = nc_put_att_float(ncid, p_id, "straight_vertical_longitude_from_pole", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 90.;
  stat = nc_put_att_float(ncid, p_id, "latitude_of_projection_origin", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 71.;
  stat = nc_put_att_float(ncid, p_id, "standard_parallel", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  
  nc_close(ncid);

}

void write_PARCA_to_net_cdf(Projection PARCA_proj, float **accum_raw) {
  int utm_id, stat, ncid;
  float tmp;

  write_to_net_cdf(PARCA_proj, NULL, NULL, NULL, accum_raw, "PARCA_unproj.nc");
  stat = nc_open("PARCA_unproj.nc", NC_WRITE, &ncid);
  check_err(stat,__LINE__,__FILE__);

  nc_redef(ncid);

  stat = nc_def_var(ncid, "transverse_mercator", NC_INT, 0, NULL, &utm_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text(ncid, utm_id, "grid_mapping_name", strlen("transverse_mercator"), "transverse_mercator");
  check_err(stat,__LINE__,__FILE__);
  tmp = 0.999600;
  stat = nc_put_att_float(ncid, utm_id, "scale_factor_at_central_meridian", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 0.;
  stat = nc_put_att_float(ncid, utm_id, "latitude_of_projection_origin", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = -39.;
  stat = nc_put_att_float(ncid, utm_id, "longitude_of_central_meridian", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 500000.;
  stat = nc_put_att_float(ncid, utm_id, "false_easting", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 0.;
  stat = nc_put_att_float(ncid, utm_id, "false_northing", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);

  nc_close(ncid);

}

void write_commo_to_net_cdf(Projection commo_proj, float **accum_com, float **bed_com, float **H_com, float **up_com) {
  int stat, p_id, ncid;
  float tmp;

  write_to_net_cdf(commo_proj, bed_com, H_com, up_com, accum_com, "common.nc");
  stat = nc_open("common.nc", NC_WRITE, &ncid);
  check_err(stat,__LINE__,__FILE__);

  nc_redef(ncid);

  stat = nc_def_var(ncid, "polar_stereographic", NC_INT, 0, NULL, &p_id);
  check_err(stat,__LINE__,__FILE__);
  stat = nc_put_att_text(ncid, p_id, "grid_mapping_name", strlen("polar_stereographic"), "polar_stereographic");
  check_err(stat,__LINE__,__FILE__);
  tmp = -40.;
  stat = nc_put_att_float(ncid, p_id, "straight_vertical_longitude_from_pole", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 72.;
  stat = nc_put_att_float(ncid, p_id, "latitude_of_projection_origin", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  tmp = 85.;
  stat = nc_put_att_float(ncid, p_id, "standard_parallel", NC_FLOAT, 1, &tmp);
  check_err(stat,__LINE__,__FILE__);
  
  nc_close(ncid);

}
