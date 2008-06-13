/* commonProjection.c
 *
 * Nathan Shemonski
 * modeled after Ed Bueler's code
 * 07/24/07
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lib_proj.h>
#include <netcdf.h>
#include <string.h>

typedef struct {
  float **lat;
  float **lon;
  int m;
  int n;
  float yll;
  float xll;
  PJ *ref;
} Projection;

#define secpera 3.1556926e7
const int GRID_SIZE = 5000;
float missing_value = -9999.;

#include "commonNetCDF.c"

int check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    exit(1);
  }
  return 0;
}

void get_projection(char **parms, Projection *pro, float grid_size) {
  PJ *ref;
  XY idata;
  LP odata;
  int i,j;
  float xll, yll;

  xll = pro->xll;
  yll = pro->yll;

  if ( ! (ref = pj_init(5/*sizeof(parms)/sizeof(char *)*/, parms)) ) {
    fprintf(stderr, "Projection initialization failed\n");
    exit(1);
  }

  pro->ref = ref;

  for (i=0; i<pro->n; ++i) {
    for (j=0; j<pro->m; j++) {
      idata.y = grid_size * j + yll;
      idata.x = grid_size * i + xll;
      odata = pj_inv(idata, ref);
      pro->lat[i][j] = odata.phi * RAD_TO_DEG;
      pro->lon[i][j] = odata.lam * RAD_TO_DEG;
      if ((pro->lat[i][j] > 90.0) || (pro->lat[i][j] < -90.0)) {
        printf("Unreasonable value for latitude in writing lat.dat,lon.dat.  Exiting.\n");
        exit(1);
      }
      if ((pro->lon[i][j] > 180.0) || (pro->lon[i][j] < -180.0)) {
        printf("Unreasonable value for longitude in writing lat.dat,lon.dat.  Exiting.\n");
        exit(1);
      }
    }
  }
}

void get_NSIDC_projection(Projection *proj) {
  static char *params[] = {"proj=stere", "ellps=WGS84", "lat_0=90", "lon_0=-39", "lat_ts=71"};
  proj->xll = -800000.;
  proj->yll = -3400000.;

  get_projection(params, proj, GRID_SIZE);
}

void get_PARCA_projection(Projection *proj) {
  //static char *params[] = {"proj=utm", "ellps=WGS84", "k_0=0.9996", "lon_0=-39", "lat_0=0"};
  static char *params[] = {"proj=utm", "ellps=WGS84", "k_0=0.9996", "lon_0=-39", "lat_0=0"};
  proj->xll = -302499.05169293;
  proj->yll = 6627499.8615037;

  get_projection(params, proj, GRID_SIZE);
}

void get_commo_projection(Projection *proj) {
  static char *params[] = {"proj=stere", "ellps=WGS84", "lat_0=72", "lon_0=-40", "lat_ts=85"};
  proj->xll = -900000.;
  proj->yll = -1400000.;

  get_projection(params, proj, GRID_SIZE);

  if (proj->n%2 == 0) {
    proj->n--;
  }
  if (proj->m%2 == 0) {
    proj->m--;
  }
}

void init_proj(Projection *proj) {
  int i;
  
  proj->lat = (float **)malloc(sizeof(float *) * proj->n);
  proj->lon = (float **)malloc(sizeof(float *) * proj->n);
  for (i=0; i<proj->n; ++i) {
    proj->lat[i] = (float *)malloc(sizeof(float) * proj->m);
    proj->lon[i] = (float *)malloc(sizeof(float) * proj->m);
  }
}

void read_PARCA_data(float **data, int m, int n, char *file_name) {
  char tmp_buff[32];
  float tmp_num;
  int i, j;
  FILE *file;
  
  file = fopen (file_name, "r");

  for (i=0; i<6; ++i) {
    fscanf(file, "%s%f\n", tmp_buff, &tmp_num);
  }

  for (j=0; j<m; j++) {
    for (i=0; i<n; ++i) {
      fscanf(file,"%f",&data[i][m-j-1]);
      data[i][m-j-1] /*g/(cm^2 yr)*/ *= (pow(100.,2.) /*cm^2/m^2*/ ) / ( 1000. /*g/kg*/ * 910. /*kg/m^3*/ * secpera /*s/yr*/);
    }
  }

}

void read_NSIDC_data(float **data, int m, int n, char *file_name) {
  char tmp_buff[32];
  float tmp_num;
  int i, j;
  FILE *file;
  
  file = fopen (file_name, "r");

  for (j=0; j<m; j++) {
    for (i=0; i<n; ++i) {
      fscanf(file,"%f",&data[i][j]);
    }
  }

}

void read_drad_data(float **data, int n, int m, char *file_name, char *var_name) {
  int stat, ncid, i, j;
  int v_id, n_dims, n_atts;
  int dims[NC_MAX_VAR_DIMS];
  float d[n*m];
  nc_type type;

  stat = nc_open(file_name, 0, &ncid);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_inq_varid (ncid, var_name, &v_id);
  check_err(stat, __LINE__, __FILE__);

  stat = nc_get_var_float(ncid, v_id, d);
  
  for (i=0; i<n; ++i) {
    for (j=0; j<m; ++j) {
      data[i][j] = d[i*m + j];
      data[i][j] /*mm/yr*/ *= 1. / ( 1000. /*mm/m*/ * secpera /*s/yr*/);
      //printf("data[%d][%d]: %f\n", i, j, data[i][j]);
    }
  }

}

void project_to(Projection *from_p, Projection *to_p, float **dataIn, float **dataOut) {
  float curr_lat, curr_lon;
  LP p;
  XY xy;
  float x,y;
  float f11, f12, f21, f22;
  int x1, x2, y1, y2;
  int i, j;
  int m,n;
  int to_m, to_n; 

  for (i=0; i<to_p->n; ++i) {
    for (j=0; j<to_p->m; ++j) {
      curr_lat = to_p->lat[i][j];
      curr_lon = to_p->lon[i][j];

      if (from_p != NULL) {
        // grab the x-y coordininates here the corresponding 
        // 4 lat-lon points will be around.      
        p.phi = curr_lat * DEG_TO_RAD;
        p.lam = curr_lon * DEG_TO_RAD;
        xy = pj_fwd(p, from_p->ref);
        x = (xy.x - from_p->xll)/GRID_SIZE;
        y = (xy.y - from_p->yll)/GRID_SIZE;
        n = from_p->n;
        m = from_p->m;
      } else {
        y = fmod(curr_lon + 360, 360);
        x = 89.5 - curr_lat;
        n = 180;
        m = 360;
      }
      if (from_p != NULL && (x < 0 || ceil(x) >= n || y < 0 || ceil(y) >= m)) {
        dataOut[i][j] = missing_value;
      } else {
        x1 = floor(x);
        x2 = ceil(x);
        y1 = ceil(y);
        y2 = floor(y);
        
        if (from_p == NULL) {
          x1 = fmod(x1, n);
          x2 = fmod(x2, n);
          y1 = fmod(y1, m);
          y2 = fmod(y2, m);
        }

        f11 = dataIn[x1][y1];
        f21 = dataIn[x2][y1];
        f22 = dataIn[x2][y2];
        f12 = dataIn[x1][y2];

        /* do the bi-linear interpolation */
        if (f11 == missing_value || f21 == missing_value || f12 == missing_value || f22 == missing_value) {
          dataOut[i][j] = missing_value;
        } else {
          if (x2 - x1 == 0 && y2 - y1 == 0) {
            //printf("Special 1\n");
            dataOut[i][j] = f11;
          }
          if (x2 - x1 == 0) {
            //printf("Special 2 (%d, %d)\n", i, j);
            dataOut[i][j] = (y2 - y) * f11 / (y2 - y1) + (y - y1) * f12 / (y2 - y1);
          } else if (y2 - y1 == 0) {
            //printf("Special 3 (%d, %d)\n", i, j);
            dataOut[i][j] = (x2 - x) * f11 / (x2 - x1) + (x - x1) * f21 / (x2 - x1);
          } else {
            dataOut[i][j] = ( ( f11 ) * (x2 - x) * (y2 - y) ) / ( (x2 - x1 ) * (y2 - y1) ) +
                            ( ( f21 ) * (x - x1) * (y2 - y) ) / ( (x2 - x1 ) * (y2 - y1) ) +
                            ( ( f12 ) * (x2 - x) * (y - y1) ) / ( (x2 - x1 ) * (y2 - y1) ) +
                            ( ( f22 ) * (x - x1) * (y - y1) ) / ( (x2 - x1 ) * (y2 - y1) );
          }
        }
      }
    }
  }

}

int main(int argc, char **argv) {
  int i, j;

  Projection PARCA_proj;
  Projection NSIDC_proj;  
  Projection commo_proj;

  /* initialize and get lat and lon for each point in data sets */
  PARCA_proj.n = 290;
  PARCA_proj.m = 541;

  NSIDC_proj.n = 301;
  NSIDC_proj.m = 561;

  commo_proj.n = 332;
  commo_proj.m = 564;

  init_proj(&PARCA_proj);
  init_proj(&NSIDC_proj);
  init_proj(&commo_proj);

  printf("About to unproject PARCA data.\n");
  get_PARCA_projection(&PARCA_proj);
  printf("About to unproject NSIDC data.\n");
  get_NSIDC_projection(&NSIDC_proj);
  printf("About generate common projection.\n");
  get_commo_projection(&commo_proj);

  /* read in the data sets */

  /* first deal with PARCA */
  float **accum_raw;

  accum_raw = (float **)malloc(sizeof(float *) * PARCA_proj.n);
  for (i=0; i<PARCA_proj.n; ++i) {
    accum_raw[i] = (float *)malloc(sizeof(float) * PARCA_proj.m);
  }

  printf("Reading file: \"PARCA/accum_corrected.txt\"\n");
  read_PARCA_data(accum_raw, PARCA_proj.m, PARCA_proj.n, "PARCA/accum_corrected.txt");

  write_PARCA_to_net_cdf(PARCA_proj, accum_raw);

  /* now read in the NSIDC and write it to NetCDF*/
  float **bed_raw, **H_raw;

  bed_raw = (float **)malloc(sizeof(float*) * NSIDC_proj.n);
  H_raw = (float **)malloc(sizeof(float*) * NSIDC_proj.n);
  for (i=0; i<NSIDC_proj.n; ++i) {
    bed_raw[i] = (float *)malloc(sizeof(float) * NSIDC_proj.m);
    H_raw[i] = (float *)malloc(sizeof(float) * NSIDC_proj.m);
  }

  printf("Reading file: \"NSIDC/5k/bed_5km_corrected\"\n");
  read_NSIDC_data(bed_raw, NSIDC_proj.m, NSIDC_proj.n, "NSIDC/5k/bed_5km_corrected");
  printf("Reading file: \"NSIDC/5k/thick_5km_corrected\"\n");
  read_NSIDC_data(H_raw, NSIDC_proj.m, NSIDC_proj.n, "NSIDC/5k/thick_5km_corrected");

  write_NSIDC_to_net_cdf(NSIDC_proj, bed_raw, H_raw);

  /* time for the drad data */
  float **up_raw;

  up_raw = (float **)malloc(sizeof(float *) * 180);
  for (i=0; i<180; ++i) {
      up_raw[i] = (float *)malloc(sizeof(float) * 360);
  }

  printf("Reading file: \"drad/drad.1grid.0.nc\"\n");
  read_drad_data(up_raw, 180, 360, "drad/drad.1grid.0.nc", "Drad_d");


  /* Since we have all the data, it's time to project it all to a common
   * projection, then write it to another nc file
   */
  float **accum_com;

  accum_com = (float **)malloc(sizeof(float *) * commo_proj.n);
  for (i=0; i<commo_proj.n; ++i) {
    accum_com[i] = (float *)malloc(sizeof(float) * commo_proj.m);
  }
  
  float **bed_com, **H_com;

  bed_com = (float **)malloc(sizeof(float*) * commo_proj.n);
  H_com = (float **)malloc(sizeof(float*) * commo_proj.n);
  for (i=0; i<commo_proj.n; ++i) {
    bed_com[i] = (float *)malloc(sizeof(float) * commo_proj.m);
    H_com[i] = (float *)malloc(sizeof(float) * commo_proj.m);
  }

  float **up_com;

  up_com = (float **)malloc(sizeof(float *) * commo_proj.n);
  for (i=0; i<commo_proj.n; ++i) {
    up_com[i] = (float *)malloc(sizeof(float) * commo_proj.m);
  }

  printf("projecting PARCA accum onto common grid\n");
  project_to(&PARCA_proj, &commo_proj, accum_raw, accum_com);
  printf("projecting NSIDC bed onto common grid\n");
  project_to(&NSIDC_proj, &commo_proj, bed_raw, bed_com);
  printf("projecting NSIDC H onto common grid\n");
  project_to(&NSIDC_proj, &commo_proj, H_raw, H_com);
  printf("projecting drad uplift onto common grid\n");
  project_to(NULL, &commo_proj, up_raw, up_com);

  for (i=0; i<commo_proj.n; ++i) {
    for (j=0; j<commo_proj.m; ++j) {
      if (H_com[i][j] == missing_value) {
        H_com[i][j] = 0;
      }
      if (accum_com[i][j] < 0) {
        accum_com[i][j] = missing_value;
      }
      if (bed_com[i][j] == missing_value) {
        bed_com[i][j] = -.1;
      }
    }
  }

  write_commo_to_net_cdf(commo_proj, accum_com, bed_com, H_com, up_com);
}

/*
 * gcc -o commonProjection commonProjection.c -lm -lnetcdf -lproj4
 */
