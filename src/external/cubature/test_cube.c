
/*
   Usage: ./test_cube <dim> <tol> <integrand> <maxeval>

   where <dim> = # dimensions, <tol> = relative tolerance,
   <integrand> is either 0/1/2 for the three test integrands (see below),
   and <maxeval> is the maximum # function evaluations (0 for none).
*/

/* author Steven G. Johnson */

#include "cubature.h"
#include <math.h>
#include <gsl/gsl_math.h>       /* M_PI */

int count = 0;
int which_integrand = 0;
const double radius = 0.50124145262344534123412; /* random */

double f_test(unsigned dim, const double *x, void *data)
{
     double val;
     unsigned i;
     ++count;
     switch (which_integrand) {
     case 0: /* discontinuous objective: volume of hypersphere */
          val = 0;
          for (i = 0; i < dim; ++i)
           val += x[i] * x[i];
          val = val < radius * radius;
          break;
     case 1: /* simple smooth (separable) objective: prod. cos(x[i]). */
          val = 1;
          for (i = 0; i < dim; ++i)
           val *= cos(x[i]);
          break;
     case 2: { /* integral of exp(-x^2), rescaled to (0,infinity) limits */
          double scale = 1.0;
          val = 0;
          for (i = 0; i < dim; ++i) {
           double z = (1 - x[i]) / x[i];
           val += z * z;
           scale *= M_2_SQRTPI / (x[i] * x[i]);
          }
          val = exp(-val) * scale;
          break;
     }

     default:
          fprintf(stderr, "unknown integrand %d\n", which_integrand);
          exit(EXIT_FAILURE);
     }
     /* if (count < 100) printf("%d: f(%g, ...) = %g\n", count, x[0], val); */
     return val;
}

/* surface area of n-dimensional unit hypersphere */
static double S(unsigned n)
{
     double val;
     int fact = 1;
     if (n % 2 == 0) { /* n even */
      val = 2 * pow(M_PI, n * 0.5);
      n = n / 2;
      while (n > 1) fact *= (n -= 1);
      val /= fact;
     }
     else { /* n odd */
      val = (1 << (n/2 + 1)) * pow(M_PI, n/2);
      while (n > 2) fact *= (n -= 2);
      val /= fact;
     }
     return val;
}

static double exact_integral(unsigned dim, const double *xmax) {
     unsigned i;
     double val;
     switch(which_integrand) {
     case 0:
          val = dim == 0 ? 1 : S(dim) * pow(radius * 0.5, dim) / dim;
          break;
     case 1:
          val = 1;
          for (i = 0; i < dim; ++i)
           val *= sin(xmax[i]);
          break;
     case 2:
          val = 1;
          break;
     default:
          fprintf(stderr, "unknown integrand %d\n", which_integrand);
              exit(EXIT_FAILURE);
     }
     return val;
}

int main(int argc, char **argv)
{
     double *xmin, *xmax;
     double tol, val, err;
     unsigned i, dim, maxEval;

     dim = argc > 1 ? atoi(argv[1]) : 2;
     tol = argc > 2 ? atof(argv[2]) : 1e-2;
     which_integrand = argc > 3 ? atoi(argv[3]) : 1;
     maxEval = argc > 4 ? atoi(argv[4]) : 0;

     xmin = (double *) malloc(dim * sizeof(double));
     xmax = (double *) malloc(dim * sizeof(double));
     for (i = 0; i < dim; ++i) {
      xmin[i] = 0;
      xmax[i] = 1 + (which_integrand == 2 ? 0 : 0.4 * sin(i));
     }

     printf("%u-dim integral, tolerance = %g, integrand = %d\n",
        dim, tol, which_integrand);
     adapt_integrate(f_test, 0, dim, xmin, xmax, maxEval, 0, tol, &val, &err);
     printf("integration val = %g, est. err = %g, true err = %g\n",
        val, err, fabs(val - exact_integral(dim, xmax)));
     printf("#evals = %d\n", count);

     free(xmax);
     free(xmin);

     return 0;
}
