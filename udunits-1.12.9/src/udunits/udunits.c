/*
 * $Id: udunits.c,v 1.3 2008/04/04 17:29:08 steve Exp $
 *
 * Test program for the Unidata udunits(3) library.
 */

#ifndef _XOPEN_SOURCE
#  define _XOPEN_SOURCE 500
#endif
#include <limits.h>	/* for _POSIX_MAX_INPUT */
#include <stddef.h>	/* for size_t */
#include <stdio.h>	/* for I/O functions, and NULL */
#include <string.h>	/* for strlen() & memmove() */
#include <errno.h>	/* for errno */
#include <ctype.h>	/* for isspace() */
#include "udunits.h"

#ifdef YYDEBUG
    extern int	utdebug;
#endif

static char	*progname;

#define	ABS(x)		((x) < 0 ? -(x) : (x))


    static void
WriteError(name)
    char	*name;
{
    (void) fprintf(stderr, "%s: Error writing to standard output: ", 
		   name);
    perror("");
}


    static void
ReadError(name)
    char	*name;
{
    (void) fprintf(stderr, "%s: Error reading from standard output: ", 
		   name);
    perror("");
}


/*
 * Get a unit specification.
 */
   static int
GetSpec(prompt, spec, size)
    char	*prompt;
    char	*spec;
    size_t	size;
{
    int		status;
    static char	me[]	= "GetSpec";

    if (fputs(prompt, stdout) == EOF) {
	WriteError(me);
	status	= -1;
    } else if (fgets(spec, size, stdin) == NULL) {
	putchar('\n');
	if (feof(stdin)) {
	    status	=  0;
	} else {
	    ReadError(progname);
	    status	= -1;
	}
    } else {
	/*
	 * Trim whitespace from the specification.
	 */
	char	*start, *stop;

	for (start = spec; *start != 0; ++start)
	    if (!isspace(*start))
		break;
	for (stop = start + strlen(start); stop > start; --stop)
	    if (!isspace(stop[-1]))
		break;
	*stop	= 0;
	(void)memmove((void*)spec, (void*)start, (size_t)(stop-start+1));

	status	= 1;
    }

    return status;
}


main(argc, argv)
    int			argc;
    char		**argv;
{
    int			status;
    int			c;
    char		HaveSpec[_POSIX_MAX_INPUT+1];
    char		WantSpec[_POSIX_MAX_INPUT+1];
    utUnit		HaveUnit, WantUnit;
    static char		*UnitsPath	= NULL;
    extern int		optind;
    extern char		*optarg;

    progname	= argv[0];

    while ((c = getopt(argc, argv, "d")) != -1)
	switch (c) {
	case 'd':
#ifdef YYDEBUG
	    utdebug	= 1;
#endif
	    break;
	default:
	    (void) fprintf(stderr, "%s: Unknown option `%c'\n", progname, c);
	    exit(1);
	}

    if (optind < argc)
	UnitsPath	= argv[optind];

    if (utInit(UnitsPath) != 0) {
	(void) fprintf(stderr, 
		       "%s: Couldn't initialize udunits(3) package\n", 
		       progname);

    } else {
	for (;;) {
	    int		i;
	    double	slope, intercept;

	    i	=  GetSpec("You have: ", HaveSpec, sizeof(HaveSpec));
	    if (i == -1)
		goto failure;
	    if (i == 0)
		goto success;
	    if (HaveSpec[0] == 0)
		continue;
	    if (utScan(HaveSpec, &HaveUnit) != 0) {
		(void) fprintf(stderr, 
			       "%s: Don't recognize `%s'\n", 
			       progname, HaveSpec);
		if (ferror(stderr))
		    goto failure;
	    } else {
		for (;;) {
		    i	= GetSpec("You want: ", WantSpec, 
				      sizeof(WantSpec));
		    if (i == -1)
			goto failure;
		    if (i == 0)
			goto success;
		    if (WantSpec[0] == 0) {
			char	*s;

			(void)utPrint(&HaveUnit, &s);
			if (printf("    Definition: \"%s\"\n", s) < 0) {
			    WriteError(progname);
			    goto failure;
			}
			if (utIsTime(&HaveUnit) && atoi(HaveSpec) != 0) {
			    int	year, month, day, hour, minute;
			    float	second;

			    (void) utCalendar(1.0, &HaveUnit, 
					      &year, &month, &day,
					      &hour, &minute, &second);
			    (void) printf(
			    "    \"%s\" is %d-%02d-%02d %02d:%02d:%g UTC\n",
					  HaveSpec, year, month, day, hour,
					  minute, second);
			}
		    } else if (utScan(WantSpec, &WantUnit) != 0) {
			(void) fprintf(stderr, 
				       "%s: Don't recognize `%s'\n", 
				       progname, WantSpec);
			if (ferror(stderr))
			    goto failure;
			continue;
		    } else {
			i	= utConvert(&HaveUnit, &WantUnit, 
					    &slope, &intercept);
			if (i == UT_ECONVERT) {
			    (void) fprintf(stderr,
					   "%s: Units are incompatible\n", 
					   progname);
			    if (ferror(stderr))
				goto failure;
			} else if (i == UT_EINVALID) {
			    (void) fprintf(stderr, 
					   "%s: A unit is corrupted\n", 
					   progname);
			    if (ferror(stderr))
				goto failure;
			} else {
			    if (intercept == 0) {
				if (printf("    <%s> = <%s>*%g\n", 
					   WantSpec, HaveSpec, slope)
				       < 0)  {
				    WriteError(progname);
				    goto failure;
				}
				if (printf("    <%s> = <%s>/%g\n", 
					   WantSpec, HaveSpec, 1./slope)
				       < 0)  {
				    WriteError(progname);
				    goto failure;
				}
			    } else {
				if (printf("    <%s> = <%s>*%g %s %g\n", 
					   WantSpec, HaveSpec, slope,
					   intercept < 0
						? "-"
						: "+",
					   ABS(intercept)) < 0) {
				    WriteError(progname);
				    goto failure;
				}
				if (printf("    <%s> = <%s>/%g %s %g\n", 
					   WantSpec, HaveSpec, 1./slope,
					   intercept < 0
						? "-"
						: "+",
					   ABS(intercept)) < 0) {
				    WriteError(progname);
				    goto failure;
				}
			    }
			}		/* successful conversion */
		    }			/* "want" decoded */
		    break;
		}			/* "want" loop */
	    }				/* "have" decoded */
	}				/* "have" loop */
    }					/* units package initialized */

    failure:
	status	= 1;
	goto exit;
    success:
	status	= 0;
	goto exit;

    exit:
	utTerm();
	return status;
}
