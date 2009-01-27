#include <stdio.h>	/* for diagnostics */
#include <stdlib.h>
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "udunits.h"

extern utUnit	*utNew();

static int
not_here(s)
char *s;
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(name, arg)
char *name;
int arg;
{
    errno = 0;
    switch (*name) {
    case 'A':
	break;
    case 'B':
	break;
    case 'C':
	break;
    case 'D':
	break;
    case 'E':
	break;
    case 'F':
	break;
    case 'G':
	break;
    case 'H':
	break;
    case 'I':
	break;
    case 'J':
	break;
    case 'K':
	break;
    case 'L':
	break;
    case 'M':
	break;
    case 'N':
	break;
    case 'O':
	break;
    case 'P':
	break;
    case 'Q':
	break;
    case 'R':
	break;
    case 'S':
	break;
    case 'T':
	break;
    case 'U':
	if (strEQ(name, "EALLOC"))
#ifdef UT_EALLOC
	    return UT_EALLOC;
#else
	    goto not_there;
#endif
	if (strEQ(name, "ECONVERT"))
#ifdef UT_ECONVERT
	    return UT_ECONVERT;
#else
	    goto not_there;
#endif
	if (strEQ(name, "EINVALID"))
#ifdef UT_EINVALID
	    return UT_EINVALID;
#else
	    goto not_there;
#endif
	if (strEQ(name, "EIO"))
#ifdef UT_EIO
	    return UT_EIO;
#else
	    goto not_there;
#endif
	if (strEQ(name, "ENOFILE"))
#ifdef UT_ENOFILE
	    return UT_ENOFILE;
#else
	    goto not_there;
#endif
	if (strEQ(name, "ENOINIT"))
#ifdef UT_ENOINIT
	    return UT_ENOINIT;
#else
	    goto not_there;
#endif
	if (strEQ(name, "ENOROOM"))
#ifdef UT_ENOROOM
	    return UT_ENOROOM;
#else
	    goto not_there;
#endif
	if (strEQ(name, "ENOTTIME"))
#ifdef UT_ENOTTIME
	    return UT_ENOTTIME;
#else
	    goto not_there;
#endif
	if (strEQ(name, "EOF"))
#ifdef UT_EOF
	    return UT_EOF;
#else
	    goto not_there;
#endif
	if (strEQ(name, "ESYNTAX"))
#ifdef UT_ESYNTAX
	    return UT_ESYNTAX;
#else
	    goto not_there;
#endif
	if (strEQ(name, "EUNKNOWN"))
#ifdef UT_EUNKNOWN
	    return UT_EUNKNOWN;
#else
	    goto not_there;
#endif
	if (strEQ(name, "MAXNUM_BASE_QUANTITIES"))
#ifdef UT_MAXNUM_BASE_QUANTITIES
	    return UT_MAXNUM_BASE_QUANTITIES;
#else
	    goto not_there;
#endif
	if (strEQ(name, "NAMELEN"))
#ifdef UT_NAMELEN
	    return UT_NAMELEN;
#else
	    goto not_there;
#endif
	break;
    case 'V':
	break;
    case 'W':
	break;
    case 'X':
	break;
    case 'Y':
	break;
    case 'Z':
	break;
    }
    errno = EINVAL;
    return 0;

not_there:
    errno = ENOENT;
    return 0;
}


MODULE = UDUNITS	PACKAGE = UDUNITS


double
constant(name,arg)
	char *		name
	int		arg


int
init(path)
    char *	path
    CODE:
    {
	RETVAL = utInit(path);
    }
    OUTPUT:
	RETVAL


void
term()
    CODE:
    {
	utTerm();
    }


utUnit *
new()
    CODE:
    {
	RETVAL = utNew();

	if (RETVAL == 0)
	    croak("Couldn't allocate %lu bytes for new unit structure",
		  sizeof(utUnit));
    }
    OUTPUT:
	RETVAL


void
scan(spec)
    char *	spec
    CODE:
    {
	int	status;
	utUnit	*unit = utNew();

	if (unit == 0)
	    croak("Couldn't allocate %lu bytes for new unit structure",
		  sizeof(utUnit));

	ST(0) = sv_newmortal();

	status = utScan(spec, unit);

	if (status == UT_ENOINIT)
	    croak("units module not initialized");
	else if (status == 0)
	    sv_setref_pv(ST(0), "utUnitPtr", (void*)unit);
    }


MODULE = UDUNITS	PACKAGE = utUnitPtr


int
istime(unit)
    utUnit *	unit
    CODE:
    {
	RETVAL = utIsTime(unit);
    }
    OUTPUT:
	RETVAL


int
hasorigin(unit)
    utUnit *	unit
    CODE:
    {
	RETVAL = utHasOrigin(unit);
    }
    OUTPUT:
	RETVAL


void
clear(unit)
    utUnit *	unit
    CODE:
    {
	(void) utClear(unit);
    }


utUnit *
dup(source)
    utUnit *	source
    CODE:
    {
	utUnit	*dest = utNew();

	if (dest == 0)
	    croak("Couldn't allocate %lu bytes for new unit structure",
		  sizeof(utUnit));

	RETVAL = utCopy(source, dest);
    }
    OUTPUT:
	RETVAL


void
shift(unit, amount)
    utUnit *	unit
    double	amount
    CODE:
    {
	(void) utShift(unit, amount, unit);
    }


void
scale(unit, coefficient)
    utUnit *	unit
    double	coefficient
    CODE:
    {
	(void) utScale(unit, coefficient, unit);
    }


void
multiply(unit, otherunit)
    utUnit *	unit
    utUnit *	otherunit
    CODE:
    {
	(void) utMultiply(unit, otherunit, unit);
    }


void
invert(unit)
    utUnit *	unit
    CODE:
    {
	(void) utInvert(unit, unit);
    }


void
divide(unit, divisor)
    utUnit *	unit
    utUnit *	divisor
    CODE:
    {
	(void) utDivide(unit, divisor, unit);
    }


void
raise(unit, power)
    utUnit *	unit
    int		power
    CODE:
    {
	(void) utRaise(unit, power, unit);
    }


void
print(unit)
    utUnit *	unit
    CODE:
    {
	char	*staticbuf;

	(void) utPrint(unit, &staticbuf);
	ST(0) = sv_newmortal();
	sv_setpv((SV*)ST(0), staticbuf);
    }


int
convert(from_unit, to_unit, slope, intercept)
    utUnit *	from_unit
    utUnit *	to_unit
    double	&slope = NO_INIT
    double	&intercept = NO_INIT
    CODE:
    {
	RETVAL = utConvert(from_unit, to_unit, &slope, &intercept);
    }
    OUTPUT:
	slope
	intercept
	RETVAL


int
valtocal(unit, value, year, month, day, hour, minute, second)
    utUnit *	unit
    double	value
    int		&year = NO_INIT
    int		&month = NO_INIT
    int		&day = NO_INIT
    int		&hour = NO_INIT
    int		&minute = NO_INIT
    float	&second = NO_INIT
    CODE:
    {
	RETVAL = utCalendar(value, unit,
			    &year, &month, &day,
			    &hour, &minute, &second);
    }
    OUTPUT:
	year
	month
	day
	hour
	minute
	second
	RETVAL


double
caltoval(unit, year, month, day, hour, minute, second)
    utUnit *	unit
    int		year
    int		month
    int		day
    int		hour
    int		minute
    double	second
    CODE:
    {
	int	status = utInvCalendar(year, month, day,
				       hour, minute, second,
				       unit, &RETVAL);

	if (status == UT_EINVALID)
	    croak("not a unit of time");
	if (status == UT_ENOINIT)
	    croak("units module not initialized");
    }
    OUTPUT:
	RETVAL


void
DESTROY(unit)
    utUnit *	unit
    CODE:
    {
	free((char*)unit);
    }
