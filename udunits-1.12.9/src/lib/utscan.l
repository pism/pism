%{

/*
 * lex(1) specification for tokens for the Unidata units package, udunits.
 */

#ifndef _XOPEN_SOURCE
#   define _XOPEN_SOURCE
#endif
#ifndef _ANSI_C_SOURCE
#   define _ANSI_C_SOURCE
#endif

#include <udposix.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include "udunits.h"
#include "utprivate.h"	/* private definitions */
#include "utparse.h"	/* #include "udunits.h" before this one */
#include "utscan.h"

%}

space		([ \t])
integer		([+-]?[0-9]+)
exp		([eE][+-]?[0-9]+)
real		([+-]?[0-9]*("."[0-9]*{exp}?|{exp}))
year		([+-]?[0-9]{1,4})
month		(0?[1-9]|1[0-2])
day		(0?[1-9]|[1-2][0-9]|30|31)
hour		([0-1]?[0-9]|2[0-3])
minute		([0-5]?[0-9])
second		(({minute}|60)("."[0-9]*)?)
name		(%|[a-zA-Z][a-zA-Z_]*([0-9]+[a-zA-Z_]+)*)

%Start		shift_seen date_seen time_seen

%%

"#".* {
    /* comment */ ;
}

{space}*("@"|after|AFTER|from|FROM|since|SINCE|ref|REF){space}* {
    BEGIN shift_seen;
    return SHIFT;
}

{space}*(per|PER|"/"){space}* {
    BEGIN 0;
    return DIVIDE;
}

"-"|"."|"*"|{space}+ {
    BEGIN 0;
    return MULTIPLY;
}

"^"|"**" {
    BEGIN 0;
    return EXPONENT;
}

<shift_seen>{year}"-"{month}"-"{day}T?{space}* {
    int		year;
    int		month;
    int		day;

    (void) sscanf((char*)yytext, "%d-%d-%d", &year, &month, &day);
    yylval.rval	= utencdate(year, month, day);

    BEGIN date_seen;
    return DATE;
}

<shift_seen>now{space}* {
    time_t	now;
    struct tm	*tp;

    errno	= 0;
    if (time(&now) == -1) {
	(void) fprintf(stderr, "udunits(3): Can't get current time. %s\n",
		       strerror(errno));
	return ERR;
    }

    if ((tp = gmtime(&now)) == NULL) {
	(void) fprintf(stderr, "udunits(3): can't get UTC time\n");
	return ERR;
    }
    yylval.rval	= utencdate(1900+tp->tm_year, 1+tp->tm_mon, tp->tm_mday) +
		  utencclock(tp->tm_hour, tp->tm_min, (double)tp->tm_sec);

    BEGIN 0;
    return REAL;
}

<date_seen>{hour}(":"{minute}(":"{second})?)?{space}* {
    int		hour;
    int		minute	= 0;
    double	second	= 0.0;

    (void) sscanf((char*)yytext, "%d:%d:%lf", &hour, &minute, &second);
    yylval.rval	= utencclock(hour, minute, second);

    BEGIN time_seen;
    return TIME;
}

<date_seen>[0-2][0-9][0-5][0-9]{space}* {
    int		hour;
    int		minute;

    (void) sscanf((char*)yytext, "%2d%2d", &hour, &minute);
    if (hour >= 24 || minute >= 60)
	return ERR;
    yylval.rval	= utencclock(hour, minute, 0.0);

    BEGIN time_seen;
    return TIME;
}

<date_seen>([0-9][0-5][0-9]){space}* {
    int		hour;
    int		minute;

    (void) sscanf((char*)yytext, "%1d%2d", &hour, &minute);
    if (minute >= 60)
	return ERR;
    yylval.rval	= utencclock(hour, minute, 0.0);

    BEGIN time_seen;
    return TIME;
}

<time_seen>UTC{space}* {
    yylval.rval	= 0;

    BEGIN 0;
    return ZONE;
}

<time_seen>([+-]?{hour}":"{minute}){space}* {
    int		hour;
    int		minute;

    (void) sscanf((char*)yytext, "%d:%d", &hour, &minute);
    if (hour < 0)
	minute	= -minute;

    yylval.rval	= utencclock(hour, minute, 0.0);

    BEGIN 0;
    return ZONE;
}

<time_seen>([+-]?{hour}":"?){space}* {
    int		hour;

    (void) sscanf((char*)yytext, "%d", &hour);

    yylval.rval	= utencclock(hour, 0, 0.0);

    BEGIN 0;
    return ZONE;
}

<time_seen>([+-]{hour}{minute}){space}* {
    int		zone;
    int		hour;
    int		minute;

    (void) sscanf((char*)yytext, "%d", &zone);
    minute	= zone < 0
		    ? -(-zone % 100)
		    :    zone % 100;
    hour	= (zone - minute) / 100;

    yylval.rval	= utencclock(hour, minute, 0.0);

    BEGIN 0;
    return ZONE;
}

{name} {
    if (strlen((char*)yytext) > UT_NAMELEN-1)
	(void) fprintf(stderr, "udunits(3): Name `%s' too long; truncated\n",
		       yytext);
    (void) strncpy(yylval.name, (char*)yytext, UT_NAMELEN-1);
    yylval.name[UT_NAMELEN-1]	= 0;

    BEGIN 0;
    return NAME;
}

{integer} {
    int		status;

    errno	= 0;
    yylval.ival = atol((char*)yytext);
    if (errno == 0) {
	status	= INT;
    } else {
	(void) fputs("udunits(3): Invalid integer\n", stderr);
	status	= ERR;
    }

    BEGIN 0;
    return status;
}

{real} {
    int		status;

    errno	= 0;
    yylval.rval = atof((char*)yytext);
    if (errno == 0) {
	status	= REAL;
    } else {
	(void) fputs("udunits(3): Invalid real\n", stderr);
	status	= ERR;
    }

    BEGIN 0;
    return status;
}

")" {
    BEGIN 0;
    return yytext[0];
}

. {
    return yytext[0];
}

%%

    void
utLexReset()
{
    BEGIN 0;
}
