/*
 * strftime.c
 *
 * Public-domain relatively quick-and-dirty implemenation of
 * ANSI library routine for System V Unix systems.
 *
 * It's written in old-style C for maximal portability.
 * However, since I'm used to prototypes, I've included them too.
 *
 * If you want stuff in the System V ascftime routine, add the SYSV_EXT define.
 * For stuff needed to implement the P1003.2 date command, add POSIX2_DATE.
 * For complete POSIX semantics, add POSIX_SEMANTICS.
 *
 * The code for %c, %x, and %X is my best guess as to what's "appropriate".
 * This version ignores LOCALE information.
 * It also doesn't worry about multi-byte characters.
 * So there.
 *
 * Arnold Robbins
 * January, February, March, 1991
 * Updated March 1992
 *
 * Fixes from ado@elsie.nci.nih.gov
 * February 1991
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>

#ifndef __STDC__
#define const	/**/

extern void *malloc();
extern void *realloc();
extern void tzset();
extern char *strchr();
extern char *getenv();
static int weeknumber();
#else
extern void *malloc(unsigned count);
extern void *realloc(void *ptr, unsigned count);
extern void tzset(void);
extern char *strchr(const char *str, int ch);
extern char *getenv(const char *v);
static int weeknumber(const struct tm *timeptr, int firstweekday);
#endif

#ifdef __GNUC__
#define inline	__inline__
#else
#define inline	/**/
#endif

#define range(low, item, hi)	max(low, min(item, hi))

extern char *tzname[2];
extern int daylight;

#define SYSV_EXT	1	/* stuff in System V ascftime routine */
#define POSIX2_DATE	1	/* stuff in Posix 1003.2 date command */
#define VMS_EXT		1	/* include %V for VMS date format */
#define POSIX_SEMANTICS 1	/* call tzset() if TZ changes */

#if defined(POSIX2_DATE) && ! defined(SYSV_EXT)
#define SYSV_EXT	1
#endif

/* min --- return minimum of two numbers */

#ifndef __STDC__
static inline int
min(a, b)
int a, b;
#else
static inline int
min(int a, int b)
#endif
{
	return (a < b ? a : b);
}

/* max --- return maximum of two numbers */

#ifndef __STDC__
static inline int
max(a, b)
int a, b;
#else
static inline int
max(int a, int b)
#endif
{
	return (a > b ? a : b);
}

/* strftime --- produce formatted time */

#ifndef __STDC__
size_t
strftime(s, maxsize, format, timeptr)
char *s;
size_t maxsize;
const char *format;
const struct tm *timeptr;
#else
size_t
strftime(char *s, size_t maxsize, const char *format, const struct tm *timeptr)
#endif
{
	char *endp = s + maxsize;
	char *start = s;
	char tbuf[100];
	int i;
	static short first = 1;
	static char *savetz = NULL, *oldtz = NULL;
	static int savetzlen = 0;
	char *tz;
	int tzlen;

	/* various tables, useful in North America */
	static char *days_a[] = {
		"Sun", "Mon", "Tue", "Wed",
		"Thu", "Fri", "Sat",
	};
	static char *days_l[] = {
		"Sunday", "Monday", "Tuesday", "Wednesday",
		"Thursday", "Friday", "Saturday",
	};
	static char *months_a[] = {
		"Jan", "Feb", "Mar", "Apr", "May", "Jun",
		"Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
	};
	static char *months_l[] = {
		"January", "February", "March", "April",
		"May", "June", "July", "August", "September",
		"October", "November", "December",
	};
	static char *ampm[] = { "AM", "PM", };

	if (s == NULL || format == NULL || timeptr == NULL || maxsize == 0)
		return 0;

	if (strchr(format, '%') == NULL && strlen(format) + 1 >= maxsize)
		return 0;

#ifndef POSIX_SEMANTICS
	if (first) {
		tzset();
		first = 0;
	}
#else	/* POSIX_SEMANTICS */
	tz = getenv("TZ");
	tzlen = strlen(tz);
	if (first) {
		if (tz != NULL) {
			savetz = (char *) malloc(tzlen + 1);
			if (savetz != NULL) {
				savetzlen = tzlen + 1;
				strcpy(savetz, tz);
			}
		}
		tzset();
		first = 0;
	}
	/* if we have a saved TZ, and it is different, recapture and reset */
	if (tz && savetz && (tz[0] != savetz[0] || strcmp(tz, savetz) != 0)) {
		i = strlen(tz) + 1;
		if (i > savetzlen) {
			savetz = (char *) realloc(savetz, i);
			if (savetz) {
				savetzlen = i;
				strcpy(savetz, tz);
			}
		} else
			strcpy(savetz, tz);
		tzset();
	}
#endif /* POSIX_SEMANTICS */

	for (; *format && s < endp - 1; format++) {
		tbuf[0] = '\0';
		if (*format != '%') {
			*s++ = *format;
			continue;
		}
	again:
		switch (*++format) {
		case '\0':
			*s++ = '%';
			goto out;

		case '%':
			*s++ = '%';
			continue;

		case 'a':	/* abbreviated weekday name */
			if (timeptr->tm_wday < 0 || timeptr->tm_wday > 6)
				strcpy(tbuf, "?");
			else
				strcpy(tbuf, days_a[timeptr->tm_wday]);
			break;

		case 'A':	/* full weekday name */
			if (timeptr->tm_wday < 0 || timeptr->tm_wday > 6)
				strcpy(tbuf, "?");
			else
				strcpy(tbuf, days_l[timeptr->tm_wday]);
			break;

#ifdef SYSV_EXT
		case 'h':	/* abbreviated month name */
#endif
		case 'b':	/* abbreviated month name */
			if (timeptr->tm_mon < 0 || timeptr->tm_mon > 11)
				strcpy(tbuf, "?");
			else
				strcpy(tbuf, months_a[timeptr->tm_mon]);
			break;

		case 'B':	/* full month name */
			if (timeptr->tm_mon < 0 || timeptr->tm_mon > 11)
				strcpy(tbuf, "?");
			else
				strcpy(tbuf, months_l[timeptr->tm_mon]);
			break;

		case 'c':	/* appropriate date and time representation */
			sprintf(tbuf, "%s %s %2d %02d:%02d:%02d %d",
				days_a[range(0, timeptr->tm_wday, 6)],
				months_a[range(0, timeptr->tm_mon, 11)],
				range(1, timeptr->tm_mday, 31),
				range(0, timeptr->tm_hour, 23),
				range(0, timeptr->tm_min, 59),
				range(0, timeptr->tm_sec, 61),
				timeptr->tm_year + 1900);
			break;

		case 'd':	/* day of the month, 01 - 31 */
			i = range(1, timeptr->tm_mday, 31);
			sprintf(tbuf, "%02d", i);
			break;

		case 'H':	/* hour, 24-hour clock, 00 - 23 */
			i = range(0, timeptr->tm_hour, 23);
			sprintf(tbuf, "%02d", i);
			break;

		case 'I':	/* hour, 12-hour clock, 01 - 12 */
			i = range(0, timeptr->tm_hour, 23);
			if (i == 0)
				i = 12;
			else if (i > 12)
				i -= 12;
			sprintf(tbuf, "%02d", i);
			break;

		case 'j':	/* day of the year, 001 - 366 */
			sprintf(tbuf, "%03d", timeptr->tm_yday + 1);
			break;

		case 'm':	/* month, 01 - 12 */
			i = range(0, timeptr->tm_mon, 11);
			sprintf(tbuf, "%02d", i + 1);
			break;

		case 'M':	/* minute, 00 - 59 */
			i = range(0, timeptr->tm_min, 59);
			sprintf(tbuf, "%02d", i);
			break;

		case 'p':	/* am or pm based on 12-hour clock */
			i = range(0, timeptr->tm_hour, 23);
			if (i < 12)
				strcpy(tbuf, ampm[0]);
			else
				strcpy(tbuf, ampm[1]);
			break;

		case 'S':	/* second, 00 - 61 */
			i = range(0, timeptr->tm_sec, 61);
			sprintf(tbuf, "%02d", i);
			break;

		case 'U':	/* week of year, Sunday is first day of week */
			sprintf(tbuf, "%d", weeknumber(timeptr, 0));
			break;

		case 'w':	/* weekday, Sunday == 0, 0 - 6 */
			i = range(0, timeptr->tm_wday, 6);
			sprintf(tbuf, "%d", i);
			break;

		case 'W':	/* week of year, Monday is first day of week */
			sprintf(tbuf, "%d", weeknumber(timeptr, 1));
			break;

		case 'x':	/* appropriate date representation */
			sprintf(tbuf, "%s %s %2d %d",
				days_a[range(0, timeptr->tm_wday, 6)],
				months_a[range(0, timeptr->tm_mon, 11)],
				range(1, timeptr->tm_mday, 31),
				timeptr->tm_year + 1900);
			break;

		case 'X':	/* appropriate time representation */
			sprintf(tbuf, "%02d:%02d:%02d",
				range(0, timeptr->tm_hour, 23),
				range(0, timeptr->tm_min, 59),
				range(0, timeptr->tm_sec, 61));
			break;

		case 'y':	/* year without a century, 00 - 99 */
			i = timeptr->tm_year % 100;
			sprintf(tbuf, "%d", i);
			break;

		case 'Y':	/* year with century */
			sprintf(tbuf, "%d", 1900 + timeptr->tm_year);
			break;

		case 'Z':	/* time zone name or abbrevation */
			i = 0;
			if (daylight && timeptr->tm_isdst)
				i = 1;
			strcpy(tbuf, tzname[i]);
			break;

#ifdef SYSV_EXT
		case 'n':	/* same as \n */
			tbuf[0] = '\n';
			tbuf[1] = '\0';
			break;

		case 't':	/* same as \t */
			tbuf[0] = '\t';
			tbuf[1] = '\0';
			break;

		case 'D':	/* date as %m/%d/%y */
			strftime(tbuf, sizeof tbuf, "%m/%d/%y", timeptr);
			break;

		case 'e':	/* day of month, blank padded */
			sprintf(tbuf, "%2d", range(1, timeptr->tm_mday, 31));
			break;

		case 'r':	/* time as %I:%M:%S %p */
			strftime(tbuf, sizeof tbuf, "%I:%M:%S %p", timeptr);
			break;

		case 'R':	/* time as %H:%M */
			strftime(tbuf, sizeof tbuf, "%H:%M", timeptr);
			break;

		case 'T':	/* time as %H:%M:%S */
			strftime(tbuf, sizeof tbuf, "%H:%M:%S", timeptr);
			break;
#endif


#ifdef VMS_EXT
		case 'V':	/* date as dd-bbb-YYYY */
			sprintf(tbuf, "%2d-%3.3s-%4d",
				range(1, timeptr->tm_mday, 31),
				months_a[range(0, timeptr->tm_mon, 11)],
				timeptr->tm_year + 1900);
			for (i = 3; i < 6; i++)
				if (islower(tbuf[i]))
					tbuf[i] = toupper(tbuf[i]);
			break;
#endif


#ifdef POSIX2_DATE
		case 'C':
			sprintf(tbuf, "%02d", (timeptr->tm_year + 1900) / 100);
			break;


		case 'E':
		case 'O':
			/* POSIX locale extensions, ignored for now */
			goto again;
#endif
		default:
			tbuf[0] = '%';
			tbuf[1] = *format;
			tbuf[2] = '\0';
			break;
		}
		i = strlen(tbuf);
		if (i)
			if (s + i < endp - 1) {
				strcpy(s, tbuf);
				s += i;
			} else
				return 0;
	}
out:
	if (s < endp && *format == '\0') {
		*s = '\0';
		return (s - start);
	} else
		return 0;
}

/* weeknumber --- figure how many weeks into the year */

/* With thanks and tip of the hatlo to ado@elsie.nci.nih.gov */

#ifndef __STDC__
static int
weeknumber(timeptr, firstweekday)
const struct tm *timeptr;
int firstweekday;
#else
static int
weeknumber(const struct tm *timeptr, int firstweekday)
#endif
{
	if (firstweekday == 0)
		return (timeptr->tm_yday + 7 - timeptr->tm_wday) / 7;
	else
		return (timeptr->tm_yday + 7 -
			(timeptr->tm_wday ? (timeptr->tm_wday - 1) : 6)) / 7;
}

#if 0
/* ADR --- I'm loathe to mess with ado's code ... */

Date:         Wed, 24 Apr 91 20:54:08 MDT
From: Michal Jaegermann <audfax!emory!vm.ucs.UAlberta.CA!NTOMCZAK>
To: arnold@audiofax.com

Hi Arnold,
in a process of fixing of strftime() in libraries on Atari ST I grabbed
some pieces of code from your own strftime.  When doing that it came
to mind that your weeknumber() function compiles a little bit nicer
in the following form:
/*
 * firstweekday is 0 if starting in Sunday, non-zero if in Monday
 */
{
    return (timeptr->tm_yday - timeptr->tm_wday +
	    (firstweekday ? (timeptr->tm_wday ? 8 : 1) : 7)) / 7;
}
How nicer it depends on a compiler, of course, but always a tiny bit.

   Cheers,
   Michal
   ntomczak@vm.ucs.ualberta.ca
#endif
