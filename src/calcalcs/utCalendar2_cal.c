/*
    The CalCalcs routines, a set of C-language routines to perform
    calendar calculations.

    Version 1.01, released 9 January 2010

    Copyright (C) 2010 David W. Pierce, dpierce@ucsd.edu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "udunits2.h"
#include "calcalcs.h"

#include "utCalendar2_cal.h"

static int have_initted=0;
static calcalcs_cal *cal_std=NULL;

/* Following controls the rounding precision of the routines. I.e., if we end up with a value
   such as 59.99999999 seconds, we are going to round it to 60 seconds. The value is given
   in seconds, so, for example, 1.e-3 means round up to 1 second if the value is 0.999 seconds or greater,
   and 1.e-6 means round up to 1 second if the value is 0.999999 seconds or greater.
*/
static double sec_rounding_value = 1.e-8;


/* Internal to this file only */
static void 	initialize( ut_system *units_system );
static void 	get_origin( ut_unit *dataunits, int *y0, int *mon0, int *d0, int *h0, int *min0, double *s0 );
static cv_converter *get_day_to_user_converter( ut_unit *uu, int y0, int mon0, int d0, int h0, int min0, double s0 );
static cv_converter *get_user_to_day_converter( ut_unit *uu, int y0, int mon0, int d0, int h0, int min0, double s0 );
static calcalcs_cal *getcal( const char *name );
static void unknown_cal_emit_warning( const char *calendar_name );

/* Stores previuosly initialized calendars and their names */
static int		maxcals_known=100;
static int		ncals_known=0;
static calcalcs_cal	**known_cal;		/* ptr to array of calcals_cal ptrs */
static char		**known_cal_name;

/* Stores previously emitted "unknown calendar" warnings */
#define UTC2_MAX_UNKCAL_WARNS 	1000
static	char *unknown_cal_emitted_warning_for[ UTC2_MAX_UNKCAL_WARNS ];
static int n_unkcal=0;

/*========================================================================================
 * Turns the passed value into a Y/M/D date 
 */
int utCalendar2_cal( double val, ut_unit *dataunits, int *year, int *month, int *day, int *hour,
                                int *minute, double *second, const char *calendar_name )
{

	int	jdnew, ndinc, ierr, iorig, iround;
	double	fdays, extra_seconds, tot_extra_seconds;
	int	ndays;
	calcalcs_cal	*cal2use;

	/* Following vars are saved between invocations and reused if the
	 * passed units are the same as last time.
	 */
	static	ut_unit *prev_units=NULL;
	static	cv_converter *conv_user_units_to_days=NULL;
	static	int	y0, mon0, d0, h0, min0, jday0;
	static	double  s0, extra_seconds0;
	static	char *prev_calendar=NULL;

	if( dataunits == NULL ) {
		fprintf( stderr, "Error, utCalendar2 passed a NULL units\n" );
		return( UT_ENOINIT );
		}

	if( have_initted == 0 ) 
		initialize( ut_get_system(dataunits) );

	/* Get the calendar we will be using, based on the passed name
	 */
	cal2use = getcal( calendar_name );
	if( cal2use == NULL ) {
		unknown_cal_emit_warning( calendar_name );
		cal2use = getcal( "Standard" );
		}

	/* See if we are being passed the same units and calendar as last time.  If so,
	 * we can optimize by not recomputing all this junk 
	 */
	if( (prev_units != NULL) && (prev_calendar != NULL) 
			&& (strcmp(prev_calendar,cal2use->name)==0) 
			&& (ut_compare( prev_units, dataunits ) == 0)) {
		/* Units and calendar are same as used last call */
		}
	else
		{
		/* Units/calendar combo are different from previously saved units, must redo calcs */

		if( prev_units != NULL ) 
			ut_free( prev_units );

		if( prev_calendar != NULL )
			free( prev_calendar );

		/* Get origin day of the data units */
		get_origin( dataunits, &y0, &mon0, &d0, &h0, &min0, &s0 );	/* Note: static vars */

		/* Number of seconds into the specified origin day */
		extra_seconds0 = h0*3600.0 + min0*60.0 + s0;			/* Note: static vars */

		/* Convert the origin day to Julian Day number in the specified calendar */
		if( (ierr = ccs_date2jday( cal2use, y0, mon0, d0, &jday0 )) != 0 ) {
			fprintf( stderr, "Error in utCalendar2: %s\n", ccs_err_str(ierr) );
			return( UT_EINVALID );
			}

		/* Get converter from user-specified units to "days" */
		if( conv_user_units_to_days != NULL ) 
			cv_free( conv_user_units_to_days );
		conv_user_units_to_days = get_user_to_day_converter( dataunits, y0, mon0, d0, h0, min0, s0 );

		/* Save these units so we can reuse our time-consuming
		 * calculations next time if they are the same units
		 */
		prev_units = ut_clone( dataunits );
		if( ut_compare( prev_units, dataunits ) != 0 ) {
			fprintf( stderr, "error, internal error in udunits2 library found in routine utCalendar2: a clone of the user's units does not equal the original units!\n" );
			exit(-1);
			}

		prev_calendar = (char *)malloc( sizeof(char) * (strlen(cal2use->name)+1 ));
		strcpy( prev_calendar, cal2use->name );
		}

	/* Convert user value of offset to floating point days */
	fdays = cv_convert_double( conv_user_units_to_days, val );

	/* Get integer number of days and seconds past that */
	ndays = fdays;	
	extra_seconds = (fdays - ndays)*86400.0;

	/* Get new Julian day */
	jdnew = jday0 + ndays;

	/* Handle the sub-day part */
	tot_extra_seconds = extra_seconds0 + extra_seconds;
	ndinc = tot_extra_seconds / 86400.0;
	jdnew += ndinc;
	tot_extra_seconds -= ndinc*86400.0;
	if( tot_extra_seconds < 0.0 ) {
		tot_extra_seconds += 86400.0;
		jdnew--;
		}

	/* Convert to a date */
	if( (ierr = ccs_jday2date( cal2use, jdnew, year, month, day )) != 0 ) {
		fprintf( stderr, "Error in utCalendar2: %s\n", ccs_err_str(ierr) );
		return( UT_EINVALID );
		}

	*hour = tot_extra_seconds / 3600.0;
	tot_extra_seconds -= *hour * 3600.0;
	*minute = tot_extra_seconds / 60.0;
	tot_extra_seconds -= *minute * 60.0;
	*second = tot_extra_seconds;

	/* Handle the rouding issues */
	iorig  = *second;			/* Integer conversion */
	iround = *second + sec_rounding_value;	
	if( iround > iorig ) {
		/* printf( "rounding alg invoked, orig date: %04d-%02d-%02d %02d:%02d:%.20lf\n", *year, *month, *day, *hour, *minute, *second ); */
		*second = (double)iround;
		if( *second >= 60.0 ) {
			*second -= 60.0;
			*minute += 1.0;
			if( *minute >= 60.0 ) {
				*minute -= 60.0;
				*hour += 1.0;
				if( *hour >= 24.0 ) {
					*hour -= 24.0;
					if( (ierr = ccs_jday2date( cal2use, jdnew+1, year, month, day )) != 0 ) {
						fprintf( stderr, "Error in utCalendar2: %s\n", ccs_err_str(ierr) );
						return( UT_EINVALID );
						}
					}
				}
			}
		/* printf( "after rounding alg here is new date: %04d-%02d-%02d %02d:%02d:%02.20lf\n", *year, *month, *day, *hour, *minute, *second ); */
		}

	return(0);
}

/*========================================================================================
 * Turn the passed Y/M/D date into a value in the user's units
 */
int utInvCalendar2_cal( int year, int month, int day, int hour, int minute, double second,
	ut_unit *user_unit, double *value, const char *calendar_name )
{
	int	jday, ierr, diff_in_days;
	double	fdiff_in_days, val_days, val_partdays, fdiff_in_partdays, fpartday;
	calcalcs_cal *cal2use;

	/* Following vars are static and retained between invocations for efficiency */
	static	ut_unit *prev_units=NULL;
	static int	y0, mon0, d0, h0, min0, jday0;
	static double	s0, fpartday0;
	static	cv_converter *conv_days_to_user_units=NULL;
	static char 	*prev_calendar=NULL;

	if( have_initted == 0 ) 
		initialize( ut_get_system(user_unit) );

	/* Get the calendar we will be using, based on the passed name
	 */
	cal2use = getcal( calendar_name );
	if( cal2use == NULL ) {
		unknown_cal_emit_warning( calendar_name );
		cal2use = getcal( "Standard" );
		}

	if( (prev_units != NULL) && (prev_calendar != NULL) 
			&& (strcmp(prev_calendar,cal2use->name)==0) 
			&& (ut_compare( prev_units, user_unit ) == 0)) {
		/* Units are same as used last call */
		}
	else
		{
		if( prev_units != NULL )
			ut_free( prev_units );

		if( prev_calendar != NULL )
			free( prev_calendar );

		if( conv_days_to_user_units != NULL ) 
			cv_free( conv_days_to_user_units );

		/* Get origin day of the data units */
		get_origin( user_unit, &y0, &mon0, &d0, &h0, &min0, &s0 );	/* Note: static vars */

		/* Convert the origin day to Julian Day number in the specified calendar */
		if( (ierr = ccs_date2jday( cal2use, y0, mon0, d0, &jday0 )) != 0 ) {
			fprintf( stderr, "Error in utCalendar2: %s\n", ccs_err_str(ierr) );
			return( UT_EINVALID );
			}

		/* Get the origin's HMS in fractional (floating point) part of a Julian day */
		fpartday0 = (double)h0/24.0 + (double)min0/1440.0 + s0/86400.0;

		/* Get converter for turning days into user's units */
		conv_days_to_user_units = get_day_to_user_converter( user_unit, y0, mon0, d0, h0, min0, s0 );

		/* Save these units so we can reuse our time-consuming
		 * calculations next time if they are the same units
		 */
		prev_units = ut_clone( user_unit );
		if( ut_compare( prev_units, user_unit ) != 0 ) {
			fprintf( stderr, "error, internal error in udunits2 library found in routine utInvCalendar2: a clone of the user's units does not equal the original units!\n" );
			exit(-1);
			}

		prev_calendar = (char *)malloc( sizeof(char) * (strlen(cal2use->name)+1 ));
		strcpy( prev_calendar, cal2use->name );
		}

	/* Turn passed date into a Julian day */
	if( (ierr = ccs_date2jday( cal2use, year, month, day, &jday )) != 0 ) {
		fprintf( stderr, "Error in utInvCalendar2: %s\n", ccs_err_str(ierr) );
		return( UT_EINVALID );
		}

	/* jday and jday0 can be very large and nearly equal, so we difference
	 * them first to keep the precision high
	 */
	diff_in_days = jday - jday0;
	fdiff_in_days = (double)diff_in_days;

	/* Get the fractional (floating point) part of a Julian day difference
	 */
	fpartday = (double)hour/24.0 + (double)minute/1440.0 + second/86400.0;
	fdiff_in_partdays = fpartday - fpartday0;

	/* Convert days and partial days to user's units */
	val_days     = cv_convert_double( conv_days_to_user_units, fdiff_in_days     );
	val_partdays = cv_convert_double( conv_days_to_user_units, fdiff_in_partdays );

	/* Hopefully this will minimize the roundoff errors */
	*value = val_days + val_partdays;

	return(0);
}

/*==============================================================================================
 * Get a converter that turns the user's units into days 
 */
static cv_converter *get_user_to_day_converter( ut_unit *uu, int y0, int mon0, int d0, int h0, int min0, double s0 ) 
{
	char		daystr[1024];
	ut_unit 	*udu_days;
	cv_converter	*conv_user_units_to_days;

	sprintf( daystr, "days since %04d-%02d-%02d %02d:%02d:%f",
		y0, mon0, d0, h0, min0, s0 );
		
	udu_days = ut_parse( ut_get_system(uu), daystr, UT_ASCII );
	if( udu_days == NULL ) {
		fprintf( stderr, "internal error in utCalendar2/conv_to_days: failed to parse following string: \"%s\"\n",
			daystr );
		exit(-1);
		}
	conv_user_units_to_days = ut_get_converter( uu, udu_days );
	if( conv_user_units_to_days == NULL ) {
		fprintf( stderr, "internal error in utCalendar2/conv_to_days: cannot convert from \"%s\" to user units\n",
		 	daystr );
		exit(-1);
		}

	ut_free( udu_days );
	return( conv_user_units_to_days );
}

/*==============================================================================================
 * Get a converter that turns days into the user's units 
 */
static cv_converter *get_day_to_user_converter( ut_unit *uu, int y0, int mon0, int d0, int h0, int min0, double s0 ) 
{
	char		daystr[1024];
	ut_unit 	*udu_days;
	cv_converter	*conv_days_to_user_units;

	sprintf( daystr, "days since %04d-%02d-%02d %02d:%02d:%f",
		y0, mon0, d0, h0, min0, s0 );
		
	udu_days = ut_parse( ut_get_system(uu), daystr, UT_ASCII );
	if( udu_days == NULL ) {
		fprintf( stderr, "internal error in utCalendar2/conv_to_user_units: failed to parse following string: \"%s\"\n",
			daystr );
		exit(-1);
		}
	conv_days_to_user_units = ut_get_converter( udu_days, uu );
	if( conv_days_to_user_units == NULL ) {
		fprintf( stderr, "internal error in utCalendar2/conv_to_user_units: cannot convert from user units to \"%s\"\n",
		 	daystr );
		exit(-1);
		}

	free( udu_days );
	return( conv_days_to_user_units );
}

/*==========================================================================================
 * The user specified some origin to the time units. For example, if the units string
 * were "days since 2005-10-15", then the origin date is 2005-10-15.  This routine
 * deduces the specified origin date from the passed units structure 
 */
static void get_origin( ut_unit *dataunits, int *y0, int *mon0, int *d0, int *h0, int *min0, double *s0 )
{
	double		s0lib, rez, tval, tval_conv, resolution;
	cv_converter	*conv_user_date_to_ref_date;
	ut_unit		*tmpu;
	int		y0lib, mon0lib, d0lib, h0lib, min0lib;
	char		ustr[1024];

	static ut_unit 	*udu_ref_date=NULL;	/* Saved between invocations */

	if( udu_ref_date == NULL ) {

		/* Make a timestamp units that refers to the udunits2 library's intrinsic
		 * time origin.  The first thing we do is parse a timestampe unit
		 * (any old timestamp unit) and immediately discard the result, to account
		 * for a bug in the udunits-2 library that fails to set the library's
		 * time origin unless this step is done first.  Without this, a call to
		 * ut_decode_time with tval==0 returns year=-4713 rather than 2001.
		 */
		if( (tmpu = ut_parse( ut_get_system(dataunits), "days since 2010-01-09", UT_ASCII )) == NULL ) {
			fprintf( stderr, "Internal error in routnie utCalendar2/get_origin: failed to parse temp timestamp string\n" );
			exit(-1);
			}
		ut_free( tmpu );

		tval = 0.0;
		ut_decode_time( tval, &y0lib, &mon0lib, &d0lib, &h0lib, &min0lib, &s0lib, &rez );
		sprintf( ustr, "seconds since %04d-%02d-%02d %02d:%02d:%f",
			y0lib, mon0lib, d0lib, h0lib, min0lib, s0lib );
		udu_ref_date = ut_parse( ut_get_system(dataunits), ustr, UT_ASCII );
		if( udu_ref_date == NULL ) {	
			fprintf( stderr, "internal error in routine utCalendar2/get_origin: could not parse origin string \"%s\"\n",
				ustr );
			exit(-1);
			}
		}

	/* Get converter from passed units to the library's intrinsic time units */
	conv_user_date_to_ref_date = ut_get_converter( dataunits, udu_ref_date );

	/* convert tval=0 to ref date */
	tval = 0.0;
	tval_conv = cv_convert_double( conv_user_date_to_ref_date, tval );

	/* Now decode the converted value */
	ut_decode_time( tval_conv, y0, mon0, d0, h0, min0, s0, &resolution );

	cv_free( conv_user_date_to_ref_date );
}

/*========================================================================================*/
static void initialize( ut_system *units_system )
{
	int	i;

	/* Make space for our saved calendars */
	known_cal = (calcalcs_cal **)malloc( sizeof( calcalcs_cal * ) * maxcals_known );
	if( known_cal == NULL ) {
		fprintf( stderr, "Error in utCalendar2 routines, could not allocate internal storage\n" );
		exit(-1);
		}
	for( i=0; i<maxcals_known; i++ )
		known_cal[i] = NULL;
	known_cal_name = (char **)malloc( sizeof( char * ) * maxcals_known );
	if( known_cal_name == NULL ) {
		fprintf( stderr, "Error in utCalendar2 routines, could not allocate internal storage for calendar names\n" );
		exit(-1);
		}
	for( i=0; i<maxcals_known; i++ )
		known_cal_name[i] = NULL;

	/* Make a standard calendar */
	cal_std = ccs_init_calendar( "Standard" );

	have_initted = 1;	/* a global */
}

/*========================================================================================
 * Returns NULL if the passed calendar name is both not found and not creatable
 */
static calcalcs_cal *getcal( const char *name )
{
	int	i, new_index;
	calcalcs_cal *new_cal;

	if( cal_std == NULL ) {
		fprintf( stderr, "Coding error in utCalendar2_cal routines, cal_std is null!\n" );
		exit(-1);
		}

	if( strcasecmp( name, "standard" ) == 0 )
		return( cal_std );

	/* See if it is one of the previously known calendars */
	for( i=0; i<ncals_known; i++ ) {
		if( strcmp( known_cal_name[i], name ) == 0 ) 
			return( known_cal[i] );
		}

	/* If we get here, the cal is not known, so create it if possible */
	new_cal = ccs_init_calendar( name );
	if( new_cal == NULL ) 
		return( NULL );		/* unknown name */

	/* We now know a new calendar */
	new_index = ncals_known;
	ncals_known++;

	/* new_index is where we will be storing the new calendar */
	if( ncals_known > maxcals_known ) {
		ncals_known = maxcals_known;
		new_index = strlen(name);	/* some arbitrary value */
		if( new_index >= maxcals_known )
			new_index = 10;
		}

	/* If there was already a calendar stored in this slot 
	 * (because we might be reusing slots) then clear it out
	 */
	if( known_cal[new_index] != NULL ) 
		ccs_free_calendar( known_cal[new_index] );
	
	if( known_cal_name[new_index] != NULL ) 
		free( known_cal_name[new_index] );

	known_cal[new_index] = new_cal;
	known_cal_name[new_index] = (char *)malloc( sizeof(char) * (strlen(name)+1 ));
	strcpy( known_cal_name[new_index], name );

	return( new_cal );
}

/*=============================================================================================
 */
static void unknown_cal_emit_warning( const char *calendar_name )
{
	int	i;

	for( i=0; i<n_unkcal; i++ ) {
		if( strcmp( calendar_name, unknown_cal_emitted_warning_for[i] ) == 0 ) 
			/* Already emitted a warning for this calendar */
			return;
		}

	if( n_unkcal == UTC2_MAX_UNKCAL_WARNS )
		/* Too many warnings already given, give up */
		return;

	/* OK, emit the warning now */
	fprintf( stderr, "Error in utCalendar2_cal/utInvCalendar2_cal: unknown calendar \"%s\".  Using Standard calendar instead\n",
		calendar_name );

	/* Save the fact that we have warned for this calendar */
	unknown_cal_emitted_warning_for[ n_unkcal ] = (char *)malloc( sizeof(char) * (strlen(calendar_name)+1 ));
	if( unknown_cal_emitted_warning_for[ n_unkcal ] == NULL ) {
		fprintf( stderr, "Error in utCalendar_cal: could not allocate internal memory to store string \"%s\"\n",
			calendar_name );
		return;
		}

	strcpy( unknown_cal_emitted_warning_for[ n_unkcal ], calendar_name );
	n_unkcal++;
}
