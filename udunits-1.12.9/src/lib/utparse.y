%{
/*
 * $Id: utparse.y,v 1.2 2007/11/13 15:39:43 steve Exp $
 *
 * yacc(1)-based parser for decoding formatted units(3) specifications.
 */

/*LINTLIBRARY*/

#include <stdio.h>
#include <string.h>
#include "udunits.h"
#include "utscan.h"
#include "utprivate.h"

/*
 * Allocate storage for the package global variables, which are declared in
 * header file "utprivate.h".
 */
int	UtLineno;		/* input-file line index */
int	UnitNotFound;		/* parser didn't find unit */
utUnit	*FinalUnit;		/* fully-parsed specification */

%}

%union {
    double	rval;			/* floating-point numerical value */
    long	ival;			/* integer numerical value */
    char	name[UT_NAMELEN];	/* name of a quantity */
    utUnit	unit;			/* "unit" structure */
}

%token  <ival>	INT
%token  <ival>	ERR
%token	<ival>	SHIFT
%token  <ival>	SPACE
%token  <ival>	MULTIPLY
%token  <ival>	DIVIDE
%token  <ival>	EXPONENT
%token  <rval>	REAL
%token  <name>	NAME
%token	<rval>	DATE
%token	<rval>	TIME
%token	<rval>	ZONE

%type   <rval>	number_exp
%type   <rval>	value_exp
%type   <rval>	timestamp
%type   <rval>	time_exp
%type   <unit>	unit_exp
%type   <unit>	power_exp
%type   <unit>	factor_exp
%type   <unit>	quant_exp
%type   <unit>	origin_exp
%type	<unit>	unit_spec

%%

unit_spec:        /* nothing */			{
			YYACCEPT;
		  }
		| origin_exp			{
			(void)utCopy(&$1, FinalUnit);
			YYACCEPT;
		  }
		| error				{
			yyerrok;
			yyclearin;
			YYABORT;
		  }
		;

origin_exp:	  unit_exp			{
			(void)utCopy(&$1, &$$);
		  }
		| unit_exp SHIFT value_exp	{
			if (utIsTime(&$1)) {
			    /*
			     * The shift amount is divided by the unit scale
			     * factor in the following because the shift amount
			     * must be in the units of the first argument (e.g.
			     * 0.555556 kelvins for the fahrenheit unit) and a
			     * timestamp isn't necessarily so proportioned.
			    (void)utShift(&$1, $3, &$$);
			     */
			    (void)utShift(&$1, $3/$1.factor, &$$);
			} else {
			    (void) utShift(&$1, $3, &$$);
			}
		  }
		| unit_exp SHIFT timestamp	{
			if (utIsTime(&$1)) {
			    (void)utShift(&$1, $3/$1.factor, &$$);
			} else {
			    UnitNotFound	= 1;
			    YYERROR;
			}
		  }
		    ;

unit_exp:	  power_exp			{
			(void)utCopy(&$1, &$$);
		  }
		| unit_exp power_exp	{
			(void)utMultiply(&$1, &$2, &$$);
		  }
		| unit_exp MULTIPLY power_exp	{
			(void)utMultiply(&$1, &$3, &$$);
		  }
		| unit_exp DIVIDE power_exp	{
			(void)utDivide(&$1, &$3, &$$);
		  }
		;

power_exp:	  factor_exp			{
			(void)utCopy(&$1, &$$);
		  }
		| power_exp INT			{
			(void)utRaise(&$1, $2, &$$);
		  }
		| power_exp EXPONENT INT	{
			(void)utRaise(&$1, $3, &$$);
		  }
		;

factor_exp:	  number_exp			{
			utUnit	unit;
			(void)utScale(utClear(&unit), $1, &$$);
		  }
		| quant_exp			{
			(void)utCopy(&$1, &$$);
		  }
		| '(' origin_exp ')'		{
			(void)utCopy(&$2, &$$);
		  }
		;

quant_exp:	  NAME                          {
			utUnit     unit;
			if (utFind($1, &unit) != 0) {
			    UnitNotFound	= 1;
			    YYERROR;
			}
			(void)utCopy(&unit, &$$);
		  }
		| NAME INT			{
			utUnit     unit;
			if (utFind($1, &unit) != 0) {
			    UnitNotFound	= 1;
			    YYERROR;
			}
			(void)utRaise(&unit, $2, &$$);
		  }
		;

value_exp:	  number_exp			{ $$ = $1; }
		| '(' value_exp ')'		{ $$ = $2; }
		;

number_exp:	  INT				{ $$ = $1; }
		| REAL				{ $$ = $1; }
		;

timestamp:	  time_exp			{ $$ = $1; }
		| '(' timestamp ')'		{ $$ = $2; }
		;

time_exp:	  DATE 				{ $$ = $1; }
		| DATE TIME 			{ $$ = $1 + $2; }
		| DATE TIME ZONE 		{ $$ = $1 + ($2 - $3); }
		;

%%
