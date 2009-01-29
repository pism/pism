/*
 * $Id: udalloc.c,v 1.1.1.1 1995/06/15 22:32:00 steve Exp $
 *
 * This file implements the Unidata memory-allocation abstraction.
 */

/*LINTLIBRARY*/

#include "udposix.h"
#include <stddef.h>		/* for size_t and NULL */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "udalloc.h"

#ifdef lint
    static void	lint_malloc(n) size_t n; { ++n; }
    static void	lint_realloc(p,n) voidp p; size_t n; { n+=(size_t)p; }
#   define	malloc(n)	(lint_malloc(n), (voidp)NULL)
#   define	realloc(p,n)	(lint_realloc(p,n), (voidp)NULL)
#endif


/*
 *	Allocate storage.  If unsuccessful, print error message and return
 *	a null pointer.
 */

    voidp
udmalloc(nbytes)
    size_t	nbytes;
{
    if (nbytes > 0) {
	voidp	mem	= malloc(nbytes);

	return mem;

    } else {
	return NULL;
    }
}


/*
 *	Reallocate storage.  If unsuccessful, print error message and return
 *	a null pointer.  A null pointer is also returned if the requested 
 *	size is non-positive.  If a null pointer is returned, then the
 *	contents pointed at by the original pointer should be unchanged.
 */

    voidp
udrealloc(ptr, nbytes)
    voidp	ptr;
    size_t	nbytes;
{
    if (nbytes > 0) {
	voidp	mem	= ptr == NULL 
				?         udmalloc(nbytes) 
				: (voidp) realloc(ptr, nbytes);

	return mem;

    } else {
	return NULL;
    }
}


/*
 * Duplicate a given number of characters.  0-terminate the string.
 */
    char*
udstrndup(from, len)
    const char	*from;		/* string */
    size_t	len;		/* string length (excluding '\0' ) */
{
    char	*to;

    if (from == NULL) {
	to	= NULL;
    } else {
	to	= UD_ALLOC((size_t)(len+1), char);

	if (to != NULL) {
	    (void)strncpy(to, from, len);
	    to[len]	= 0;
	}
    }

    return to;
}


/*
 * Duplicate a string.
 */
    char*
udstrdup(from)
    const char	*from;
{
    return udstrndup(from, strlen(from));
}
