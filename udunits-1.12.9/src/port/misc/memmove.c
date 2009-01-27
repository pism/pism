/*
 * $Id: memmove.c,v 1.1.1.1 1995/06/15 22:31:57 steve Exp $
 */

/*LINTLIBRARY*/


#include "udposix.h"
#include "string.h"


/*
 * Copy bytes.  Handle overlap correctly.
 */
    voidp
memmove(s1, s2, n)
    voidp	s1;
    const voidp	s2;
    size_t	n;
{
    char	*sc1	= (char*)s1;
    const char	*sc2	= (const char*)s2;

    if (sc2 < sc1 && sc1 < sc2 + n)
	for (sc1 += n, sc2 += n; 0 < n; --n)
	    *--sc1 = *--sc2;		/* copy backwards */
    else
	for (; 0 < n; --n)
	    *sc1++	= *sc2++;	/* copy forwards */

    return s1;
}
