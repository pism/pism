/*
 * $Id: udalloc.h,v 1.1.1.1 1995/06/15 22:32:00 steve Exp $
 */

#ifndef UD_ALLOC_H_INCLUDED
#define UD_ALLOC_H_INCLUDED

#include <stddef.h>	/* for `size_t' */
#include <stdlib.h>	/* for `*alloc()' */

/*
 *	Interface to the Unidata memory-allocation abstraction:
 */
UD_EXTERN_FUNC(voidp	udmalloc,	(size_t nbytes));
UD_EXTERN_FUNC(voidp	udrealloc,	(voidp ptr, size_t nbytes));
UD_EXTERN_FUNC(char	*udstrdup,	(const char *s));
UD_EXTERN_FUNC(char	*udstrndup,	(const char *s, size_t nbytes));


/*
 *	Some macros to make life easier:
 */
#define UD_ALLOC(theNum, theType) \
	    (theType*)udmalloc((size_t)(sizeof(theType)*(theNum))) 

#define UD_REALLOC(ptr, theNum, theType) \
	    (theType*)udrealloc((voidp)(ptr), \
		    (size_t)(sizeof(theType)*(theNum)))

#define FREE(ptr)		(void)free((voidp)(ptr))

#endif	/* !UD_ALLOC_H_INCLUDED */
