/*
 * $Id: atexit.c,v 1.1.1.1 1995/06/15 22:31:57 steve Exp $
 */

/*LINTLIBRARY*/


#include "udposix.h"
#include "stdlib.h"


    int
atexit(func)
    void	(*func)UD_PROTO((void));
{
    return on_exit(func, (char*)0);
}
