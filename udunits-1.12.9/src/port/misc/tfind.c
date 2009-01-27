/*
 * Tree search generalized from Knuth (6.2.2) Algorithm T just like
 * the AT&T man page says.
 *
 * The node_t structure is for internal use only, lint doesn't grok it.
 *
 * Written by reading the System V Interface Definition, not the code.
 *
 * Totally public domain.
 */
/*LINTLIBRARY*/

#include <search.h>
#include "search-node.h"

/* find datum in search tree */
void*
tfind(key, rp, compar)
    const void 	*key;		/* key to be located */
    void	*const*rp;	/* address of tree root */
    int		(*compar)();	/* ordering function */
{
    register node *q;
    register node **rootp = (node**)rp;		/* address of tree root */

    if (rootp == (struct node_t **)0)
	return 0;
    while (*rootp != (struct node_t *)0)	/* Knuth's T1: */
    {
	int r;

	if ((r = (*compar)(key, (*rootp)->key)) == 0)	/* T2: */
	    return &(*rootp)->key;		/* we found it! */
	rootp = (r < 0) ?
	    &(*rootp)->left :			/* T3: follow left branch */
	    &(*rootp)->right;			/* T4: follow right branch */
    }
    return 0;
}
