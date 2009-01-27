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

void *tsearch(key, rp, compar)
/* find or insert datum into search tree */
const void 	*key;		/* key to be located */
void	**rp;			/* address of tree root */
int	(*compar)();		/* ordering function */
{
    register node *q;
    register node	**rootp = (node**)rp;	/* address of tree root */

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
    q = (node *) malloc(sizeof(node));		/* T5: key not found */
    if (q != (struct node_t *)0)		/* make new node */
    {
	*rootp = q;				/* link new node to old */
	q->key = key;				/* initialize new node */
	q->left = q->right = (struct node_t *)0;
    }
    return &q->key;
}
