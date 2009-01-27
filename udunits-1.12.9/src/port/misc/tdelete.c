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

void *tdelete(key, rp, compar)
/* delete node with given key */
const void	*key;			/* key to be deleted */
void **rp;	/* address of the root of tree */
int	(*compar)();		/* comparison function */
{
    node *p;
    register node *q;
    register node *r;
    int cmp;
    register node	**rootp = (node**)rp;	/* address of the root of tree */

    if (rootp == (struct node_t **)0 || (p = *rootp) == (struct node_t *)0)
	return ((struct node_t *)0);
    while ((cmp = (*compar)(key, (*rootp)->key)) != 0)
    {
	p = *rootp;
	rootp = (cmp < 0) ?
	    &(*rootp)->left :		/* follow left branch */
	    &(*rootp)->right;		/* follow right branch */
	if (*rootp == (struct node_t *)0)
	    return ((struct node_t *)0);	/* key not found */
    }
    r = (*rootp)->right;			/* D1: */
    if ((q = (*rootp)->left) == (struct node_t *)0)	/* Left (struct node_t *)0? */
	q = r;
    else if (r != (struct node_t *)0)		/* Right link is null? */
    {
	if (r->left == (struct node_t *)0)	/* D2: Find successor */
	{
	    r->left = q;
	    q = r;
	}
	else
	{			/* D3: Find (struct node_t *)0 link */
	    for (q = r->left; q->left != (struct node_t *)0; q = r->left)
		r = q;
	    r->left = q->right;
	    q->left = (*rootp)->left;
	    q->right = (*rootp)->right;
	}
    }
    free((struct node_t *) *rootp);	/* D4: Free node */
    *rootp = q;				/* link parent to new node */
    return(p);
}
