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

static void trecurse(root, action, level)
/* Walk the nodes of a tree */
register const node	*root;		/* Root of the tree to be walked */
register void	(*action)();	/* Function to be called at each node */
register int	level;
{
    if (root->left == (struct node_t *)0 && root->right == (struct node_t *)0)
	(*action)(&root->key, leaf, level);
    else
    {
	(*action)(&root->key, preorder, level);
	if (root->left != (struct node_t *)0)
	    trecurse(root->left, action, level + 1);
	(*action)(&root->key, postorder, level);
	if (root->right != (struct node_t *)0)
	    trecurse(root->right, action, level + 1);
	(*action)(&root->key, endorder, level);
    }
}

void twalk(rt, action)		/* Walk the nodes of a tree */
const void	*rt;		/* Root of the tree to be walked */
void	(*action)();		/* Function to be called at each node */
{
    const node	*root = (const node*)rt;/* Root of the tree to be walked */

    if (root != (node *)0 && action != (void(*)())0)
	trecurse(root, action, 0);
}
