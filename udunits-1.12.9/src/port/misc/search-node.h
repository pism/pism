#ifndef	UD_SEARCH_NODE_H
#define	UD_SEARCH_NODE_H


typedef struct node_t
{
    const void	  *key;
    struct node_t *left, *right;
}
node;


#endif	/* header-file lockout */
