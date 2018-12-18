#include "cs.h"
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
/* (KI) : ignore entry with zero value */
csi cs_entry (cs *T, csi i, csi j, double x)
{
    if (!CS_TRIPLET (T) || i < 0 || j < 0) return (0) ;     /* check inputs */
    if (T->nz >= T->nzmax && !cs_sprealloc (T,2*(T->nzmax))) return (0) ;
	/* change amount of rows and columns, if need */
    T->m = CS_MAX (T->m, i+1) ;
    T->n = CS_MAX (T->n, j+1) ;
	/* do not add zero entry */
    if (T->x && x == 0.0) return (1) ;
	/* add entry at the end */
    if (T->x) T->x [T->nz] = x ;
    T->i [T->nz] = i ;
    T->p [T->nz++] = j ;
    return (1) ;
}
