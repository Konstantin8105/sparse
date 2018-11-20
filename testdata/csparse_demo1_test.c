#include "cs.h"

csi test_cs_print (const cs *A, csi brief)
{
    csi p, j, m, n, nzmax, nz, *Ap, *Ai ;
    double *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    printf ("CSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        printf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %10e\n", (double) m,
            (double) n, (double) nzmax, (double) (Ap [n]), cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    col %g : locations %g to %g\n", (double) j, 
                (double) (Ap [j]), (double) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                printf ("      %g : %10e\n", (double) (Ai [p]), Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        printf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double) m,
            (double) n, (double) nzmax, (double) nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            printf ("    %g %g : %10e\n", (double) (Ai [p]), (double) (Ap [p]),
                Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}

int main(void)
{
	// load data
	cs *T = cs_load(stdin);

	// print data
	// cs_print(T,0);

	// A compressed-column form of T
	cs *A = cs_compress (T);

	// print data
	// cs_print(A,0);

	// Transpose
	cs * AT = cs_transpose (A,1);

	// print data
	// cs_print(AT,0);

	// m = # of rows of A
    int m = A ? A->m : 0 ;

	// triplet identify
    T = cs_spalloc (m, m, m, 1, 1) ;
    for (int i = 0 ; i < m ; i++) cs_entry (T, i, i, 1) ;
    cs *Eye = cs_compress (T) ;
	// cs_print(Eye,0);

	// C = A*A'
	cs * C = cs_multiply(A,AT);
	// cs_print(C,0);

	// D = C + Eye*norm(C,1)
	cs * D = cs_add(C,Eye,1,cs_norm(C));
	test_cs_print(D,0);

	return 0;
}
