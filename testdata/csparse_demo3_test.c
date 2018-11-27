#include "csparse_demo_test.h"

// free workspace for demo3 */
static csi done3 (csi ok, css *S, csn *N, double *y, cs *W, cs *E, csi *p)
{
    cs_sfree (S) ;
    cs_nfree (N) ;
    cs_free (y) ;
    cs_spfree (W) ;
    cs_spfree (E) ;
    cs_free (p) ;
    return (ok) ;
}

// Cholesky update/downdate */
csi demo3 (problem *Prob)
{
    cs *A, *C, *W = NULL, *WW, *WT, *E = NULL, *W2 ;
    csi n, k, *Li, *Lp, *Wi, *Wp, p1, p2, *p = NULL, ok ;
    double *b, *x, *resid, *y = NULL, *Lx, *Wx, s,  t, t1 ;
    css *S = NULL ;
    csn *N = NULL ;
    if (!Prob || !Prob->sym || Prob->A->n == 0) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; resid = Prob->resid;
    n = A->n ;
    if (!Prob->sym || n == 0) return (1) ;
    rhs (x, b, n) ;                             // compute right-hand side */
    printf ("\nchol then update/downdate ") ;
    print_order (1) ;
    y = cs_malloc (n, sizeof (double)) ;
    t = tic () ;
    S = cs_schol (1, C) ;                       // symbolic Chol, amd(A+A') */
    printf ("\nsymbolic chol time %8.2f\n", toc (t)) ;
    t = tic () ;
    N = cs_chol (C, S) ;                        // numeric Cholesky */
    printf ("numeric  chol time %8.2f\n", toc (t)) ;
    if (!S || !N || !y) return (done3 (0, S, N, y, W, E, p)) ;
    t = tic () ;
    cs_ipvec (S->pinv, b, y, n) ;               // y = P*b */
    cs_lsolve (N->L, y) ;                       // y = L\y */
    cs_ltsolve (N->L, y) ;                      // y = L'\y */
    cs_pvec (S->pinv, y, x, n) ;                // x = P'*y */
    printf ("solve    chol time %8.2f\n", toc (t)) ;
    printf ("original: ") ;
    print_resid (1, C, x, b, resid) ;           // print residual */
    k = n/2 ;                                   // construct W  */
    W = cs_spalloc (n, 1, n, 1, 0) ;
    if (!W) return (done3 (0, S, N, y, W, E, p)) ;
    Lp = N->L->p ; Li = N->L->i ; Lx = N->L->x ;
    Wp = W->p ; Wi = W->i ; Wx = W->x ;
    Wp [0] = 0 ;
    p1 = Lp [k] ;
    Wp [1] = Lp [k+1] - p1 ;
    s = Lx [p1] ;
    // srand (1) ;
	double counter = 0.001;
    for ( ; p1 < Lp [k+1] ; p1++)
    {
        p2 = p1 - Lp [k] ;
        Wi [p2] = Li [p1] ;
        Wx [p2] = s * counter; // rand () / ((double) RAND_MAX) ;
		counter *= 1.05;
    }
    t = tic () ;
    ok = cs_updown (N->L, +1, W, S->parent) ;   // update: L*L'+W*W' */
    t1 = toc (t) ;
    printf ("update:   time: %8.2f\n", t1) ;
    if (!ok) return (done3 (0, S, N, y, W, E, p)) ;
    t = tic () ;
    cs_ipvec (S->pinv, b, y, n) ;               // y = P*b */
    cs_lsolve (N->L, y) ;                       // y = L\y */
    cs_ltsolve (N->L, y) ;                      // y = L'\y */
    cs_pvec (S->pinv, y, x, n) ;                // x = P'*y */
    t = toc (t) ;
    p = cs_pinv (S->pinv, n) ;
    W2 = cs_permute (W, p, NULL, 1) ;           // E = C + (P'W)*(P'W)' */
	WT = cs_transpose (W2,1) ;
    WW = cs_multiply (W2, WT) ;
    cs_spfree (WT) ;
    cs_spfree (W2) ;

    E = cs_add (C, WW, 1, 1) ;
    cs_spfree (WW) ;
    if (!E || !p) return (done3 (0, S, N, y, W, E, p)) ;
    printf ("update:   time: %8.2f (incl solve) ", t1+t) ;
    print_resid (1, E, x, b, resid) ;           // print residual */
    cs_nfree (N) ;                              // clear N */
    t = tic () ;
    N = cs_chol (E, S) ;                        // numeric Cholesky */

	test_cs_print(E,0);
    if (!N) return (done3 (0, S, N, y, W, E, p)) ;
    cs_ipvec (S->pinv, b, y, n) ;               // y = P*b */
    cs_lsolve (N->L, y) ;                       // y = L\y */
    cs_ltsolve (N->L, y) ;                      // y = L'\y */
    cs_pvec (S->pinv, y, x, n) ;                // x = P'*y */
    t = toc (t) ;
    printf ("rechol:   time: %8.2f (incl solve) ", t) ;
    print_resid (1, E, x, b, resid) ;           // print residual */
    t = tic () ;
    ok = cs_updown (N->L, -1, W, S->parent) ;   // downdate: L*L'-W*W' */
    t1 = toc (t) ;
    if (!ok) return (done3 (0, S, N, y, W, E, p)) ;
    printf ("downdate: time: %8.2f\n", t1) ;
    t = tic () ;
    cs_ipvec (S->pinv, b, y, n) ;               // y = P*b */
    cs_lsolve (N->L, y) ;                       // y = L\y */
    cs_ltsolve (N->L, y) ;                      // y = L'\y */
    cs_pvec (S->pinv, y, x, n) ;                // x = P'*y */
    t = toc (t) ;
    printf ("downdate: time: %8.2f (incl solve) ", t1+t) ;
    print_resid (1, C, x, b, resid) ;           // print residual */
    return (done3 (1, S, N, y, W, E, p)) ;
} 


// cs_demo3: read a matrix and test Cholesky update/downdate */
int main (void)
{
    problem *Prob = get_problem (stdin, 1e-14) ;
    int result = demo3 (Prob) ;
	printf("Result demo3 : %d\n", result);
    free_problem (Prob) ;
    return (0) ;
}
