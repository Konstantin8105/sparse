#include "test_function.h"

// solve a linear system using Cholesky, LU, and QR, with various orderings */
csi demo2(problem* Prob)
{
    cs *A, *C;
    double *b, *x, *resid, t, tol;
    csi k, m, n, ok, order, nb, ns, *r, *s, *rr, sprank;
    csd* D;
    if (!Prob)
        return (0);
    A = Prob->A;
    C = Prob->C;
    b = Prob->b;
    x = Prob->x;
    resid = Prob->resid;
    m = A->m;
    n = A->n;
    tol = Prob->sym ? 0.001 : 1; // partial pivoting tolerance */
    D = cs_dmperm(C, 1); // randomized dmperm analysis */
    if (!D) {
        printf("D is nil\n");
        return (0);
    }

    nb = D->nb;
    r = D->r;
    s = D->s;
    rr = D->rr;
    sprank = rr[3];
    for (ns = 0, k = 0; k < nb; k++) {
        ns += ((r[k + 1] == r[k] + 1) && (s[k + 1] == s[k] + 1));
    }
#ifdef PRINT
    printf("blocks: %g singletons: %g structural rank: %g\n",
        (double)nb, (double)ns, (double)sprank);
#endif // PRINT
    cs_dfree(D);

    for (order = 0; order <= 3; order += 3) // natural and amd(A'*A) */
    {
#ifdef PRINT
        printf("Order : %d\n", order);
        printf("M is : %d\n", m);
        // if (!order && m > 1000) continue ; */
        printf("Start order : %d\n", order);
        printf("QR   ");
        print_order(order);
#endif // PRINT
        rhs(x, b, m); // compute right-hand side */
        t = tic();
        ok = cs_qrsol(order, C, x); // min norm(Ax-b) with QR */
#ifdef PRINT
        printf("time: %8.2f ", toc(t));
        print_resid(ok, C, x, b, resid); // print residual */
        if (ok) {
            for (int r = 0; r < m; r++)
                printf("x[%d] = %10e\n", r, (double)x[r]);
        }
#endif // PRINT
    }

#ifdef PRINT
    printf("m,n,sprank : %d:%d:%d\n", m, n, sprank);
#endif // PRINT

    if (m != n || sprank < n)
        return (1); // return if rect. or singular*/
    for (order = 0; order <= 3; order++) // try all orderings */
    {
#ifdef PRINT
        printf("Order : %d\n", order);
        printf("M is : %d\n", m);
        // if (!order && m > 1000) continue ; */
        printf("Start order : %d\n", order);
        printf("LU   ");
        print_order(order);
#endif // PRINT
        rhs(x, b, m); // compute right-hand side */
        t = tic();
        ok = cs_lusol(order, C, x, tol); // solve Ax=b with LU */
#ifdef PRINT
        printf("time: %8.2f ", toc(t));
        print_resid(ok, C, x, b, resid); // print residual */
        if (ok) {
            for (int r = 0; r < m; r++)
                printf("x[%d] = %10e\n", r, (double)x[r]);
        }
#endif // PRINT
    }

#ifdef PRINT
    printf("Problem sym is : %d\n", (int)Prob->sym);
#endif // PRINT

    if (!Prob->sym)
        return (1);
    for (order = 0; order <= 1; order++) // natural and amd(A+A') */
    {
#ifdef PRINT
        printf("Order : %d\n", order);
        printf("M is : %d\n", m);
        // if (!order && m > 1000) continue ; */
        printf("Start order : %d\n", order);
        printf("Chol ");
        print_order(order);
#endif // PRINT
        rhs(x, b, m); // compute right-hand side */
        t = tic();
        ok = cs_cholsol(order, C, x); // solve Ax=b with Cholesky */
#ifdef PRINT
        printf("time: %8.2f ", toc(t));
        print_resid(ok, C, x, b, resid); // print residual */
        if (ok) {
            for (int r = 0; r < m; r++)
                printf("x[%d] = %10e\n", r, (double)x[r]);
        }
#endif // PRINT
    }
    return (1);
}

// cs_demo2: read a matrix and solve a linear system
int main(void)
{
    problem* Prob = get_problem(stdin, 1e-14);
    // print_problem(Prob) ; */
    int result = demo2(Prob);
	printf("Result demo2 : %d\n", result);
    free_problem(Prob);
    return (0);
}
