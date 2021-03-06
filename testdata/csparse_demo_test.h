#include "cs.h"
#include <time.h>

csi test_cs_print(const cs* A, csi brief)
{
#ifdef PRINT
    csi p, j, m, n, nzmax, nz, *Ap, *Ai;
    double* Ax;
    if (!A) {
        printf("(null)\n");
        return (0);
    }
    m = A->m;
    n = A->n;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;
    nzmax = A->nzmax;
    nz = A->nz;
    printf("Sparse\n");
    if (nz < 0) {
        printf("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %10e\n", (double)m,
            (double)n, (double)nzmax, (double)(Ap[n]), cs_norm(A));
        for (j = 0; j < n; j++) {
            printf("    col %g : locations %g to %g\n", (double)j,
                (double)(Ap[j]), (double)(Ap[j + 1] - 1));
            for (p = Ap[j]; p < Ap[j + 1]; p++) {
                printf("      %g : %10e\n", (double)(Ai[p]), Ax ? Ax[p] : 1);
                if (brief && p > 20) {
                    printf("  ...\n");
                    return (1);
                }
            }
        }
    } else {
        printf("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (double)m,
            (double)n, (double)nzmax, (double)nz);
        for (p = 0; p < nz; p++) {
            printf("    %g %g : %10e\n", (double)(Ai[p]), (double)(Ap[p]),
                Ax ? Ax[p] : 1);
            if (brief && p > 20) {
                printf("  ...\n");
                return (1);
            }
        }
    }
#endif // PRINT
    return (1);
}

typedef struct problem_struct {
    cs* A;
    cs* C;
    csi sym;
    double* x;
    double* b;
    double* resid;
} problem;

void print_problem(problem* P)
{
#ifdef PRINT
    printf("Matrix A:\n");
    if (P->A) {
        test_cs_print(P->A, 0);
    };
    printf("Matrix C:\n");
    if (P->C) {
        test_cs_print(P->C, 0);
    };
    printf("sym = %d\n", (int)P->sym);

    printf("Vector x\n");
    for (int i = 0; i < P->A->n; i++)
        printf("x[%d] = %f\n", i, (double)P->x[i]);
    for (int i = 0; i < P->A->n; i++)
        printf("b[%d] = %f\n", i, (double)P->b[i]);
    for (int i = 0; i < P->A->n; i++)
        printf("resid[%d] = %f\n", i, (double)P->resid[i]);
#endif // PRINT
}

// infinity-norm of x */
static double norm(double* x, csi n)
{
    csi i;
    double normx = 0;
    for (i = 0; i < n; i++)
        normx = CS_MAX(normx, fabs(x[i]));
    return (normx);
}

// create a right-hand side */
static void rhs(double* x, double* b, csi m)
{
    csi i;
    for (i = 0; i < m; i++)
        b[i] = 1 + ((double)i) / m;
    for (i = 0; i < m; i++)
        x[i] = b[i];
}

// 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static csi is_sym(cs* A)
{
    csi is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i;
    if (m != n)
        return (0);
    is_upper = 1;
    is_lower = 1;
    for (j = 0; j < n; j++) {
        for (p = Ap[j]; p < Ap[j + 1]; p++) {
            if (Ai[p] > j)
                is_upper = 0;
            if (Ai[p] < j)
                is_lower = 0;
        }
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0));
}

// true for off-diagonal entries */
static csi dropdiag(csi i, csi j, double aij, void* other) { return (i != j); }

// C = A + triu(A,1)' */
static cs* make_sym(cs* A)
{
    cs *AT, *C;
    AT = cs_transpose(A, 1); // AT = A' */
    cs_fkeep(AT, &dropdiag, NULL); // drop diagonal entries from AT */
    C = cs_add(A, AT, 1, 1); // C = A+AT */
    cs_spfree(AT);
    return (C);
}

// free a problem */
problem* free_problem(problem* Prob)
{
    if (!Prob)
        return (NULL);
    cs_spfree(Prob->A);
    if (Prob->sym)
        cs_spfree(Prob->C);
    cs_free(Prob->b);
    cs_free(Prob->x);
    cs_free(Prob->resid);
    return (cs_free(Prob));
}

// compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)) */
static void print_resid(csi ok, cs* A, double* x, double* b, double* resid)
{
    csi i, m, n;
    if (!ok) {
        printf("    (failed)\n");
        return;
    }
    m = A->m;
    n = A->n;
    for (i = 0; i < m; i++)
        resid[i] = -b[i]; // resid = -b */
    cs_gaxpy(A, x, resid); // resid = resid + A*x  */
    // printf ("resid: %8.2e\n", norm (resid,m) / ((n == 0) ? 1 : */
    //     (cs_norm (A) * norm (x,n) + norm (b,m)))) ; */
    printf("\n");
}

static void print_order(csi order)
{
#ifdef PRINT
    switch (order) {
    case 0:
        printf("natural    ");
        break;
    case 1:
        printf("amd(A+A')  ");
        break;
    case 2:
        printf("amd(S'*S)  ");
        break;
    case 3:
        printf("amd(A'*A)  ");
        break;
    }
#endif // PRINT
}

// read a problem from a file; use %g for integers to avoid csi conflicts */
problem* get_problem(FILE* f, double tol)
{
    cs *T, *A, *C;
    csi sym, m, n, mn, nz1, nz2;
    problem* Prob;
    Prob = cs_calloc(1, sizeof(problem));
    if (!Prob)
        return (NULL);
    T = cs_load(f); // load triplet matrix T from a file */
    Prob->A = A = cs_compress(T); // A = compressed-column form of T */
    cs_spfree(T); // clear T */
    if (!cs_dupl(A))
        return (free_problem(Prob)); // sum up duplicates */
    Prob->sym = sym = is_sym(A); // determine if A is symmetric */
    m = A->m;
    n = A->n;
    mn = CS_MAX(m, n);
    nz1 = A->p[n];
    cs_dropzeros(A); // drop zero entries */
	nz2 = A->p[n];

#ifdef PRINT
	printf("n   = %d\n",(int) n  );
	printf("nz1 = %d\n",(int) nz1);
	printf("nz2 = %d\n",(int) nz2);
	printf("A->p[%d] = %d\n", (int)n,(int)A->p[(int)n]);

	printf("A before drop\n");
	test_cs_print(A,0);

	printf("tol = %.5e\n", tol);
#endif // PRINT
    if (tol > 0){
        int ok = cs_droptol(A, tol); // drop tiny entries (just to test) */
#ifdef PRINT
		printf("droptol = %d\n", ok);
#endif // PRINT
	}

#ifdef PRINT
	printf("A before make_sym\n");
	test_cs_print(A,0);
#endif // PRINT


    Prob->C = C = sym ? make_sym(A) : A; // C = A + triu(A,1)', or C=A */
    if (!C)
        return (free_problem(Prob));
#ifdef PRINT
    printf("\n--- Matrix: %g-by-%g, nnz: %g (sym: %g: nnz %g), norm: %8.2e\n",
        (double)m, (double)n, (double)(A->p[n]), (double)sym,
        (double)(sym ? C->p[n] : 0), 
		(double) 0);
		// cs_norm(C));
    if (nz1 != nz2)
        printf("zero entries dropped: %g\n", (double)(nz1 - nz2));
	
	test_cs_print(A,0);

	printf("nz2 = %d\n",(int) nz2);
	printf("A->p[%d] = %d\n", (int) n,(int) A->p[(int) n]);
    if (nz2 != A->p[n])
        printf("tiny entries dropped: %g\n",
            (double)(nz2 - A->p[n]));
#endif // PRINT
    Prob->b = cs_malloc(mn, sizeof(double));
    Prob->x = cs_malloc(mn, sizeof(double));
    Prob->resid = cs_malloc(mn, sizeof(double));
    // initialization
    for (int pos = 0; pos < mn; pos++) {
        Prob->b[pos] = 0;
        Prob->x[pos] = 0;
        Prob->resid[pos] = 0;
    }
    return ((!Prob->b || !Prob->x || !Prob->resid) ? free_problem(Prob) : Prob);
}

static double tic(void) { return 0; } //  (clock () / (double) CLOCKS_PER_SEC) ; }
static double toc(double t) { return 0; } // double s = tic () ; return (CS_MAX (0, s-t)) ; }

