#include "test_function.h"

int main(void)
{
    // load data
    cs* T = cs_load(stdin);

    // A compressed-column form of T
    cs* A = cs_compress(T);

    // Transpose
    cs* AT = cs_transpose(A, 1);

    // m = # of rows of A
    int m = A ? A->m : 0;

    // triplet identify
    T = cs_spalloc(m, m, m, 1, 1);
    for (int i = 0; i < m; i++)
        cs_entry(T, i, i, 1);
    cs* Eye = cs_compress(T);

    // C = A*A'
    cs* C = cs_multiply(A, AT);

    // D = C + Eye*norm(C,1)
    cs* D = cs_add(C, Eye, 1, cs_norm(C));
    test_cs_print(D, 0);

    return 0;
}
