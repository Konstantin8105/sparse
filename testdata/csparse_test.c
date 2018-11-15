#include "cs.h"
int main(void)
{
	// load data
	cs *T = cs_load(stdin);

	// print data
	cs_print(T,0);

	// A compressed-column form of T
	cs *A = cs_compress (T);

	// print data
	cs_print(A,0);

	// Transpose
	cs * AT = cs_transpose (A,1);

	// print data
	cs_print(AT,0);

	// m = # of rows of A
    int m = A ? A->m : 0 ;

	// triplet identify
    T = cs_spalloc (m, m, m, 1, 1) ;
    for (int i = 0 ; i < m ; i++) cs_entry (T, i, i, 1) ;
    cs *Eye = cs_compress (T) ;
	cs_print(Eye,0);


	return 0;
}
