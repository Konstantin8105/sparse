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

	return 0;
}
