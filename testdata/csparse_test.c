#include "cs.h"
int main(void)
{
	// load data
	cs *T = cs_load(stdin);

	// print data
	cs_print(T,0);

	// A compressed-column form of T
	cs *AC = cs_compress (T);

	// print data
	cs_print(AC,0);

	return 0;
}
