#include "cs.h"
int main(void)
{
	// load data
	cs *T = cs_load(stdin);

	// print data
	cs_print(T,0);

	return 0;
}
