# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "eispack.h"

int main()
{
  timestamp ( );
  double a[4*4] = {
    5.0, 4.0, 1.0, 1.0,
    4.0, 5.0, 1.0, 1.0,
    1.0, 1.0, 4.0, 2.0,
    1.0, 1.0, 2.0, 4.0 };
  double a2[4*4];
  int i;
  int ierr;
  int j;
  int k;
  int matz;
  int n = 4;
  double w[4];
  double x[4*4];
/*
  Save a copy of the matrix.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a2[i+j*n] = a[i+j*n];
    }
  }
  printf ( "\n" );
  printf ( "TEST06 (KI)\n" );
  printf ( "  RS computes the eigenvalues and eigenvectors\n" );
  printf ( "  of a real symmetric matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order = %d\n", n );

  r8mat_print ( n, n, a, "  The matrix A:" );

  matz = 0;

  ierr = rs ( n, a, w, matz, x );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "TEST06 - Warning!\n" );
    printf ( "  The error return flag IERR = %d\n", ierr );
    return -1;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  timestamp ( );
  return 0;
}
