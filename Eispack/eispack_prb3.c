# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "eispack.h"

int main()
{
  timestamp ( );

  double *a;
  double *a2;
  int i;
  int ierr;
  int j;
  int matz;
  int mb = 2;
  int n = 5;
  double *r;
  double *w;
  double *x;

  a = ( double * ) malloc ( n * mb * sizeof ( double ) );

  for ( j = 0; j < mb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = 0.0;
    }
  }
  j = mb - 1;
  for ( i = 0; i < n; i++ )
  {
    a[i+j*n] = 2.0;
  }
  j = 0;
  for ( i = 1; i < n; i++ )
  {
    a[i+j*n] = -1.0;
  }

  a2 = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a2[i+j*n] = 2.0;
      }
      else if ( abs ( i - j ) == 1 )
      {
        a2[i+j*n] = - 1.0;
      }
      else
      {
        a2[i+j*n] = 0.0;
      }
    }
  }

  printf ( "\n" );
  printf ( "TEST07 (KI)\n" );
  printf ( "  RSB computes the eigenvalues and eigenvectors\n" );
  printf ( "  of a real symmetric band matrix.\n" );
  printf ( "\n" );
  printf ( "  Matrix order = %d\n", n );

  r8mat_print ( n, n, a2, "  The matrix A:" );

  w = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * n * sizeof ( double ) );
  matz = 0;

  ierr = rsb ( n, mb, a, w, matz, x );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "TEST07 - Warning!\n" );
    printf ( "  The error return flag IERR = %d\n", ierr );
    return -1;
  }

  r8vec_print ( n, w, "  The eigenvalues Lambda:" );

  timestamp ( );
  return 0;
}
