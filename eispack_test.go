package sparse

import (
	"fmt"
	"math"
	"os"
	"testing"
)

func TestEispack(t *testing.T) {
	// compile eispack

	// gcc -o eispack eispack_prb1.c eispack.c eispack.h -lm
	// ./eispack

	// generate output

	// run transpiled code

	// compare output
}

// r8mat_print - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1307
//
//
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title.
//
func r8mat_print(m int, n int, a []float64, title string) {
	r8mat_print_some(m, n, a, 1, 1, m, n, title)
}

// r8mat_print_some - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1351
//
//
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title.
//
func r8mat_print_some(m int, n int, a []float64, ilo int, jlo int, ihi int, jhi int, title string) {
	var i int
	var i2hi int
	var i2lo int
	var j int
	var j2hi int
	var j2lo int
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "%s\n", title)
	if m <= 0 || n <= 0 {
		fmt.Fprintf(os.Stdout, "\n")
		fmt.Fprintf(os.Stdout, "  (None)\n")
		return
	}
	{
		//
		//  Print the columns of the matrix, in strips of 5.
		//
		for j2lo = jlo; j2lo <= jhi; j2lo = j2lo + 5 {
			j2hi = j2lo + 5 - 1
			j2hi = i4_min(j2hi, n)
			j2hi = i4_min(j2hi, jhi)
			fmt.Fprintf(os.Stdout, "\n")
			//
			//  For each column J in the current range...
			//
			//  Write the header.
			//
			fmt.Fprintf(os.Stdout, "  Col:  ")
			for j = j2lo; j <= j2hi; j++ {
				fmt.Fprintf(os.Stdout, "  %7d     ", j-1)
			}
			fmt.Fprintf(os.Stdout, "\n")
			fmt.Fprintf(os.Stdout, "  Row\n")
			fmt.Fprintf(os.Stdout, "\n")
			//
			//  Determine the range of the rows in this strip.
			//
			i2lo = i4_max(ilo, 1)
			i2hi = i4_min(ihi, m)
			for i = i2lo; i <= i2hi; i++ {
				//
				//  Print out (up to) 5 entries in row I, that lie in the current strip.
				//
				fmt.Fprintf(os.Stdout, "%5d:", i-1)
				for j = j2lo; j <= j2hi; j++ {
					fmt.Fprintf(os.Stdout, "  %14f", a[i-1+(j-1)*m])
				}
				fmt.Fprintf(os.Stdout, "\n")
			}
		}
	}
}

// r8mat_uniform_01_new - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1459
//
//
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
func r8mat_uniform_01_new(m int, n int, seed *int) []float64 {
	var i int
	var j int
	var k int
	var r []float64
	r = make([]float64, uint32(m*n)*8*1/8)
	for j = 0; j < n; j++ {
		for i = 0; i < m; i++ {
			k = *seed / 127773
			*seed = 16807*(*seed-k*127773) - k*2836
			if *seed < 0 {
				*seed = *seed + 2147483647
			}
			r[i+j*m] = float64(*seed) * 4.656612875e-10
		}
	}
	return r
}

// r8vec_print - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1549
//
//
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, char *TITLE, a title.
//
func r8vec_print(n int, a []float64, title string) {
	var i int
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "%s\n", title)
	fmt.Fprintf(os.Stdout, "\n")
	for i = 0; i < n; i++ {
		fmt.Fprintf(os.Stdout, "  %8d: %14f\n", i, a[i])
	}
}

// transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:28
//
//  Purpose:
//
//    MAIN is the main program for EISPACK_PRB1.
//
//  Discussion:
//
//    EISPACK_PRB1 calls the EISPACK sample programs.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
func intenalEispackTest(t *testing.T) {
	// timestamp()
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "EISPACK_PRB\n")
	fmt.Fprintf(os.Stdout, "  C version.\n")
	fmt.Fprintf(os.Stdout, "  Test the EISPACK library.\n")
	//
	//  test01 ( );
	//  test02 ( );
	//  test03 ( );
	//  test04 ( );
	//  test05 ( );
	//
	test06()
	test065()
	test07()
	//
	//  test08 ( );
	//  test09 ( );
	//  test10 ( );
	//
	//  test11 ( );
	//  test12 ( );
	//  test13 ( );
	//  test14 ( );
	//  test15 ( );
	//  test16 ( );
	//
	//
	//  Terminate.
	//
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "EISPACK_PRB1\n")
	fmt.Fprintf(os.Stdout, "  Normal end of execution.\n")
	fmt.Fprintf(os.Stdout, "\n")
	// timestamp()
	return
}

// test06 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:95
//
func test06() {
	var a []float64 = []float64{5, 4, 1, 1, 4, 5, 1, 1, 1, 1, 4, 2, 1, 1, 2, 4}
	var a2 []float64 = make([]float64, 16)
	var i int
	var ierr int
	var j int
	// var k int
	var matz int
	var n int = 4
	var r []float64
	var w []float64 = make([]float64, 4)
	var x []float64 = make([]float64, 16)
	{
		//
		//
		//  Purpose:
		//
		//    TEST06 tests RS.
		//
		//  Licensing:
		//
		//    This code is distributed under the GNU LGPL license.
		//
		//  Modified:
		//
		//    08 November 2012
		//
		//  Author:
		//
		//    John Burkardt
		//
		//
		//  Save a copy of the matrix.
		//
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				a2[i+j*n] = a[i+j*n]
			}
		}
	}
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "TEST06\n")
	fmt.Fprintf(os.Stdout, "  RS computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(os.Stdout, "  of a real symmetric matrix.\n")
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a, "  The matrix A:")
	matz = 1
	ierr = rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(os.Stdout, "\n")
		fmt.Fprintf(os.Stdout, "TEST06 - Warning!\n")
		fmt.Fprintf(os.Stdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
	if matz != 0 {
		r8mat_print(n, n, x, "  The eigenvector matrix:")
		r = r8mat_mm_new(n, n, n, a2, x)
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				r[i+j*n] = r[i+j*n] - w[j]*x[i+j*n]
			}
		}
		r8mat_print(n, n, r, "  The residual (A-Lambda*I)*X:")
	}
	_ = r
}

// test065 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:187
//
//
//
//  Purpose:
//
//    TEST065 tests RS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
func test065() {
	var a []float64
	var a2 []float64 = make([]float64, 9)
	var i int
	var ierr int
	var j int
	// var k int
	var matz int
	var n int = 3
	var r []float64
	var seed int
	var t float64
	var w []float64 = make([]float64, 3)
	var x []float64 = make([]float64, 9)
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "TEST065\n")
	fmt.Fprintf(os.Stdout, "  RS computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(os.Stdout, "  of a real symmetric matrix.\n")
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "  Matrix order = %d\n", n)
	seed = 123456789
	a = r8mat_uniform_01_new(n, n, &seed)
	for i = 0; i < n-1; i++ {
		for j = i + 1; j < n; j++ {
			t = 0.5 * (a[i+j*n] + a[j+i*n])
			a[i+j*n] = t
			a[j+i*n] = t
		}
	}
	{
		//
		//  Save a copy of the matrix.
		//
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				a2[i+j*n] = a[i+j*n]
			}
		}
	}
	r8mat_print(n, n, a, "  The matrix A:")
	matz = 1
	ierr = rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(os.Stdout, "\n")
		fmt.Fprintf(os.Stdout, "TEST065 - Warning!\n")
		fmt.Fprintf(os.Stdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
	if matz != 0 {
		r8mat_print(n, n, x, "  The eigenvector matrix:")
		r = r8mat_mm_new(n, n, n, a2, x)
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				r[i+j*n] = r[i+j*n] - w[j]*x[i+j*n]
			}
		}
		r8mat_print(n, n, r, "  The residual (A-Lambda*I)*X:")
		_ = r
	}
	_ = a
}

// test07 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:293
//
//
//
//  Purpose:
//
//    TEST07 tests RSB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2012
//
//  Author:
//
//    John Burkardt
//
func test07() {
	var a []float64
	var a2 []float64
	var i int
	var ierr int
	var j int
	var matz int
	var mb int = 2
	var n int = 5
	var r []float64
	var w []float64
	var x []float64
	a = make([]float64, uint32(n*mb)*8*1/8)
	for j = 0; j < mb; j++ {
		for i = 0; i < n; i++ {
			a[i+j*n] = 0
		}
	}
	j = mb - 1
	for i = 0; i < n; i++ {
		a[i+j*n] = 2
	}
	j = 0
	for i = 1; i < n; i++ {
		a[i+j*n] = -1
	}
	a2 = make([]float64, uint32(n*n)*8*1/8)
	for j = 0; j < n; j++ {
		for i = 0; i < n; i++ {
			if i == j {
				a2[i+j*n] = 2
			} else if int(math.Abs(float64(i-j))) == 1 {
				a2[i+j*n] = -1
			} else {
				a2[i+j*n] = 0
			}
		}
	}
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "TEST07\n")
	fmt.Fprintf(os.Stdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(os.Stdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(os.Stdout, "\n")
	fmt.Fprintf(os.Stdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, uint32(n)*8*1/8)
	x = make([]float64, uint32(n*n)*8*1/8)
	matz = 1
	ierr = rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(os.Stdout, "\n")
		fmt.Fprintf(os.Stdout, "TEST07 - Warning!\n")
		fmt.Fprintf(os.Stdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
	if matz != 0 {
		r8mat_print(n, n, x, "  The eigenvector matrix X:")
		r = r8mat_mm_new(n, n, n, a2, x)
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				r[i+j*n] = r[i+j*n] - x[i+j*n]*w[j]
			}
		}
		r8mat_print(n, n, r, "  The residual (A-Lambda*I)*X:")
		_ = r
	}
	_ = a
	_ = a2
	_ = w
	_ = x
}
