package sparse

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"os/exec"
	"strings"
	"testing"
)

func runEispack(source string, t *testing.T) string {
	{
		// compile eispack
		//
		// build testdata application :
		//
		// clang -o eispack eispack_prb1.c eispack.c eispack.h -lm
		var args []string
		args = append(args, []string{
			source,
			"Eispack/eispack.c",
		}...)
		args = append(args, "-lm")
		args = append(args, "-o")
		args = append(args, "Eispack/eispack")

		cmd := exec.Command("clang", args...)
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		err := cmd.Run()
		if err != nil {
			t.Fatalf("cmd.Run() failed with %s.\n%s\n%s\n",
				err,
				stderr.String(),
				stdout.String(),
			)
		}
	}

	var out string
	{
		// ./eispack
		cmd := exec.Command("./Eispack/eispack")
		var stdout, stderr bytes.Buffer
		cmd.Stdout = &stdout
		cmd.Stderr = &stderr
		err := cmd.Run()
		if err != nil {
			t.Fatalf("cmd.Run() failed with %s.\n%s\n%s\n",
				err,
				stderr.String(),
				stdout.String(),
			)
		}
		out += stdout.String()
	}

	return out
}

func TestEispack(t *testing.T) {
	tcs := []struct {
		source string
		f      func()
	}{
		{
			source: "Eispack/eispack_prb1.c",
			f:      intenalEispack1,
		},
		{
			source: "Eispack/eispack_prb2.c",
			f:      intenalEispack2,
		},
		{
			source: "Eispack/eispack_prb3.c",
			f:      intenalEispack3,
		},
		{
			source: "Eispack/eispack_prb4.c",
			f:      intenalEispack4,
		},
		{
			source: "Eispack/eispack_prb5.c",
			f:      intenalEispack5,
		},
		{
			source: "Eispack/eispack_prb6.c",
			f:      intenalEispack6,
		},
		{
			source: "Eispack/eispack_prb7.c",
			f:      intenalEispack7,
		},
		{
			source: "Eispack/eispack_prb8.c",
			f:      intenalEispack8,
		},
	}

	for i := range tcs {
		t.Run(tcs[i].source, func(t *testing.T) {
			out := runEispack(tcs[i].source, t)

			// remove first and last line
			{
				lines := strings.Split(out, "\n")
				out = strings.Join(lines[1:len(lines)-2], "\n")
				out = strings.TrimSpace(out)
			}

			// generate output
			tmpfile, err := ioutil.TempFile("", "demo2")
			if err != nil {
				t.Fatal(err)
			}
			old := osStdout
			osStdout = tmpfile
			defer func() {
				osStdout = old
			}()

			// run transpiled code
			tcs[i].f()

			// compare output
			filename := tmpfile.Name()
			defer func() { _ = os.Remove(filename) }()
			err = tmpfile.Close()
			if err != nil {
				t.Fatal(err)
			}
			cb2, err := ioutil.ReadFile(filename)
			if err != nil {
				t.Fatal(err)
			}
			out2 := string(cb2)
			out2 = strings.TrimSpace(out2)

			if out != out2 {
				t.Fatal(ShowDiff(out, out2))
			}

			// t.Log(out2)
		})
	}
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "%s\n", title)
	if m <= 0 || n <= 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "  (None)\n")
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
			fmt.Fprintf(osStdout, "\n")
			//
			//  For each column J in the current range...
			//
			//  Write the header.
			//
			fmt.Fprintf(osStdout, "  Col:  ")
			for j = j2lo; j <= j2hi; j++ {
				fmt.Fprintf(osStdout, "  %7d     ", j-1)
			}
			fmt.Fprintf(osStdout, "\n")
			fmt.Fprintf(osStdout, "  Row\n")
			fmt.Fprintf(osStdout, "\n")
			//
			//  Determine the range of the rows in this strip.
			//
			i2lo = i4_max(ilo, 1)
			i2hi = i4_min(ihi, m)
			for i = i2lo; i <= i2hi; i++ {
				//
				//  Print out (up to) 5 entries in row I, that lie in the current strip.
				//
				fmt.Fprintf(osStdout, "%5d:", i-1)
				for j = j2lo; j <= j2hi; j++ {
					if math.Abs(a[i-1+(j-1)*m]) < 1e-8 {
						a[i-1+(j-1)*m] = 0.0
					}
					fmt.Fprintf(osStdout, "  %14f", a[i-1+(j-1)*m])
				}
				fmt.Fprintf(osStdout, "\n")
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
	r = make([]float64, m*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "%s\n", title)
	fmt.Fprintf(osStdout, "\n")
	for i = 0; i < n; i++ {
		fmt.Fprintf(osStdout, "  %8d: %14f\n", i, a[i])
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
func intenalEispack1() {
	// timestamp()
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "EISPACK_PRB\n")
	fmt.Fprintf(osStdout, "  C version.\n")
	fmt.Fprintf(osStdout, "  Test the EISPACK library.\n")
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "EISPACK_PRB1\n")
	fmt.Fprintf(osStdout, "  Normal end of execution.\n")
	fmt.Fprintf(osStdout, "\n")
	// timestamp()
	return
}

func intenalEispack2() {
	var a []float64 = []float64{5, 4, 1, 1, 4, 5, 1, 1, 1, 1, 4, 2, 1, 1, 2, 4}
	var a2 []float64 = make([]float64, 16)
	var i int
	var j int
	var matz int
	var n int = 4
	var w []float64 = make([]float64, 4)
	var x []float64 = make([]float64, 16)
	for j = 0; j < n; j++ {
		for i = 0; i < n; i++ {
			a2[i+j*n] = a[i+j*n]
		}
	}
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST06 (KI)\n")
	fmt.Fprintf(osStdout, "  RS computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a, "  The matrix A:")
	matz = 0
	ierr := rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST06 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

func intenalEispack3() {
	var a []float64
	var a2 []float64
	var i int
	// var ierr int
	var j int
	var matz int
	var mb int = 2
	var n int = 5
	var w []float64
	var x []float64
	a = make([]float64, n*mb)
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
	a2 = make([]float64, n*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST07 (KI)\n")
	fmt.Fprintf(osStdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, n)
	x = make([]float64, n*n)
	matz = 0
	ierr := rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST07 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

func intenalEispack4() {
	var a []float64 = []float64{15, 4, 1, 115, 4, 5, 5, 1, 1, 1, 4, 2, 1, 1, 2, 4}
	var a2 []float64 = make([]float64, 16)
	var i int
	// var ierr int
	var j int
	var matz int
	var n int = 4
	var w []float64 = make([]float64, 4)
	var x []float64 = make([]float64, 16)
	for j = 0; j < n; j++ {
		for i = 0; i < n; i++ {
			a2[i+j*n] = a[i+j*n]
		}
	}
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST06 (KI: not symmetrical)\n")
	fmt.Fprintf(osStdout, "  RS computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real not symmetric matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a, "  The matrix A:")
	matz = 0
	ierr := rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST06 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

func intenalEispack5() {
	var a []float64
	var a2 []float64
	var i int
	// var ierr int
	var j int
	var matz int
	var mb int = 1
	var n int = 5
	var w []float64
	var x []float64
	a = make([]float64, n*mb)
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
	a2 = make([]float64, n*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST07 (KI)\n")
	fmt.Fprintf(osStdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, n)
	x = make([]float64, n*n)
	matz = 0
	ierr := rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST07 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

func intenalEispack6() {
	var a []float64
	var a2 []float64
	var i int
	// var ierr int
	var j int
	var matz int
	var mb int = 6
	var n int = 100
	var w []float64
	var x []float64
	a = make([]float64, n*mb)
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
	a2 = make([]float64, n*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST07 (KI)\n")
	fmt.Fprintf(osStdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, n)
	x = make([]float64, n*n)
	matz = 0
	ierr := rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST07 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

func intenalEispack7() {
	var a []float64
	var a2 []float64
	var i int
	// var ierr int
	var j int
	var matz int
	var mb int = 6
	var n int = 100
	var w []float64
	var x []float64
	a = make([]float64, n*mb)
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
	a2 = make([]float64, n*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST07 (KI)\n")
	fmt.Fprintf(osStdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, n)
	x = make([]float64, n*n)
	matz = 1
	ierr := rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST07 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

func intenalEispack8() {
	var a []float64
	var a2 []float64
	var i int
	// var ierr int
	var j int
	var matz int
	var mb int = 6
	var n int = 200
	var w []float64
	var x []float64
	a = make([]float64, n*mb)
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
	a2 = make([]float64, n*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST07 (KI)\n")
	fmt.Fprintf(osStdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, n)
	x = make([]float64, n*n)
	matz = 1
	ierr := rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST07 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
		return
	}
	r8vec_print(n, w, "  The eigenvalues Lambda:")
}

// test06 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:95
//
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
func test06() {
	var a []float64 = []float64{5, 4, 1, 1, 4, 5, 1, 1, 1, 1, 4, 2, 1, 1, 2, 4}
	var a2 []float64 = make([]float64, 16)
	var i int
	// var ierr int
	var j int
	// var k int
	var matz int
	var n int = 4
	var r []float64
	var w []float64 = make([]float64, 4)
	var x []float64 = make([]float64, 16)
	for j = 0; j < n; j++ {
		for i = 0; i < n; i++ {
			a2[i+j*n] = a[i+j*n]
		}
	}
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST06\n")
	fmt.Fprintf(osStdout, "  RS computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a, "  The matrix A:")
	matz = 1
	ierr := rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST06 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
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
	// var ierr int
	var j int
	// var k int
	var matz int
	var n int = 3
	var r []float64
	var seed int
	var t float64
	var w []float64 = make([]float64, 3)
	var x []float64 = make([]float64, 9)
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST065\n")
	fmt.Fprintf(osStdout, "  RS computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
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
	ierr := rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST065 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
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
	// var ierr int
	var j int
	var matz int
	var mb int = 2
	var n int = 5
	var r []float64
	var w []float64
	var x []float64
	a = make([]float64, n*mb)
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
	a2 = make([]float64, n*n)
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
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "TEST07\n")
	fmt.Fprintf(osStdout, "  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Fprintf(osStdout, "  of a real symmetric band matrix.\n")
	fmt.Fprintf(osStdout, "\n")
	fmt.Fprintf(osStdout, "  Matrix order = %d\n", n)
	r8mat_print(n, n, a2, "  The matrix A:")
	w = make([]float64, n)
	x = make([]float64, n*n)
	matz = 1
	ierr := rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Fprintf(osStdout, "\n")
		fmt.Fprintf(osStdout, "TEST07 - Warning!\n")
		fmt.Fprintf(osStdout, "  The error return flag IERR = %d\n", ierr)
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

func BenchmarkWithOrWithoutEigenVector(b *testing.B) {
	var a []float64 = []float64{5, 4, 1, 1, 4, 5, 1, 1, 1, 1, 4, 2, 1, 1, 2, 4}
	var n int = 4
	var w []float64 = make([]float64, 4)
	var x []float64 = make([]float64, 16)

	a2 := make([]float64, 16)

	b.Run("with    eigenvector", func(b *testing.B) {
		matz := 1
		for i := 0; i < n; i++ {
			w[i] = 0.0
		}
		for i := 0; i < n*n; i++ {
			x[i] = 0.0
		}
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			for i := 0; i < n*n; i++ {
				a2[i] = a[i]
			}
			ierr := rs(n, a2, w, matz, x)
			if ierr != 0 {
				panic(fmt.Errorf("ierr = %v", ierr))
			}
		}
	})

	b.Run("without eigenvector", func(b *testing.B) {
		matz := 0
		for i := 0; i < n; i++ {
			w[i] = 0.0
		}
		for i := 0; i < n*n; i++ {
			x[i] = 0.0
		}
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			for i := 0; i < n*n; i++ {
				a2[i] = a[i]
			}
			ierr := rs(n, a2, w, matz, x)
			if ierr != 0 {
				panic(fmt.Errorf("ierr = %v", ierr))
			}
		}
	})
}
