//
//	Package main - transpiled by c4go
//
//	If you have found any issues, please raise an issue at:
//	https://github.com/Konstantin8105/c4go/
//

package main

import "unsafe"
import "github.com/Konstantin8105/c4go/noarch"
import "fmt"

// main - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:28
//
func main() {
	defer noarch.AtexitRun()
	//
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
	timestamp()
	fmt.Printf("\n")
	fmt.Printf("EISPACK_PRB\n")
	fmt.Printf("  C version.\n")
	fmt.Printf("  Test the EISPACK library.\n")
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
	fmt.Printf("\n")
	fmt.Printf("EISPACK_PRB1\n")
	fmt.Printf("  Normal end of execution.\n")
	fmt.Printf("\n")
	timestamp()
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
	var k int
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
	fmt.Printf("\n")
	fmt.Printf("TEST06\n")
	fmt.Printf("  RS computes the eigenvalues and eigenvectors\n")
	fmt.Printf("  of a real symmetric matrix.\n")
	fmt.Printf("\n")
	noarch.Printf([]byte("  Matrix order = %d\n\x00"), n)
	r8mat_print(n, n, a, []byte("  The matrix A:\x00"))
	matz = 1
	ierr = rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Printf("\n")
		fmt.Printf("TEST06 - Warning!\n")
		noarch.Printf([]byte("  The error return flag IERR = %d\n\x00"), ierr)
		return
	}
	r8vec_print(n, w, []byte("  The eigenvalues Lambda:\x00"))
	if matz != 0 {
		r8mat_print(n, n, x, []byte("  The eigenvector matrix:\x00"))
		r = r8mat_mm_new(n, n, n, a2, x)
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				r[i+j*n] = r[i+j*n] - w[j]*x[i+j*n]
			}
		}
		r8mat_print(n, n, r, []byte("  The residual (A-Lambda*I)*X:\x00"))
	}
	_ = r
}

// test065 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:187
//
func test065() {
	var a []float64
	var a2 []float64 = make([]float64, 9)
	var i int
	var ierr int
	var j int
	var k int
	var matz int
	var n int = 3
	var r []float64
	var seed int
	var t float64
	var w []float64 = make([]float64, 3)
	var x []float64 = make([]float64, 9)
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
	fmt.Printf("\n")
	fmt.Printf("TEST065\n")
	fmt.Printf("  RS computes the eigenvalues and eigenvectors\n")
	fmt.Printf("  of a real symmetric matrix.\n")
	fmt.Printf("\n")
	noarch.Printf([]byte("  Matrix order = %d\n\x00"), n)
	seed = 123456789
	a = r8mat_uniform_01_new(n, n, (*[100000000]int)(unsafe.Pointer(&seed))[:])
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
	r8mat_print(n, n, a, []byte("  The matrix A:\x00"))
	matz = 1
	ierr = rs(n, a, w, matz, x)
	if ierr != 0 {
		fmt.Printf("\n")
		fmt.Printf("TEST065 - Warning!\n")
		noarch.Printf([]byte("  The error return flag IERR = %d\n\x00"), ierr)
		return
	}
	r8vec_print(n, w, []byte("  The eigenvalues Lambda:\x00"))
	if matz != 0 {
		r8mat_print(n, n, x, []byte("  The eigenvector matrix:\x00"))
		r = r8mat_mm_new(n, n, n, a2, x)
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				r[i+j*n] = r[i+j*n] - w[j]*x[i+j*n]
			}
		}
		r8mat_print(n, n, r, []byte("  The residual (A-Lambda*I)*X:\x00"))
		_ = r
	}
	_ = a
}

// test07 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack_prb1.c:293
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
			} else if noarch.Abs(i-j) == 1 {
				a2[i+j*n] = -1
			} else {
				a2[i+j*n] = 0
			}
		}
	}
	fmt.Printf("\n")
	fmt.Printf("TEST07\n")
	fmt.Printf("  RSB computes the eigenvalues and eigenvectors\n")
	fmt.Printf("  of a real symmetric band matrix.\n")
	fmt.Printf("\n")
	noarch.Printf([]byte("  Matrix order = %d\n\x00"), n)
	r8mat_print(n, n, a2, []byte("  The matrix A:\x00"))
	w = make([]float64, uint32(n)*8*1/8)
	x = make([]float64, uint32(n*n)*8*1/8)
	matz = 1
	ierr = rsb(n, mb, a, w, matz, x)
	if ierr != 0 {
		fmt.Printf("\n")
		fmt.Printf("TEST07 - Warning!\n")
		noarch.Printf([]byte("  The error return flag IERR = %d\n\x00"), ierr)
		return
	}
	r8vec_print(n, w, []byte("  The eigenvalues Lambda:\x00"))
	if matz != 0 {
		r8mat_print(n, n, x, []byte("  The eigenvector matrix X:\x00"))
		r = r8mat_mm_new(n, n, n, a2, x)
		for j = 0; j < n; j++ {
			for i = 0; i < n; i++ {
				r[i+j*n] = r[i+j*n] - x[i+j*n]*w[j]
			}
		}
		r8mat_print(n, n, r, []byte("  The residual (A-Lambda*I)*X:\x00"))
		_ = r
	}
	_ = a
	_ = a2
	_ = w
	_ = x
}