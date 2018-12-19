package sparse

import (
	"math"
)

// bakvec - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:10
//
//
//
//  Purpose:
//
//    BAKVEC determines eigenvectors by reversing the FIGI transformation.
//
//  Discussion:
//
//    This subroutine forms the eigenvectors of a nonsymmetric tridiagonal
//    matrix by back transforming those of the corresponding symmetric
//    matrix determined by FIGI.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double T[N*3], contains the nonsymmetric matrix.  Its
//    subdiagonal is stored in the positions 2:N of the first column,
//    its diagonal in positions 1:N of the second column,
//    and its superdiagonal in positions 1:N-1 of the third column.
//    T(1,1) and T(N,3) are arbitrary.
//
//    Input/output, double E[N].  On input, E(2:N) contains the
//    subdiagonal elements of the symmetric matrix.  E(1) is arbitrary.
//    On output, the contents of E have been destroyed.
//
//    Input, int M, the number of eigenvectors to be back
//    transformed.
//
//    Input/output, double Z[N*M], contains the eigenvectors.
//    On output, they have been transformed as requested.
//
//    Output, int BAKVEC, an error flag.
//    0, for normal return,
//    2*N+I, if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
//    In this case, the symmetric matrix is not similar
//    to the original matrix, and the eigenvectors
//    cannot be found by this program.
//
// func bakvec(n int, t []float64, e []float64, m int, z []float64) int {
// 	var i int
// 	var ierr int
// 	var j int
// 	ierr = 0
// 	if m == 0 {
// 		return ierr
// 	}
// 	e[0] = 1
// 	if n == 1 {
// 		return ierr
// 	}
// 	for i = 1; i < n; i++ {
// 		if e[i] == 0 {
// 			if t[i+0*3] != 0 || t[i-1+2*3] != 0 {
// 				ierr = 2*n + (i + 1)
// 				return ierr
// 			}
// 			e[i] = 1
// 		} else {
// 			e[i] = e[i-1] * e[i] / t[i-1+2*3]
// 		}
// 	}
// 	for j = 0; j < m; j++ {
// 		for i = 1; i < n; i++ {
// 			z[i+j*n] = z[i+j*n] * e[i]
// 		}
// 	}
// 	return ierr
// }

// balbak - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:129
//
//
//  Purpose:
//
//    BALBAK determines eigenvectors by undoing the BALANC transformation.
//
//  Discussion:
//
//    This subroutine forms the eigenvectors of a real general matrix by
//    back transforming those of the corresponding balanced matrix
//    determined by BALANC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    Parlett and Reinsch,
//    Numerische Mathematik,
//    Volume 13, pages 293-304, 1969.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int LOW, IGH, column indices determined by BALANC.
//
//    Input, double SCALE[N], contains information determining
//    the permutations and scaling factors used by BALANC.
//
//    Input, int M, the number of columns of Z to be
//    back-transformed.
//
//    Input/output, double Z[N*M], contains the real and imaginary
//    parts of the eigenvectors, which, on return, have been back-transformed.
//
// func balbak(n int, low int, igh int, scale []float64, m int, z []float64) {
// 	var i int
// 	var ii int
// 	var j int
// 	var k int
// 	// var s float64
// 	var t float64
// 	if m <= 0 {
// 		return
// 	}
// 	if igh != low {
// 		for i = low - 1; i <= igh-1; i++ {
// 			for j = 0; j < m; j++ {
// 				z[i+j*n] = z[i+j*n] * scale[i]
// 			}
// 		}
// 	}
// 	for ii = 1; ii <= n; ii++ {
// 		i = ii
// 		if i < low || igh < i {
// 			if i < low {
// 				i = low - ii
// 			}
// 			k = int((scale[i-1]))
// 			if k != i {
// 				for j = 0; j < m; j++ {
// 					t = z[i-1+j*n]
// 					z[i-1+j*n] = z[k-1+j*n]
// 					z[k-1+j*n] = t
// 				}
// 			}
// 		}
// 	}
// }

// bandr - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:246
//
//
//
//  Purpose:
//
//    BANDR reduces a symmetric band matrix to symmetric tridiagonal form.
//
//  Discussion:
//
//    This subroutine reduces a real symmetric band matrix
//    to a symmetric tridiagonal matrix using and optionally
//    accumulating orthogonal similarity transformations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MB, is the (half) band width of the matrix,
//    defined as the number of adjacent diagonals, including the principal
//    diagonal, required to specify the non-zero portion of the
//    lower triangle of the matrix.
//
//    Input/output, double A[N*MB].  On input, contains the lower
//    triangle of the symmetric band input matrix stored as an N by MB array.
//    Its lowest subdiagonal is stored in the last N+1-MB positions of the first
//    column, its next subdiagonal in the last N+2-MB positions of the second
//    column, further subdiagonals similarly, and finally its principal diagonal
//    in the N positions of the last column.  Contents of storages not part of
//    the matrix are arbitrary.  On output, A has been destroyed, except for
//    its last two columns which contain a copy of the tridiagonal matrix.
//
//    Output, double D[N], the diagonal elements of the tridiagonal
//    matrix.
//
//    Output, double E[N], the subdiagonal elements of the tridiagonal
//    matrix in E(2:N).  E(1) is set to zero.
//
//    Output, double E2[N], contains the squares of the corresponding
//    elements of E.  E2 may coincide with E if the squares are not needed.
//
//    Input, logical MATZ, should be set to TRUE if the transformation matrix is
//    to be accumulated, and to FALSE otherwise.
//
//    Output, double Z[N*N], the orthogonal transformation matrix
//    produced in the reduction if MATZ has been set to TRUE.  Otherwise, Z is
//    not referenced.
//
func bandr(n int, mb int, a []float64, d []float64, e []float64, e2 []float64, matz int, z []float64) {
	var b1 float64
	var b2 float64
	var c2 float64
	var dmin float64
	var dminrt float64
	var f1 float64
	var f2 float64
	var g float64
	var i int
	var i1 int
	var i2 int
	var j int
	var j1 int
	var j2 int
	var jj int
	var k int
	var kr int
	var l int
	var m1 int
	var maxl int
	var maxr int
	var mr int
	var r int
	var r1 int
	var s2 float64
	var u float64
	var ugl int
	dmin = r8_epsilon()
	dminrt = math.Sqrt(dmin)
	{
		//
		//  Initialize the diagonal scaling matrix.
		//
		for i = 0; i < n; i++ {
			d[i] = 1
		}
	}
	if matz != 0 {
		r8mat_identity(n, z)
	}
	if mb == 1 {
		{
			//
			//  Is input matrix diagonal?
			//
			for i = 0; i < n; i++ {
				d[i] = a[i+(mb-1)*n]
				e[i] = 0
				e2[i] = 0
			}
		}
		return
	}
	m1 = mb - 1
	if m1 != 1 {
		for k = 1; k <= n-2; k++ {
			maxr = i4_min(m1, n-k)
			for r1 = 2; r1 <= maxr; r1++ {
				r = maxr + 2 - r1
				kr = k + r
				mr = mb - r
				g = a[kr-1+(mr-1)*n]
				a[kr-2+0*n] = a[kr-2+mr*n]
				ugl = k
				for j = kr; j <= n; j = j + m1 {
					j1 = j - 1
					j2 = j1 - 1
					if g == 0 {
						break
					}
					b1 = a[j1-1+0*n] / g
					b2 = b1 * d[j1-1] / d[j-1]
					s2 = 1 / (1 + b1*b2)
					if s2 < 0.5 {
						b1 = g / a[j1-1+0*n]
						b2 = b1 * d[j-1] / d[j1-1]
						c2 = 1 - s2
						d[j1-1] = c2 * d[j1-1]
						d[j-1] = c2 * d[j-1]
						f1 = 2 * a[j-1+(m1-1)*n]
						f2 = b1 * a[j1-1+(mb-1)*n]
						a[j-1+(m1-1)*n] = -b2*(b1*a[j-1+(m1-1)*n]-a[j-1+(mb-1)*n]) - f2 + a[j-1+(m1-1)*n]
						a[j1-1+(mb-1)*n] = b2*(b2*a[j-1+(mb-1)*n]+f1) + a[j1-1+(mb-1)*n]
						a[j-1+(mb-1)*n] = b1*(f2-f1) + a[j-1+(mb-1)*n]
						for l = ugl; l <= j2; l++ {
							i2 = mb - j + l
							u = a[j1-1+i2*n] + b2*a[j-1+(i2-1)*n]
							a[j-1+(i2-1)*n] = -b1*a[j1-1+i2*n] + a[j-1+(i2-1)*n]
							a[j1-1+i2*n] = u
						}
						ugl = j
						a[j1-1+0*n] = a[j1+0*n] + b2*g
						if j != n {
							maxl = i4_min(m1, n-j1)
							for l = 2; l <= maxl; l++ {
								i1 = j1 + l
								i2 = mb - l
								u = a[i1-1+(i2-1)*n] + b2*a[i1-1+i2*n]
								a[i1-1+i2*n] = -b1*a[i1-1+(i2-1)*n] + a[i1-1+i2*n]
								a[i1-1+(i2-1)*n] = u
							}
							i1 = j + m1
							if i1 <= n {
								g = b2 * a[i1-1+0*n]
							}
						}
						if matz != 0 {
							for l = 1; l <= n; l++ {
								u = z[l-1+(j1-1)*n] + b2*z[l-1+(j-1)*n]
								z[l-1+(j-1)*n] = -b1*z[l-1+(j1-1)*n] + z[l-1+(j-1)*n]
								z[l-1+(j1-1)*n] = u
							}
						}
					} else {
						u = d[j1-1]
						d[j1-1] = s2 * d[j-1]
						d[j-1] = s2 * u
						f1 = 2 * a[j-1+(m1-1)*n]
						f2 = b1 * a[j-1+(mb-1)*n]
						u = b1*(f2-f1) + a[j1-1+(mb-1)*n]
						a[j-1+(m1-1)*n] = b2*(b1*a[j-1+(m1-1)*n]-a[j1-1+(mb-1)*n]) + f2 - a[j-1+(m1-1)*n]
						a[j1-1+(mb-1)*n] = b2*(b2*a[j1-1+(mb-1)*n]+f1) + a[j-1+(mb-1)*n]
						a[j-1+(mb-1)*n] = u
						for l = ugl; l <= j2; l++ {
							i2 = mb - j + l
							u = b2*a[j1-1+i2*n] + a[j-1+(i2-1)*n]
							a[j-1+(i2-1)*n] = -a[j1-1+i2*n] + b1*a[j-1+(i2-1)*n]
							a[j1-1+i2*n] = u
						}
						ugl = j
						a[j1-1+0*n] = b2*a[j1-1+0*n] + g
						if j != n {
							maxl = i4_min(m1, n-j1)
							for l = 2; l <= maxl; l++ {
								i1 = j1 + l
								i2 = mb - l
								u = b2*a[i1-1+(i2-1)*n] + a[i1-1+i2*n]
								a[i1-1+i2*n] = -a[i1-1+(i2-1)*n] + b1*a[i1-1+i2*n]
								a[i1-1+(i2-1)*n] = u
							}
							i1 = j + m1
							if i1 <= n {
								g = a[i1-1+0*n]
								a[i1-1+0*n] = b1 * a[i1-1+0*n]
							}
						}
						if matz != 0 {
							for l = 1; l <= n; l++ {
								u = b2*z[l-1+(j1-1)*n] + z[l-1+(j-1)*n]
								z[l-1+(j-1)*n] = -z[l-1+(j1-1)*n] + b1*z[l-1+(j-1)*n]
								z[l-1+(j1-1)*n] = u
							}
						}
					}
				}
			}
			if k%64 == 0 {
				{
					//
					//  Rescale to avoid underflow or overflow.
					//
					for j = k; j <= n; j++ {
						if d[j-1] < dmin {
							maxl = i4_max(1, mb+1-j)
							for jj = maxl; jj <= m1; jj++ {
								a[j-1+(jj-1)*n] = dminrt * a[j-1+(jj-1)*n]
							}
							if j != n {
								maxl = i4_min(m1, n-j)
								for l = 1; l <= maxl; l++ {
									i1 = j + l
									i2 = mb - l
									a[i1-1+(i2-1)*n] = dminrt * a[i1-1+(i2-1)*n]
								}
							}
							if matz != 0 {
								for i = 1; i <= n; i++ {
									z[i-1+(j-1)*n] = dminrt * z[i-1+(j-1)*n]
								}
							}
							a[j-1+(mb-1)*n] = dmin * a[j-1+(mb-1)*n]
							d[j-1] = d[j-1] / dmin
						}
					}
				}
			}
		}
	}
	{
		//
		//  Form square root of scaling matrix.
		//
		for i = 1; i < n; i++ {
			e[i] = math.Sqrt(d[i])
		}
	}
	if matz != 0 {
		for j = 1; j < n; j++ {
			for i = 0; i < n; i++ {
				z[i+j*n] = z[i+j*n] * e[j]
			}
		}
	}
	u = 1
	for j = 1; j < n; j++ {
		a[j+(m1-1)*n] = u * e[j] * a[j+(m1-1)*n]
		u = e[j]
		e2[j] = a[j+(m1-1)*n] * a[j+(m1-1)*n]
		a[j+(mb-1)*n] = d[j] * a[j+(mb-1)*n]
		d[j] = a[j+(mb-1)*n]
		e[j] = a[j+(m1-1)*n]
	}
	d[0] = a[0+(mb-1)*n]
	e[0] = 0
	e2[0] = 0
}

// cbabk2 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:606
//
//
//
//  Purpose:
//
//    CBABK2 finds eigenvectors by undoing the CBAL transformation.
//
//  Discussion:
//
//    This subroutine forms the eigenvectors of a complex general
//    matrix by back transforming those of the corresponding
//    balanced matrix determined by CBAL.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int LOW, IGH, values determined by CBAL.
//
//    Input, double SCALE[N], information determining the permutations
//    and scaling factors used by CBAL.
//
//    Input, int M, the number of eigenvectors to be back
//    transformed.
//
//    Input/output, double ZR[N*M], ZI[N*M].  On input, the real
//    and imaginary parts, respectively, of the eigenvectors to be back
//    transformed in their first M columns.  On output, the transformed
//    eigenvectors.
//
// func cbabk2(n int, low int, igh int, scale []float64, m int, zr []float64, zi []float64) {
// 	var i int
// 	var ii int
// 	var j int
// 	var k int
// 	var s float64
// 	if m == 0 {
// 		return
// 	}
// 	if igh != low {
// 		for i = low; i <= igh; i++ {
// 			s = scale[i]
// 			for j = 0; j < m; j++ {
// 				zr[i+j*n] = zr[i+j*n] * s
// 				zi[i+j*n] = zi[i+j*n] * s
// 			}
// 		}
// 	}
// 	for ii = 0; ii < n; ii++ {
// 		i = ii
// 		if i < low || igh < i {
// 			if i < low {
// 				i = low - ii
// 			}
// 			k = int(scale[i])
// 			if k != i {
// 				for j = 0; j < m; j++ {
// 					s = zr[i+j*n]
// 					zr[i+j*n] = zr[k+j*n]
// 					zr[k+j*n] = s
// 					s = zi[i+j*n]
// 					zi[i+j*n] = zi[k+j*n]
// 					zi[k+j*n] = s
// 				}
// 			}
// 		}
// 	}
// }

// csroot - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:724
//
//
//
//  Purpose:
//
//    CSROOT computes the complex square root of a complex quantity.
//
//  Discussion:
//
//    The branch of the square function is chosen so that
//      0.0 <= YR
//    and
//      sign ( YI ) == sign ( XI )
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, double XR, XI, the real and imaginary parts of the
//    quantity whose square root is desired.
//
//    Output, double *YR, *YI, the real and imaginary parts of the
//    square root.
//
// func csroot(xr float64, xi float64, yr []float64, yi []float64) {
// 	var (
// 		tr = xr
// 		ti = xi
// 		s  = math.Sqrt(0.5 * (pythag(tr, ti) + math.Abs(tr)))
// 	)
// 	if 0 <= tr {
// 		yr[0] = s
// 	}
// 	if ti < 0 {
// 		s = -s
// 	}
// 	if tr <= 0 {
// 		yi[0] = s
// 	}
// 	if tr < 0 {
// 		yr[0] = 0.5 * (ti / yi[0])
// 	} else if 0 < tr {
// 		yi[0] = 0.5 * (ti / yr[0])
// 	}
// }

// i4_max - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:814
//
//
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
func i4_max(i1 int, i2 int) int {
	var value int
	if i2 < i1 {
		value = i1
	} else {
		value = i2
	}
	return value
}

// i4_min - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:855
//
//
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
func i4_min(i1 int, i2 int) int {
	var value int
	if i1 < i2 {
		value = i1
	} else {
		value = i2
	}
	return value
}

// pythag - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:896
//
//
//
//  Purpose:
//
//    PYTHAG computes SQRT ( A * A + B * B ) carefully.
//
//  Discussion:
//
//    The formula
//
//      PYTHAG = sqrt ( A * A + B * B )
//
//    is reasonably accurate, but can fail if, for example, A^2 is larger
//    than the machine overflow.  The formula can lose most of its accuracy
//    if the sum of the squares is very large or very small.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Modified:
//
//    08 November 2012
//
//  Parameters:
//
//    Input, double A, B, the two legs of a right triangle.
//
//    Output, double PYTHAG, the length of the hypotenuse.
//
func pythag(a float64, b float64) (p float64) {
	p = math.Max(math.Abs(a), math.Abs(b))
	if p != 0 {
		r := math.Min(math.Abs(a), math.Abs(b)) / p
		r = r * r
		for {
			t := 4 + r
			if t == 4 {
				break
			}
			s := r / t
			u := 1 + 2*s
			p = u * p
			r = s / u * (s / u) * r
		}
	}
	return p
}

// r8_abs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:988
// removed

// r8_epsilon - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1029
//
//
//
//  Purpose:
//
//    R8_EPSILON returns the R8 round off unit.
//
//  Discussion:
//
//    R8_EPSILON is a number R which is a power of 2 with the property that,
//    to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
func r8_epsilon() float64 {
	return 2.220446049250313e-16
}

// r8_max - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1068
// removed

// r8_min - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1109
// removed

// r8_sign - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1150
//
//
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
func r8_sign(x float64) float64 {
	var value float64
	if x < 0 {
		value = -1
	} else {
		value = +1
	}
	return value
}

// r8mat_identity - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1191
//
//
//
//  Purpose:
//
//    R8MAT_IDENTITY sets an R8MAT to the identity matrix.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, double A[N*N], the N by N identity matrix.
//
func r8mat_identity(n int, a []float64) {
	k := 0
	for j := 0; j < n; j++ {
		for i := 0; i < n; i++ {
			if i == j {
				a[k] = 1
			} else {
				a[k] = 0
			}
			k = k + 1
		}
	}
}

// r8mat_mm_new - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1248
//
//
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
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
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
//
func r8mat_mm_new(n1, n2, n3 int, a, b []float64) []float64 {
	c := make([]float64, n1*n3)
	for i := 0; i < n1; i++ {
		for j := 0; j < n3; j++ {
			c[i+j*n1] = 0.0
			for k := 0; k < n2; k++ {
				c[i+j*n1] = c[i+j*n1] + a[i+k*n1]*b[k+j*n2]
			}
		}
	}
	return c
}

// rs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1596
//
//
//
//  Purpose:
//
//    RS computes eigenvalues and eigenvectors of real symmetric matrix.
//
//  Discussion:
//
//    This subroutine calls the recommended sequence of
//    subroutines from the eigensystem subroutine package (eispack)
//    to find the eigenvalues and eigenvectors (if desired)
//    of a real symmetric matrix.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the real symmetric matrix.
//
//    Input, int MATZ, is zero if only eigenvalues are desired,
//    and nonzero if both eigenvalues and eigenvectors are desired.
//
//    Output, double W[N], the eigenvalues in ascending order.
//
//    Output, double Z[N*N], contains the eigenvectors, if MATZ
//    is nonzero.
//
//    Output, int RS, is set equal to an error
//    completion code described in the documentation for TQLRAT and TQL2.
//    The normal completion code is zero.
//
func rs(n int, a []float64, w []float64, matz int, z []float64) (ierr int) {
	if matz == 0 {
		fv1 := make([]float64, n)
		fv2 := make([]float64, n)
		tred1(n, a, w, fv1, fv2)
		ierr = tqlrat(n, w, fv2)
		return
	}

	fv1 := make([]float64, n)
	tred2(n, a, w, fv1, z)
	ierr = tql2(n, w, fv1, z)

	return
}

// rsb - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1690
//
//
//
//  Purpose:
//
//    RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.
//
//  Discussion:
//
//    This subroutine calls the recommended sequence of
//    subroutines from the eigensystem subroutine package (eispack)
//    to find the eigenvalues and eigenvectors (if desired)
//    of a real symmetric band matrix.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, int MB, the half band width of the matrix,
//    defined as the number of adjacent diagonals, including the principal
//    diagonal, required to specify the non-zero portion of the lower triangle
//    of the matrix.
//
//    Input, double A[N*MB], contains the lower triangle of the real
//    symmetric band matrix.  Its lowest subdiagonal is stored in the last N+1-MB
//    positions of the first column, its next subdiagonal in the last
//    N+2-MB positions of the second column, further subdiagonals similarly,
//    and finally its principal diagonal in the N positions of the last
//    column.  Contents of storages not part of the matrix are arbitrary.
//
//    Input, int MATZ, is zero if only eigenvalues are desired,
//    and nonzero if both eigenvalues and eigenvectors are desired.
//
//    Output, double W[N], the eigenvalues in ascending order.
//
//    Output, double Z[N*N], contains the eigenvectors, if MATZ
//    is nonzero.
//
//    Output, int BANDR, is set to an error
//    completion code described in the documentation for TQLRAT and TQL2.
//    The normal completion code is zero.
//
func rsb(n, mb int, a, w []float64, matz int, z []float64) (ierr int) {
	if mb <= 0 {
		ierr = 12 * n
		return ierr
	}
	if n < mb {
		ierr = 12 * n
		return ierr
	}
	if matz == 0 {
		fv1 := make([]float64, n)
		fv2 := make([]float64, n)
		tf := 0
		bandr(n, mb, a, w, fv1, fv2, tf, z)
		ierr = tqlrat(n, w, fv2)
		return
	}

	fv1 := make([]float64, n)
	tf := 1
	bandr(n, mb, a, w, fv1, fv1, tf, z)
	ierr = tql2(n, w, fv1, z)

	return ierr
}

//  // timestamp - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1810
// 	//
// 	//
// 	//  Purpose:
// 	//
// 	//    TIMESTAMP prints the current YMDHMS date as a time stamp.
// 	//
// 	//  Example:
// 	//
// 	//    31 May 2001 09:45:54 AM
// 	//
// 	//  Licensing:
// 	//
// 	//    This code is distributed under the GNU LGPL license.
// 	//
// 	//  Modified:
// 	//
// 	//    24 September 2003
// 	//
// 	//  Author:
// 	//
// 	//    John Burkardt
// 	//
// 	//  Parameters:
// 	//
// 	//    None
// 	//
// func timestamp() {
// 	var time_buffer []byte = make([]byte, 40)
// 	var tm []noarch.Tm
// 	var len uint
// 	var now noarch.TimeT
// 	now = noarch.Time(nil)
// 	tm = noarch.LocalTime((*[100000000]noarch.TimeT)(unsafe.Pointer(&now))[:])
// 	len = strftime(time_buffer, 40, "%d %B %Y %I:%M:%S %p"), tm)
// 	fmt.Fprintf(os.Stdout, "%s\n", time_buffer)
// }

// tql2 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1858
//
//
//
//  Purpose:
//
//    TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This subroutine finds the eigenvalues and eigenvectors of a symmetric
//    tridiagonal matrix by the QL method.  The eigenvectors of a full
//    symmetric matrix can also be found if TRED2 has been used to reduce this
//    full matrix to tridiagonal form.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    Bowdler, Martin, Reinsch, Wilkinson,
//    TQL2,
//    Numerische Mathematik,
//    Volume 11, pages 293-306, 1968.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N].  On input, the diagonal elements of
//    the matrix.  On output, the eigenvalues in ascending order.  If an error
//    exit is made, the eigenvalues are correct but unordered for indices
//    1,2,...,IERR-1.
//
//    Input/output, double E[N].  On input, E(2:N) contains the
//    subdiagonal elements of the input matrix, and E(1) is arbitrary.
//    On output, E has been destroyed.
//
//    Input, double Z[N*N].  On input, the transformation matrix
//    produced in the reduction by TRED2, if performed.  If the eigenvectors of
//    the tridiagonal matrix are desired, Z must contain the identity matrix.
//    On output, Z contains the orthonormal eigenvectors of the symmetric
//    tridiagonal (or full) matrix.  If an error exit is made, Z contains
//    the eigenvectors associated with the stored eigenvalues.
//
//    Output, int TQL2, error flag.
//    0, normal return,
//    J, if the J-th eigenvalue has not been determined after
//    30 iterations.
//
func tql2(n int, d []float64, e []float64, z []float64) (ierr int) {
	if n == 1 {
		return ierr
	}

	for i := 1; i < n; i++ {
		e[i-1] = e[i]
	}

	var (
		f    = 0.0
		tst1 = 0.0
	)
	e[n-1] = 0
	for l := 0; l < n; l++ {
		j := 0
		h := math.Abs(d[l]) + math.Abs(e[l])
		tst1 = math.Max(tst1, h)

		//
		//  Look for a small sub-diagonal element.
		//
		var m int
		for m = l; m < n; m++ {
			tst2 := tst1 + math.Abs(e[m])
			if tst2 == tst1 {
				break
			}
		}

		if m != l {
			for {
				if 30 <= j {
					ierr = l + 1
					return ierr
				}
				j = j + 1
				//
				//  Form shift.
				//
				l1 := l + 1
				l2 := l1 + 1
				g := d[l]
				p := (d[l1] - g) / (2 * e[l])
				r := pythag(p, 1)
				d[l] = e[l] / (p + r8_sign(p)*math.Abs(r))
				d[l1] = e[l] * (p + r8_sign(p)*math.Abs(r))
				dl1 := d[l1]
				h = g - d[l]
				for i := l2; i < n; i++ {
					d[i] = d[i] - h
				}
				f = f + h
				//
				//  QL transformation.
				//
				p = d[m]
				var (
					c   = 1.0
					c2  = c
					el1 = e[l1]
					s   = 0.0
					s2  = 0.0
					c3  = 0.0
				)
				for ii := 1; ii <= m-l; ii++ {
					c3 = c2
					c2 = c
					s2 = s
					i := m - ii
					g = c * e[i]
					h = c * p
					r = pythag(p, e[i])
					e[i+1] = s * r
					s = e[i] / r
					c = p / r
					p = c*d[i] - s*g
					d[i+1] = h + s*(c*g+s*d[i])

					//
					//  Form vector.
					//
					for k := 0; k < n; k++ {
						h = z[k+(i+1)*n]
						z[k+(i+1)*n] = s*z[k+i*n] + c*h
						z[k+i*n] = c*z[k+i*n] - s*h
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1
				e[l] = s * p
				d[l] = c * p
				tst2 := tst1 + math.Abs(e[l])
				if tst2 <= tst1 {
					break
				}
			}
		}
		d[l] += f
	}

	//
	//  Order eigenvalues and eigenvectors.
	//
	for ii := 1; ii < n; ii++ {
		i := ii - 1
		k := i
		p := d[i]
		for j := ii; j < n; j++ {
			if d[j] < p {
				k = j
				p = d[j]
			}
		}
		if k != i {
			d[k] = d[i]
			d[i] = p
			for j := 0; j < n; j++ {
				// swap
				z[j+i*n], z[j+k*n] = z[j+k*n], z[j+i*n]
			}
		}
	}

	return ierr
}

// tqlrat - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:2102
//
//
//
//  Purpose:
//
//    TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This subroutine finds the eigenvalues of a symmetric
//    tridiagonal matrix by the rational QL method.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    Christian Reinsch,
//    Algorithm 464, TQLRAT,
//    Communications of the ACM,
//    Volume 16, page 689, 1973.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N].  On input, D contains the diagonal
//    elements of the matrix.  On output, D contains the eigenvalues in ascending
//    order.  If an error exit was made, then the eigenvalues are correct
//    in positions 1 through IERR-1, but may not be the smallest eigenvalues.
//
//    Input/output, double E2[N], contains in positions 2 through N
//    the squares of the subdiagonal elements of the matrix.  E2(1) is
//    arbitrary.  On output, E2 has been overwritten by workspace
//    information.
//
//    Output, int TQLRAT, error flag.
//    0, for no error,
//    J, if the J-th eigenvalue could not be determined after 30 iterations.
//
func tqlrat(n int, d []float64, e2 []float64) int {
	var b float64
	var c float64
	var f float64
	var g float64
	var h float64
	var i int
	var ierr int
	var ii int
	var j int
	var l int
	var l1 int
	var m int
	var mml int
	var p float64
	var r float64
	var s float64
	var t float64
	ierr = 0
	if n == 1 {
		return ierr
	}
	for i = 1; i < n; i++ {
		e2[i-1] = e2[i]
	}
	f = 0
	t = 0
	e2[n-1] = 0
	for l = 0; l < n; l++ {
		j = 0
		h = math.Abs(d[l]) + math.Sqrt(e2[l])
		if t <= h {
			t = h
			b = math.Abs(t) * r8_epsilon()
			c = b * b
		}

		//
		//  Look for small squared sub-diagonal element.
		//
		for m = l; m < n; m++ {
			if e2[m] <= c {
				break
			}
		}

		if m != l {
			for {
				if 30 <= j {
					ierr = l + 1
					return ierr
				}
				j = j + 1
				//
				//  Form shift.
				//
				l1 = l + 1
				s = math.Sqrt(e2[l])
				g = d[l]
				p = (d[l1] - g) / (2 * s)
				r = pythag(p, 1)
				d[l] = s / (p + math.Abs(r)*r8_sign(p))
				h = g - d[l]
				for i = l1; i < n; i++ {
					d[i] = d[i] - h
				}
				f = f + h
				//
				//  Rational QL transformation.
				//
				g = d[m]
				if g == 0 {
					g = b
				}
				h = g
				s = 0
				mml = m - l
				for ii = 1; ii <= mml; ii++ {
					i = m - ii
					p = g * h
					r = p + e2[i]
					e2[i+1] = s * r
					s = e2[i] / r
					d[i+1] = h + s*(h+d[i])
					g = d[i] - e2[i]/g
					if g == 0 {
						g = b
					}
					h = g * p / r
				}
				e2[l] = s * g
				d[l] = h
				if h == 0 {
					//
					//  Guard against underflow in convergence test.
					//
					break
				}
				if math.Abs(e2[l]) <= math.Abs(c/h) {
					break
				}
				e2[l] = h * e2[l]
				if e2[l] == 0 {
					break
				}
			}
		}
		p = d[l] + f

		//
		//  Order the eigenvalues.
		//
		for i = l; 0 <= i; i-- {
			if i == 0 {
				d[i] = p
				break
			} else if d[i-1] <= p {
				d[i] = p
				break
			}
			d[i] = d[i-1]
		}

	}
	return ierr
}

// tred1 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:2328
//
//
//
//  Purpose:
//
//    TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.
//
//  Discussion:
//
//    The routine reduces a real symmetric matrix to a symmetric
//    tridiagonal matrix using orthogonal similarity transformations.
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
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    Martin, Reinsch, Wilkinson,
//    TRED1,
//    Numerische Mathematik,
//    Volume 11, pages 181-195, 1968.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input/output, double A[N*N], on input, contains the real
//    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
//    On output, A contains information about the orthogonal transformations
//    used in the reduction in its strict lower triangle.
//    The full upper triangle of A is unaltered.
//
//    Output, double D[N], contains the diagonal elements of the
//    tridiagonal matrix.
//
//    Output, double E[N], contains the subdiagonal elements of the
//    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.
//
//    Output, double E2[N], contains the squares of the corresponding
//    elements of E.  E2 may coincide with E if the squares are not needed.
//
func tred1(n int, a []float64, d []float64, e []float64, e2 []float64) {
	var f float64
	var g float64
	var h float64
	var i int
	// var ii int
	var j int
	var k int
	var l int
	var scale float64

	for j = 0; j < n; j++ {
		d[j] = a[n-1+j*n]
	}

	for i = 0; i < n; i++ {
		a[n-1+i*n] = a[i+i*n]
	}
	for i = n - 1; 0 <= i; i-- {
		l = i - 1
		h = 0
		//
		//  Scale row.
		//
		scale = 0
		for k = 0; k <= l; k++ {
			scale = scale + math.Abs(d[k])
		}
		if scale == 0 {
			for j = 0; j <= l; j++ {
				d[j] = a[l+j*n]
				a[l+j*n] = a[i+j*n]
				a[i+j*n] = 0
			}
			e[i] = 0
			e2[i] = 0
			continue
		}
		for k = 0; k <= l; k++ {
			d[k] = d[k] / scale
		}
		for k = 0; k <= l; k++ {
			h = h + d[k]*d[k]
		}
		e2[i] = h * scale * scale
		f = d[l]
		g = -math.Sqrt(h) * r8_sign(f)
		e[i] = scale * g
		h = h - f*g
		d[l] = f - g
		if 0 <= l {

			//
			//  Form A * U.
			//
			for k = 0; k <= l; k++ {
				e[k] = 0
			}

			for j = 0; j <= l; j++ {
				f = d[j]
				g = e[j] + a[j+j*n]*f
				for k = j + 1; k <= l; k++ {
					g = g + a[k+j*n]*d[k]
					e[k] = e[k] + a[k+j*n]*f
				}
				e[j] = g
			}
			//
			//  Form P.
			//
			f = 0
			for j = 0; j <= l; j++ {
				e[j] = e[j] / h
				f = f + e[j]*d[j]
			}
			h = f / (h + h)

			//
			//  Form Q.
			//
			for j = 0; j <= l; j++ {
				e[j] = e[j] - h*d[j]
			}

			//
			//  Form reduced A.
			//
			for j = 0; j <= l; j++ {
				f = d[j]
				g = e[j]
				for k = j; k <= l; k++ {
					a[k+j*n] = a[k+j*n] - f*e[k] - g*d[k]
				}
			}

		}
		for j = 0; j <= l; j++ {
			f = d[j]
			d[j] = a[l+j*n]
			a[l+j*n] = a[i+j*n]
			a[i+j*n] = f * scale
		}
	}
}

// tred2 - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:2527
//
//
//
//  Purpose:
//
//    TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
//
//  Discussion:
//
//    This subroutine reduces a real symmetric matrix to a
//    symmetric tridiagonal matrix using and accumulating
//    orthogonal similarity transformations.
//
//    A and Z may coincide, in which case a single storage area is used
//    for the input of A and the output of Z.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    Martin, Reinsch, Wilkinson,
//    TRED2,
//    Numerische Mathematik,
//    Volume 11, pages 181-195, 1968.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the real symmetric input matrix.  Only the
//    lower triangle of the matrix need be supplied.
//
//    Output, double D[N], the diagonal elements of the tridiagonal
//    matrix.
//
//    Output, double E[N], contains the subdiagonal elements of the
//    tridiagonal matrix in E(2:N).  E(1) is set to zero.
//
//    Output, double Z[N*N], the orthogonal transformation matrix
//    produced in the reduction.
//
func tred2(n int, a []float64, d []float64, e []float64, z []float64) {

	for j := 0; j < n; j++ {
		for i := j; i < n; i++ {
			z[i+j*n] = a[i+j*n]
		}
	}

	for j := 0; j < n; j++ {
		d[j] = a[n-1+j*n]
	}

	for i := n - 1; 1 <= i; i-- {
		l := i - 1
		h := 0.0
		//
		//  Scale row.
		//
		scale := 0.0
		for k := 0; k <= l; k++ {
			scale = scale + math.Abs(d[k])
		}
		if scale == 0 {
			e[i] = d[l]
			for j := 0; j <= l; j++ {
				d[j] = z[l+j*n]
				z[i+j*n] = 0
				z[j+i*n] = 0
			}
			d[i] = 0
			continue
		}
		for k := 0; k <= l; k++ {
			d[k] = d[k] / scale
		}
		h = 0
		for k := 0; k <= l; k++ {
			h = h + d[k]*d[k]
		}
		f := d[l]
		g := -math.Sqrt(h) * r8_sign(f)
		e[i] = scale * g
		h = h - f*g
		d[l] = f - g

		//
		//  Form A*U.
		//
		for k := 0; k <= l; k++ {
			e[k] = 0
		}

		for j := 0; j <= l; j++ {
			f = d[j]
			z[j+i*n] = f
			g = e[j] + z[j+j*n]*f
			for k := j + 1; k <= l; k++ {
				g = g + z[k+j*n]*d[k]
				e[k] = e[k] + z[k+j*n]*f
			}
			e[j] = g
		}

		//
		//  Form P.
		//
		for k := 0; k <= l; k++ {
			e[k] = e[k] / h
		}

		f = 0
		for k := 0; k <= l; k++ {
			f = f + e[k]*d[k]
		}
		hh := 0.5 * f / h

		//
		//  Form Q.
		//
		for k := 0; k <= l; k++ {
			e[k] = e[k] - hh*d[k]
		}

		//
		//  Form reduced A.
		//
		for j := 0; j <= l; j++ {
			f = d[j]
			g = e[j]
			for k := j; k <= l; k++ {
				z[k+j*n] = z[k+j*n] - f*e[k] - g*d[k]
			}
			d[j] = z[l+j*n]
			z[i+j*n] = 0
		}

		d[i] = h
	}

	//
	//  Accumulation of transformation matrices.
	//
	for i := 1; i < n; i++ {
		l := i - 1
		z[n-1+l*n] = z[l+l*n]
		z[l+l*n] = 1
		h := d[i]
		if h != 0 {
			for k := 0; k <= l; k++ {
				d[k] = z[k+i*n] / h
			}
			for j := 0; j <= l; j++ {
				g := 0.0
				for k := 0; k <= l; k++ {
					g = g + z[k+i*n]*z[k+j*n]
				}
				for k := 0; k <= l; k++ {
					z[k+j*n] = z[k+j*n] - g*d[k]
				}
			}
		}
		for k := 0; k <= l; k++ {
			z[k+i*n] = 0
		}
	}

	for j := 0; j < n; j++ {
		d[j] = z[n-1+j*n]
	}

	for j := 0; j < n-1; j++ {
		z[n-1+j*n] = 0
	}

	z[n-1+(n-1)*n] = 1
	e[0] = 0
}
