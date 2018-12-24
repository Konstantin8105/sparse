package sparse

import (
	"math"
)

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
// removed

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

// timestamp - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/Eispack/eispack.c:1810
// removed

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
				r := math.Hypot(p, 1)
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
					r = math.Hypot(p, e[i])
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
	var b, c, f, g, h float64
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
				r = math.Hypot(p, 1)
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
