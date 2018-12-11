package sparse

import (
	"fmt"
	"math"
)

// PM is power method for approximating eigenvalues.
// Find `dominant eigenvalue` of matrix A.
// Vector `x` is eigenvector.
func PM(A *Matrix, x *[]float64) (err error) {
	if x == nil {
		panic("1")
	}
	if A == nil {
		panic("2")
	}
	if A.m != len(*x) {
		panic("3")
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// main operation: x(k) = A*x(k-1)
	xNext := make([]float64, A.m)

	harmonize := func(x []float64) {
		// max factor
		max := math.Abs(x[0])
		isMinus := math.Signbit(x[0])
		for i := range x {
			if math.Abs(x[i]) > max {
				max = math.Abs(x[i])
				isMinus = math.Signbit(x[0])
			}
		}
		for i := range x {
			if isMinus {
				x[i] /= -max
				continue
			}
			x[i] /= max
		}
	}
	harmonize(*x)

	zeroize := func(x []float64) {
		for i := range x {
			x[i] = 0.0
		}
	}

	𝛌, 𝛌Last := -1.0, -1.0
	eps := 1e-8
	var iter uint64
	var iterMax uint64 = 10000

	// calculation
	for iter = 0; iter < iterMax; iter++ {
		// multiplication
		for j := 0; j < n; j++ {
			for p := Ap[j]; p < Ap[j+1]; p++ {
				xNext[Ai[p]] += Ax[p] * (*x)[j]
			}
		}
		*x, xNext = xNext, *x

		zeroize(xNext)
		harmonize(*x)

		// compute the Rayleigh quotient
		// 𝛌 = (Ax · x) / (x · x)

		// Ax
		zeroize(xNext)
		for j := 0; j < n; j++ {
			for p := Ap[j]; p < Ap[j+1]; p++ {
				xNext[Ai[p]] += Ax[p] * (*x)[j]
			}
		}

		// up : Ax · x
		var up float64
		for i := range *x {
			up += (*x)[i] * xNext[i]
		}

		// down : x · x
		var down float64
		for i := range *x {
			down += (*x)[i] * (*x)[i]
		}

		//
		𝛌 = up / down

		delta := math.Abs(math.Abs(𝛌)-math.Abs(𝛌Last)) / math.Abs(𝛌)

		if math.Abs(delta) < eps {
			break
		}
		_ = eps
		𝛌Last = 𝛌
	}
	fmt.Println(">> 𝛌     : ", 𝛌)
	fmt.Println(">> iter  : ", iter)

	return nil
}
