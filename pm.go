package sparse

import (
	"fmt"
	"math"
)

// PmConfig is default property of power method
var PmConfig = struct {
	// Maximal amount iteration
	IterationMax uint64

	// Iteration tolerance
	Tolerance float64
}{
	IterationMax: 50000,
	Tolerance:    1e-8,
}

// PM is power method for approximating eigenvalues.
// Find `dominant eigenvalue` of matrix A.
// Vector `x` is eigenvector.
// Value `𝛌` is eigenvalue.
func PM(A *Matrix) (𝛌 float64, x []float64, err error) {
	if A.n < 1 || A.m < 1 {
		panic("1")
	}
	if A == nil {
		panic("2")
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// main operation: x(k) = A*x(k-1)
	x = make([]float64, A.m)
	xNext := make([]float64, A.m)
	x[0] = 1

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
	harmonize(x)

	zeroize := func(x []float64) {
		for i := range x {
			x[i] = 0.0
		}
	}

	𝛌, 𝛌Last := -1.0, -1.0
	var iter uint64

	// calculation
	for {
		// multiplication
		for j := 0; j < n; j++ {
			for p := Ap[j]; p < Ap[j+1]; p++ {
				xNext[Ai[p]] += Ax[p] * x[j]
			}
		}
		x, xNext = xNext, x

		zeroize(xNext)
		harmonize(x)

		// compute the Rayleigh quotient
		// 𝛌 = (Ax · x) / (x · x)

		// Ax
		zeroize(xNext)
		for j := 0; j < n; j++ {
			for p := Ap[j]; p < Ap[j+1]; p++ {
				xNext[Ai[p]] += Ax[p] * x[j]
			}
		}

		// up : Ax · x
		var up float64
		for i := range x {
			// TODO : check overflow
			up += x[i] * xNext[i]
		}

		// down : x · x
		var down float64
		for i := range x {
			// TODO : check overflow
			down += x[i] * x[i]
		}

		//
		𝛌 = up / down

		// TODO : if down is zero

		// TODO : if up is Nan

		// TODO : if lambda is Nan

		delta := math.Abs((math.Abs(𝛌) - math.Abs(𝛌Last)) / 𝛌)

		if math.Abs(delta) < PmConfig.Tolerance {
			break
		}
		if iter >= PmConfig.IterationMax {

			// TODO : test to iterations

			err = fmt.Errorf("Max iteration breaking: %d >= %d. Tolerance: %.5e",
				iter, PmConfig.IterationMax, delta)
			return
		}

		𝛌Last = 𝛌
		iter++
	}

	return
}
