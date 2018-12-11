package sparse

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/errors"
)

// PmConfig is default property of power method
type PmConfig struct {
	// Maximal amount iteration
	IterationMax uint64

	// Iteration tolerance
	Tolerance float64
}

// zeroize - set 0.0 in each element of slice
func zeroize(x []float64) {
	for i := range x {
		x[i] = 0.0
	}
}

// PM is power method for approximating eigenvalues.
// Find `dominant eigenvalue` of matrix A.
// Vector `x` is eigenvector.
// Value `𝛌` is eigenvalue.
func PM(A *Matrix, config *PmConfig) (𝛌 float64, x []float64, err error) {
	// check input data
	et := errors.New("Function LU.Factorize: check input data")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if A != nil && (A.n < 1 || A.m < 1) {
		_ = et.Add(fmt.Errorf("matrix A is small"))
	}

	if et.IsError() {
		err = et
		return
	}

	// minimal configuration
	if config == nil {
		config = &PmConfig{
			IterationMax: 500,
			Tolerance:    1e-5,
		}
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// workspace
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

	𝛌, 𝛌Last := -1.0, -1.0
	var iter uint64

	// calculation
	for {
		// x(k) = A*x(k-1)
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

		if math.Abs(delta) < config.Tolerance {
			break
		}
		if iter >= config.IterationMax {

			// TODO : test to iterations

			err = fmt.Errorf("Max iteration breaking: %d >= %d. Tolerance: %.5e",
				iter, config.IterationMax, delta)
			return
		}

		𝛌Last = 𝛌
		iter++
	}

	return
}
