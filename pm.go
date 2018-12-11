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

// PM power method result
type PM struct {
	a *Matrix
	E []*Eigen
}

// Eigen result
type Eigen struct {
	// eigenvalue
	𝜦 float64

	// eigenvector
	𝑿 []float64
}

// zeroize - set 0.0 in each element of slice
func zeroize(x []float64) {
	for i := range x {
		x[i] = 0.0
	}
}

// oneMax - modify slice with max value 1.0
func oneMax(x []float64) {
	// check input data
	if len(x) == 0 {
		return
	}

	// find max value
	max := x[0]
	for i := range x {
		if math.Abs(x[i]) > math.Abs(max) {
			max = x[i]
		}
	}

	// modification
	for i := range x {
		x[i] /= max
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

	// workspace
	x = make([]float64, A.m)
	xNext := make([]float64, A.m)

	// initial values
	x[0] = 1
	oneMax(x)
	𝛌, 𝛌Last := -1.0, -1.0

	var iter uint64

	// calculation
	for {
		// x(k) = A*x(k-1)
		if err = Gaxpy(A, x, xNext); err != nil {
			return
		}
		x, xNext = xNext, x

		zeroize(xNext)
		oneMax(x)

		// compute the Rayleigh quotient
		// 𝛌 = (Ax · x) / (x · x)

		// Ax
		zeroize(xNext)
		if err = Gaxpy(A, x, xNext); err != nil {
			return
		}

		// up : Ax · x
		var up float64
		for i := range x {
			up += x[i] * xNext[i]
		}

		// down : x · x
		var down float64
		for i := range x {
			down += x[i] * x[i]
		}

		if math.Abs(down) == 0.0 || math.IsNaN(up) {
			err = fmt.Errorf("Not acceptable value")
			return
		}

		// calculation eigenvalue
		𝛌 = up / down

		// check breaking
		delta := math.Abs((math.Abs(𝛌) - math.Abs(𝛌Last)) / 𝛌)

		if math.Abs(delta) < config.Tolerance {
			break
		}
		if iter >= config.IterationMax {
			err = ErrorPm{
				Iteration:    iter,
				IterationMax: config.IterationMax,
				Delta:        delta,
				Tolerance:    config.Tolerance,
			}
			return
		}

		𝛌Last = 𝛌
		iter++
	}

	return
}

// ErrorPm error of power method
type ErrorPm struct {
	Iteration    uint64
	IterationMax uint64
	Delta        float64
	Tolerance    float64
}

func (e ErrorPm) Error() string {
	et := errors.New("Power method error")
	_ = et.Add(fmt.Errorf("Iteration: %d", e.Iteration))
	_ = et.Add(fmt.Errorf("Max.Iteration: %d", e.IterationMax))
	_ = et.Add(fmt.Errorf("Delta tolerance: %5e", e.Tolerance))
	_ = et.Add(fmt.Errorf("Tolerance: %5e", e.Tolerance))
	return et.Error()
}
