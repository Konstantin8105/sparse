package sparse

import (
	"fmt"
	"math"
	"sort"

	"github.com/Konstantin8105/errors"
)

// PmConfig is default property of power method
type PmConfig struct {
	// Maximal amount iteration
	IterationMax uint64

	// Iteration tolerance
	Tolerance float64
}

// PM power method for approximating eigenvalues.
type PM struct {
	a      *Matrix
	ignore []int
	config PmConfig

	// result of calculation
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
// List `ignore` is list of ignore row and column in calculation.
func (pm *PM) Factorize(A *Matrix, config *PmConfig, ignore ...int) (err error) {
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

	// sort ignore list
	(sort.IntSlice(ignore)).Sort()
	if len(ignore) > 0 {
		if ignore[0] < 0 {
			_ = et.Add(fmt.Errorf("ignore list have index less zero: %d", ignore[0]))
		}
		if ignore[len(ignore)-1] >= A.m || ignore[len(ignore)-1] >= A.n {
			_ = et.Add(fmt.Errorf("ignore list have index outside matrix"))
		}
	}

	if et.IsError() {
		err = et
		return
	}

	// remove duplicate from `ignore` list
	if len(ignore) > 0 {
		list := append([]int{}, ignore[0])
		for i := 1; i < len(ignore); i++ {
			if ignore[i-1] == ignore[i] {
				continue
			}
			list = append(list, ignore[i])
		}
		list, ignore = ignore, list
		cs_free(list)
	}

	// coping matrix A without ignored rows and columns
	C, _ := A.Copy()
	if len(ignore) > 0 {
		_, err := Fkeep(C, func(i, j int, x float64) bool {
			var found bool
			for k := range ignore {
				if i == ignore[k] || j == ignore[k] {
					found = true
					break
				}
			}
			return !found // if not found , then keep value
		})
		if err != nil {
			return err
		}

		// remove empty column
		var pnew []int
		for j := 0; j < C.n; j++ {
			var found bool
			for k := range ignore {
				if j == ignore[k] {
					found = true
					break
				}
			}
			if found {
				continue
			}
			pnew = append(pnew, C.p[j])
		}
		C.p = append(pnew, C.p[C.n])

		// recalculate amount rows and columns
		C.n -= len(ignore)
		C.m -= len(ignore)

		// recalculate vector i
		for j := 0; j < C.n; j++ {
			for p := C.p[j]; p < C.p[j+1]; p++ {
				incr := 0
				for k := range ignore {
					if C.i[p] > ignore[k] {
						incr++
					}
				}
				C.i[p] -= incr
			}
		}
	}

	// minimal configuration
	if config == nil {
		config = &PmConfig{
			IterationMax: 500,
			Tolerance:    1e-5,
		}
	}

	// store
	pm.ignore = ignore
	pm.a = C
	pm.config = *config

	return
}

// Next calculate next `amount` eigenvalues
func (pm *PM) Next(amount int) (err error) {
	switch {
	case amount < 0:
		return fmt.Errorf("Not valid amount: %d < 0", amount)
	case amount == 0:
		return nil
	case amount > 1:
		for i := 0; i < amount; i++ {
			err = pm.Next(1)
			if err != nil {
				return
			}
		}
		return nil
	}
	// Now `amount = 1`

	// workspace
	var (
		x     = make([]float64, pm.a.m)
		xNext = make([]float64, pm.a.m)
	)

	// initial values
	x[0] = 1
	oneMax(x)
	𝛌, 𝛌Last := -1.0, -1.0

	var iter uint64

	// calculation
	for {
		// x(k) = A*x(k-1)
		if err = Gaxpy(pm.a, x, xNext); err != nil {
			return
		}
		x, xNext = xNext, x

		zeroize(xNext)
		oneMax(x)

		// compute the Rayleigh quotient
		// 𝛌 = (Ax · x) / (x · x)

		// Ax
		zeroize(xNext)
		if err = Gaxpy(pm.a, x, xNext); err != nil {
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

		if math.Abs(delta) < pm.config.Tolerance {
			break
		}
		if iter >= pm.config.IterationMax {
			err = ErrorPm{
				Iteration:    iter,
				IterationMax: pm.config.IterationMax,
				Delta:        delta,
				Tolerance:    pm.config.Tolerance,
			}
			return
		}

		𝛌Last = 𝛌
		iter++
	}
	// TODO (KI) : add ignore list

	pm.E = append(pm.E, &Eigen{
		// eigenvalue
		𝜦: 𝛌,

		// eigenvector
		𝑿: x,
	})

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
