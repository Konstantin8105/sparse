package sparse

import (
	"fmt"
	"math"
	"math/rand"
	"sort"

	"github.com/Konstantin8105/errors"
)

// PmConfig is default property of power method
type PmConfig struct {
	// Maximal amount iteration
	IterationMax uint64

	// Iteration tolerance
	// TODO (KI) : add comment about for each lambda precasion is lost
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
	min := x[0]
	for i := range x {
		if min < x[i] && x[i] < max {
			continue
		}
		if x[i] < min {
			min = x[i]
		}
		if x[i] > max {
			max = x[i]
		}
	}
	if math.Abs(min) > max {
		max = min
	}
	// if max == 1.0 {
	// 	return
	// }
	//
	// if max == 66.66 , then max = 100
	// if max == 0.006 , then max = 0.01
	// if max > 0 {
	// 	max = math.Pow(10.0, float64(int(math.Log10(max))))
	// } else {
	// 	max = -math.Pow(10.0, float64(int(math.Log10(-max))))
	// }

	// modification
	for i := range x {
		// TODO (KI) : need parallel
		x[i] /= max
	}
}

// PM is power method for approximating eigenvalues.
// Find `dominant eigenvalue` of matrix A.
// List `ignore` is list of ignore row and column in calculation.
//
//	Algorith : Power Method
//	x(0) // initial vector
//	k = 1
//	for δ < ɛ {
//		x(k) = A · x(k-1)
//		δ = || x(k) - x(k-1) ||1
//		k = k + 1
//	}
//
// See next articles:
//
//	1. Sepandar D. Kamvar, Taher H. Haveliwala, Christopher D. Manning, Gene H. Golub
//	"Extrapolation Methods for Accelerating PageRank Computations"
//
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

	// matrix A must have diagonal element
	for j := 0; j < pm.a.n; j++ {
		found := false
		for p := pm.a.p[j]; p < pm.a.p[j+1]; p++ {
			if j != pm.a.i[p] {
				continue
			}
			// only diagonal element
			found = true
		}
		if !found {
			// diagonal element is not found in that column
			pm.a.inject(j, j, 0.0)
		}
	}

	return
}

// TODO (KI) : add research PM from input values
// TODO (KI) : [7 4 1; 4 4 4; 1 4 7]
// TODO (KI) : xo = [1 2 3]T
// TODO (KI) : xo = [0 1 -1]T

// Next calculate next `amount` eigenvalues
func (pm *PM) Next(amount int) (err error) {
	switch {
	case amount < 0:
		return fmt.Errorf("Not valid amount: %d < 0", amount)
	case amount == 0:
		return nil
	case amount > 1:
		for i := 0; i < amount; i++ {
			if err = pm.Next(1); err != nil {
				return
			}
		}
		return nil
	}
	// Now `amount = 1`
	// defer func() {
	// 	if err == nil {
	// 		// TODO (KI): check that
	// 		if math.Abs(pm.E[len(pm.E)-1].𝜦) > math.Abs(pm.E[len(pm.E)-2].𝜦) {
	// 			err = fmt.Errorf("Any next eigenvalue must be less by absolute value:"+
	// 				"[%.3e,%.3e]", pm.E[len(pm.E)-1].𝜦, pm.E[len(pm.E)-2].𝜦)
	// 		}
	// 	}
	// }()

	// workspace
	var (
		x     = make([]float64, pm.a.m+len(pm.ignore))
		xNext = make([]float64, pm.a.m+len(pm.ignore))
	)
	x = x[:pm.a.m]
	xNext = xNext[:pm.a.m]

	// initial values
	// if len(pm.E) > 0 {
	// 	// suppose any next eigenvector is like last
	// 	copy(x, pm.E[len(pm.E)-1].𝑿)
	// 	// TODO (KI) : coping without ignore list
	// } else {
	for i := range x {
		x[i] = rand.Float64() - 0.5
	}
	// }
	dlast := 1.0
	oneMax(x)

	// iteration value
	var iter uint64

	// main iteration function
	iteration := func() {
		// x(k) = A*x(k-1) - (∑(𝛌(j)*[x(j)]*[x(j)Transpose])) * x(k-1),
		//        |      |   |                                       |
		//        +------+   +---------------------------------------+
		//          part1               part 2
		// where j = 0 ... k-1
		oneMax(x)
		zeroize(xNext)
		// fmt.Println("x = ", x)
		// fmt.Println("X--- ", x)
		// part 1
		// if err = Gaxpy(pm.a, x, xNext); err != nil {
		// 	return
		// }
		_, _ = Fkeep(pm.a, func(i, j int, val float64) bool {
			xNext[i] += val * x[j]
			return true
		})
		// part 2
		for i := range pm.E {
			// value pm.E.X without ignore elements
			EX := make([]float64, pm.a.m)
			if len(pm.ignore) == 0 { // ignore list is empty
				copy(EX, pm.E[i].𝑿)
			} else {
				// TODO (KI)
				panic("TODO")
			}
			for row := 0; row < pm.a.m; row++ {
				for col := 0; col < pm.a.m; col++ {
					xNext[row] -= pm.E[i].𝜦 * EX[row] * EX[col] * x[col]
				}
			}

			// debug
			// fmt.Println("------------")
			// T := make([]float64, pm.a.m*pm.a.m)
			// for row := 0; row < pm.a.m; row++ {
			// 	for col := 0; col < pm.a.m; col++ {
			// 		T[row*pm.a.m+col] -= pm.E[i].𝜦 * EX[row] * EX[col]
			// 	}
			// }
			// _, _ = Fkeep(pm.a, func(i, j int, x float64) bool {
			// 	T[i*pm.a.m+j] += x
			// 	fmt.Println("A[", i, ",", j, "]=", x)
			// 	return true
			// })
			//
			// fmt.Println("L  = ", pm.E[i].𝜦)
			// fmt.Println("EX = ", EX)
			// for row := 0; row < pm.a.m; row++ {
			// 	for col := 0; col < pm.a.m; col++ {
			// 		fmt.Printf("%.4e ", T[row*pm.a.m+col])
			// 	}
			// 	fmt.Println("")
			// }
			// fmt.Println("------------")

		}

		// fmt.Println("--> x ", xNext)
		// post work
		x, xNext = xNext, x
	}

	// calculation
	for {

		// fmt.Println("iter ", iter)

		// main iteration function
		oneMax(x)
		iteration()

		// first norm of vector
		// link: https://en.wikipedia.org/wiki/Norm_(mathematics)#Taxicab_norm_or_Manhattan_norm
		oneMax(x)
		oneMax(xNext)
		var d float64
		for i := range x {
			d += math.Abs(x[i] - xNext[i])
		}

		if math.Abs(dlast-d) < pm.config.Tolerance {
			// tolerance breaking
			break
		}
		if iter >= pm.config.IterationMax {
			err = ErrorPm{
				Iteration:    iter,
				IterationMax: pm.config.IterationMax,
				Delta:        dlast - d,
				Tolerance:    pm.config.Tolerance,
				err:          fmt.Errorf("iteration limit"),
			}
			return
		}

		dlast = d
		iter++
	}

	// compute the Rayleigh quotient
	// 𝛌 = (Ax · x) / (x · x)

	// calculation of Ax is ignore and takes value x(k-1)

	// up : Ax · x
	// zeroize(xNext)
	// if err = Gaxpy(pm.a, x, xNext); err != nil {
	// 	return
	// }
	iteration()

	// fmt.Println(">>>>>>>>>>>>>>. ", x, xNext)

	var up float64
	for i := range x {
		up += x[i] * xNext[i]
	}

	// down : x · x
	oneMax(x)
	var down float64
	for i := range x {
		down += x[i] * x[i]
	}

	if math.Abs(down) == 0.0 || math.IsNaN(up) {
		return fmt.Errorf("Not acceptable value")
	}

	// calculation eigenvalue
	𝛌 := up / down
	// fmt.Println("--> ", 𝛌)

	// change workspace with length of ignore list
	x = append(x, make([]float64, len(pm.ignore))...)

	// ignore list
	if len(pm.ignore) > 0 {
		// decomperess vector `x`
		// short x vector: x = [1 2]
		// ignore list   :     [0 2 4]
		// result        : x = [0 1 0 2 0]
		counter := 0
		for i := len(x) - 1; i >= 0; i-- {
			var found bool
			for k := range pm.ignore {
				if i == pm.ignore[k] {
					found = true
					break
				}
			}
			if found {
				x[i] = 0
				continue
			}
			counter++
			x[i] = x[len(x)-len(pm.ignore)-counter]
		}
	}

	oneMax(x)

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
	err          error
}

func (e ErrorPm) Error() string {
	et := errors.New("Power method error")
	_ = et.Add(fmt.Errorf("Iteration: %d", e.Iteration))
	_ = et.Add(fmt.Errorf("Max.Iteration: %d", e.IterationMax))
	_ = et.Add(fmt.Errorf("Delta tolerance: %5e", e.Tolerance))
	_ = et.Add(fmt.Errorf("Tolerance: %5e", e.Tolerance))
	_ = et.Add(e.err)
	return et.Error()
}
