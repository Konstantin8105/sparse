package sparse

import (
	"fmt"
	"sort"

	"github.com/Konstantin8105/errors"
)

// is_sym - 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
func is_sym(A *Matrix) int {
	n, m := A.n, A.m
	Ap, Ai := A.p, A.i
	if m != n {
		return 0
	}
	is_upper := true
	is_lower := true
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			if Ai[p] > j {
				is_upper = false
			}
			if Ai[p] < j {
				is_lower = false
			}
		}
	}
	if is_upper {
		return 1
	}
	if is_lower {
		return -1
	}
	return 0
}

// LU is a type for creating and using the LU factorization of a matrix.
type LU struct {
	order Order
	a     *Matrix
	s     *css
	n     *csn
}

func (lu *LU) Order(order Order) {
	lu.order = order
}

// Factorize computes the LU factorization of the matrix a and stores
// the result. Input matrix A is not checked on singular error.
// List `ignore` is list of ignore row and column in calculation.
func (lu *LU) Factorize(A *Matrix, ignore ...int) error {
	// check input data
	et := errors.New("Function LU.Factorize: check input data")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
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
		return et
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
		var pnew []int // append([]int{}, 0)
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

	// store
	lu.a = C

	// partial pivoting tolerance
	tol := func() float64 {
		if is_sym(lu.a) == 1 {
			return 0.001
		}
		return 1
	}()
	// ordering and symbolic analysis
	lu.s = cs_sqr(lu.order, lu.a, false)
	if lu.s == nil {
		return fmt.Errorf("matrix S in LU decomposition is nil")
	}
	// numeric LU factorization
	lu.n = cs_lu(lu.a, lu.s, tol)
	if lu.n == nil {
		return fmt.Errorf("matrix N in LU decomposition is nil")
	}

	return nil
}

// Solve solves a system of linear equations using the LU decomposition of a matrix
//
//	A * x = b
//
func (lu *LU) Solve(b []float64) (x []float64, _ error) {
	// check input data
	et := errors.New("Function LU.Solve: check input data")
	if b == nil {
		_ = et.Add(fmt.Errorf("vector b is nil"))
	}
	if lu.s == nil {
		_ = et.Add(fmt.Errorf("matrix S in LU decomposition is nil"))
	}
	if lu.n == nil {
		_ = et.Add(fmt.Errorf("matrix N in LU decomposition is nil"))
	}

	if et.IsError() {
		return nil, et
	}

	// initialization
	n := lu.a.n

	// get workspace
	x = make([]float64, n)

	bCopy := make([]float64, n)
	for i := range b {
		bCopy[i] = b[i]
	}
	defer func() {
		for i := range b {
			b[i] = bCopy[i]
		}
	}()

	// x = b(p)
	cs_ipvec(lu.n.pinv, b, x, n)
	// x = L\x
	cs_lsolve(lu.n.L, x)
	// x = U\x
	cs_usolve(lu.n.U, x)
	// b(q) = x
	cs_ipvec(lu.s.q, x, b, n)

	return x, nil
}
