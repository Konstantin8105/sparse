package sparse

import (
	"fmt"

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
func (lu *LU) Factorize(A *Matrix) error {
	// check input data
	et := errors.New("Function LU.Factorize: check input data")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}

	if et.IsError() {
		return et
	}

	lu.a = A

	// partial pivoting tolerance
	tol := func() float64 {
		if is_sym(lu.a) == 1 {
			return 0.001
		}
		return 1
	}()
	// ordering and symbolic analysis
	lu.s = cs_sqr(lu.order, lu.a, false)
	// numeric LU factorization
	lu.n = cs_lu(lu.a, lu.s, tol)

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
