package sparse

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/errors"
)

// IsSingular return true if square matrix have zero in main diagonal
func IsSingular(A *Matrix) (bool, error) {
	// check input data
	et := errors.New("Function IsSingular: check input data")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if A != nil {
		if A.n != A.m {
			_ = et.Add(fmt.Errorf("matrix A is not square: %d != %d", A.n, A.m))
		}
	}

	if et.IsError() {
		return true, et
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// calculation
	for j := 0; j < n; j++ {
		if Ap[j+1]-Ap[j] == 0 {
			// zero column
			return true, nil
		}
		var found bool
		for p := Ap[j]; p < Ap[j+1]; p++ {
			i := Ai[p]
			x := Ax[p]
			if i == j {
				found = true
			}
			if i == j && x == 0 {
				return true, nil
			}
		}
		if !found {
			return true, nil
		}
	}

	return false, nil
}

// Zeroize - put zero in column and row with position `pos` and
// put value `center` in matrix position [pos,pos].
func Zeroize(A *Matrix, pos int, center float64) error {
	// check input data
	et := errors.New("Function Zeroize: check input data")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if A != nil {
		if pos < 0 {
			_ = et.Add(fmt.Errorf("position is less zero"))
		}
		if pos >= A.n {
			_ = et.Add(fmt.Errorf("position is more columns"))
		}
		if pos >= A.m {
			_ = et.Add(fmt.Errorf("position is more rows"))
		}
	}
	if math.IsNaN(center) {
		_ = et.Add(fmt.Errorf("center value is Nan value"))
	}
	if math.IsInf(center, 0) {
		_ = et.Add(fmt.Errorf("center value is infinity value"))
	}

	if et.IsError() {
		return et
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// calculation
	var found bool
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			i := Ai[p] // row
			j := j     // column
			if i == pos || j == pos {
				Ax[p] = 0.0
			}
			if i == pos && j == pos {
				found = true
				Ax[p] = center
			}
		}
	}
	if found {
		return Dupl(A)
	}

	panic("Not found")

	// TODO: not found

	return nil
}
