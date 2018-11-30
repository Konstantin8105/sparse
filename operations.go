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
		if pos >= A.n {
			_ = et.Add(fmt.Errorf("position is more columns"))
		}
		if pos >= A.m {
			_ = et.Add(fmt.Errorf("position is more rows"))
		}
	}
	if pos < 0 {
		_ = et.Add(fmt.Errorf("position is less zero"))
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
	var nz int = 0
	for j := 0; j < n; j++ {
		// get current location of col j
		p := Ap[j]
		// record new location of col j
		Ap[j] = nz
		for ; p < Ap[j+1]; p++ {
			i := Ai[p] // row
			j := j     // column
			if (i != pos && j != pos) || (i == pos && j == pos) {
				if Ax != nil {
					// keep A(i,j)
					Ax[nz] = Ax[p]
					if i == pos && j == pos {
						found = true
						Ax[nz] = center
					}
				}
				Ai[nz] = Ai[p]
				nz++
			}
		}
	}
	// finalize A
	Ap[n] = nz
	// remove extra space from A
	cs_sprealloc(A, 0)
	if found {
		return nil
	}

	return fmt.Errorf("cannot found center entry")
}
