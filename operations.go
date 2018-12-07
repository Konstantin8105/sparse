package sparse

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/errors"
)

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

// IsSym return true if matrix A is symmetrical
func IsSym(A *Matrix) (ok bool, err error) {
	// check input data
	et := errors.New("Function Limits: check input data")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if A != nil {
		if A.n == 0 && A.m == 0 {
			_ = et.Add(fmt.Errorf("matrix A is empty"))
		}
		if A.n != A.m {
			_ = et.Add(fmt.Errorf("matrix A is not square"))
		}
	}

	if et.IsError() {
		return false, et
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// check amount upper and lower entries
	{
		nzUpper := 0
		nzLower := 0
		for j := 0; j < n; j++ {
			for p := Ap[j]; p < Ap[j+1]; p++ {
				i := Ai[p] // row
				j := j     // column
				if i == j {
					// ignore diagonal
					continue
				}
				if i > j { // lower
					nzLower++
					continue
				}
				// upper
				nzUpper++
			}
		}
		if nzUpper != nzLower {
			return false, fmt.Errorf(
				"amount entries on up and low diagonal is not same: %d != %d",
				nzUpper, nzLower)
		}
	}

	defer func() {
		// unmark all
		for j := 0; j < n; j++ {
			for p := Ap[j]; p < Ap[j+1]; p++ {
				Ai[p] = unflip(Ai[p])
			}
		}
	}()

	// check values
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			if marked(Ai, p) {
				// ignore marked
				continue
			}

			i := Ai[p] // row
			j := j     // column
			x := Ax[p] // value

			if i == j {
				continue
			}

			// finding
			var found bool
			for p2 := Ap[i]; p2 < Ap[i+1]; p2++ {
				i2 := Ai[p2] // row
				j2 := i      // column
				x2 := Ax[p2] // value
				if marked(Ai, p2) || i != j2 || j != i2 || i2 == j2 {
					// ignore marked
					// or
					// coordinates of entries is not same
					// or
					// ignore diagonal entries
					continue
				}
				found = true
				if x != x2 {
					return false, fmt.Errorf(
						"matrix is not symmetric. Value[%d,%d] = %.5e. Value[%d,%d] = %5e.",
						i, j, x,
						i2, j2, x2)
				}
				if found {
					// mark
					mark(Ai, p2)
					mark(Ai, p)

					break
				}
			}
			if !found {
				return false, fmt.Errorf(
					"matrix is not symmetric. Cannot find entry [%d,%d]",
					i, j)
			}
		}
	}

	return true, nil
}
