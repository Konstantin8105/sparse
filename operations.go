package sparse

import (
	"fmt"

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
		for p := Ap[j]; p < Ap[j+1]; p++ {
			i := Ai[p]
			x := Ax[p]
			if i == j && x == 0 {
				return true, nil
			}
		}
	}

	return false, nil
}

// IsSymmetrical return true is square matrix A is symmetrical
// func IsSymmetrical(A *Matrix) (bool, error) {
// 	// check input data
// 	et := errors.New("Function IsSingular: check input data")
// 	if A == nil {
// 		_ = et.Add(fmt.Errorf("matrix A is nil"))
// 	}
// 	if A != nil && A.nz != -1 {
// 		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
// 	}
// 	if A != nil {
// 		if A.n != A.m {
// 			_ = et.Add(fmt.Errorf("matrix A is not square: %d != %d", A.n, A.m))
// 		}
// 	}
//
// 	if et.IsError() {
// 		return true, et
// 	}
//
// 	// initialization
// 	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x
//
// 	// calculation
// 	for j := 0; j < n; j++ {
// 		if Ap[j+1]-Ap[j] == 0 {
// 			// zero column
// 			return true, nil
// 		}
// 		for p := Ap[j]; p < Ap[j+1]; p++ {
// 			i := Ai[p]
// 			x := Ax[p]
// 			if i == j && x == 0 {
// 				return true, nil
// 			}
// 		}
// 	}
//
// 	return false, nil
// }
