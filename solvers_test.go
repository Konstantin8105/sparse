package sparse_test

import (
	"fmt"
	"os"
	"testing"

	"github.com/Konstantin8105/sparse"
)

func ExampleLU() {
	// Solve next:
	// [ 1 5 ] [ x1 ] = [ 11 ]
	// [ 2 3 ] [ x2 ]   [  8 ]
	T, err := sparse.NewTriplet()
	if err != nil {
		panic(err)
	}
	// storage
	errs := []error{
		sparse.Entry(T, 0, 0, 1),
		sparse.Entry(T, 0, 1, 5),
		sparse.Entry(T, 1, 0, 2),
		sparse.Entry(T, 1, 1, 3),
	}
	for i := range errs {
		if errs[i] != nil {
			panic(errs[i])
		}
	}
	T.Print(false)

	// compress
	A, err := sparse.Compress(T)
	if err != nil {
		panic(err)
	}
	A.Print(false)

	// singinal check
	is, err := sparse.IsSingular(A)
	if err != nil {
		panic(err)
	}
	if is {
		panic("singular matrix")
	}

	// solving
	lu := new(sparse.LU)
	err = lu.Factorize(A)
	if err != nil {
		panic(err)
	}
	b := []float64{11, 8}

	x, err := lu.Solve(b)
	if err != nil {
		panic(err)
	}
	fmt.Fprintf(os.Stdout, "%v", x)

	// Output:
	// [1 2]
}

func TestLU(t *testing.T) {

}
