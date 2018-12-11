package sparse_test

import (
	"math"
	"testing"

	"github.com/Konstantin8105/sparse"
)

func TestPM(t *testing.T) {
	// tolerance
	eps := 1e-7

	t.Run("2x2", func(t *testing.T) {
		T, err := sparse.NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			sparse.Entry(T, 0, 0, 2),
			sparse.Entry(T, 0, 1, -12),
			sparse.Entry(T, 1, 0, 1),
			sparse.Entry(T, 1, 1, -5),
		}
		for i := range errs {
			if errs[i] != nil {
				panic(errs[i])
			}
		}

		// compress
		A, err := sparse.Compress(T)
		if err != nil {
			panic(err)
		}

		_, x, err := sparse.PM(A)
		if err != nil {
			panic(err)
		}

		// result checking
		xExpect := []float64{1.0, 1.0 / 3.0}
		for i := range x {
			if math.Abs(x[i]-xExpect[i]) > eps {
				t.Errorf("Not correct %d : %e != %e", i, x[i], xExpect[i])
			}
		}
	})
	t.Run("3x3", func(t *testing.T) {
		T, err := sparse.NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			sparse.Entry(T, 0, 0, 1),
			sparse.Entry(T, 0, 1, 2),

			sparse.Entry(T, 1, 0, -2),
			sparse.Entry(T, 1, 1, 1),
			sparse.Entry(T, 1, 2, 2),

			sparse.Entry(T, 2, 0, 1),
			sparse.Entry(T, 2, 1, 3),
			sparse.Entry(T, 2, 2, 1),
		}
		for i := range errs {
			if errs[i] != nil {
				panic(errs[i])
			}
		}

		// compress
		A, err := sparse.Compress(T)
		if err != nil {
			panic(err)
		}

		_, x, err := sparse.PM(A)
		if err != nil {
			panic(err)
		}

		// result checking
		xExpect := []float64{1.0 / 2.0, 1.0 / 2.0, 1.0}
		for i := range x {
			if math.Abs(x[i]-xExpect[i]) > eps {
				t.Errorf("Not correct %d : %e != %e", i, x[i], xExpect[i])
			}
		}
	})
}
