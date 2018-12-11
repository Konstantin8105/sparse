package sparse_test

import (
	"math"
	"testing"

	"github.com/Konstantin8105/sparse"
)

func TestPM(t *testing.T) {
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

		x := make([]float64, 2)
		x[0] = 10
		x[1] = 1

		_, x, err = sparse.PM(A)
		if err != nil {
			panic(err)
		}

		// result checking
		eps := 1e-7
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

		x := make([]float64, 3)
		x[0] = 1
		x[1] = 1
		x[2] = 0

		_, x, err = sparse.PM(A)
		if err != nil {
			panic(err)
		}

		// result checking
		eps := 1e-7
		xExpect := []float64{1.0 / 2.0, 1.0 / 2.0, 1.0}
		for i := range x {
			if math.Abs(x[i]-xExpect[i]) > eps {
				t.Errorf("Not correct %d : %e != %e", i, x[i], xExpect[i])
			}
		}
	})
}
