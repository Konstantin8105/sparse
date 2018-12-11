package sparse

import (
	"math"
	"testing"
)

func TestPM(t *testing.T) {
	// tolerance
	eps := 1e-7

	t.Run("2x2", func(t *testing.T) {
		T, err := NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			Entry(T, 0, 0, 2),
			Entry(T, 0, 1, -12),
			Entry(T, 1, 0, 1),
			Entry(T, 1, 1, -5),
		}
		for i := range errs {
			if errs[i] != nil {
				panic(errs[i])
			}
		}

		// compress
		A, err := Compress(T)
		if err != nil {
			panic(err)
		}

		_, x, err := PM(A, &PmConfig{
			IterationMax: 1000,
			Tolerance:    1e-8,
		})
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
		T, err := NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			Entry(T, 0, 0, 1),
			Entry(T, 0, 1, 2),

			Entry(T, 1, 0, -2),
			Entry(T, 1, 1, 1),
			Entry(T, 1, 2, 2),

			Entry(T, 2, 0, 1),
			Entry(T, 2, 1, 3),
			Entry(T, 2, 2, 1),
		}
		for i := range errs {
			if errs[i] != nil {
				panic(errs[i])
			}
		}

		// compress
		A, err := Compress(T)
		if err != nil {
			panic(err)
		}

		_, x, err := PM(A, &PmConfig{
			IterationMax: 1000,
			Tolerance:    1e-8,
		})
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

	t.Run("2x2: IterationMax error", func(t *testing.T) {
		// too small amount iterations

		T, err := NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			Entry(T, 0, 0, 2),
			Entry(T, 0, 1, -12),
			Entry(T, 1, 0, 1),
			Entry(T, 1, 1, -5),
		}
		for i := range errs {
			if errs[i] != nil {
				panic(errs[i])
			}
		}

		// compress
		A, err := Compress(T)
		if err != nil {
			panic(err)
		}

		_, _, err = PM(A, &PmConfig{IterationMax: 1})
		if err == nil {
			t.Errorf("cannot check max iter")
		}
	})
	t.Run("oneMax: error", func(t *testing.T) {
		defer func() {
			if r := recover(); r != nil {
				t.Fatalf("Not valid")
			}
		}()
		oneMax(nil)
	})
}