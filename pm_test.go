package sparse

import (
	"bytes"
	"math"
	"testing"
)

func BenchmarkPM(b *testing.B) {
	T, err := NewTriplet()
	if err != nil {
		panic(err)
	}
	// storage
	errs := []error{
		Entry(T, 0, 3, 0.9090),
		Entry(T, 0, 4, 0.0910),

		Entry(T, 1, 0, 0.0002),
		Entry(T, 1, 2, 0.9998),

		Entry(T, 2, 1, 0.9998),
		Entry(T, 2, 3, 0.0002),

		Entry(T, 3, 0, 0.6690),
		Entry(T, 3, 4, 0.3310),

		Entry(T, 4, 0, 0.9989),
		Entry(T, 4, 1, 0.0011),
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

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		var pm PM
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 10000,
			Tolerance:    1e-6,
		})
		if err != nil {
			panic(err)
		}

		err = pm.Next(1)
		if err != nil {
			panic(err)
		}
	}
}

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
			Entry(T, 0, 0, 2000),

			Entry(T, 1, 1, 2),
			Entry(T, 1, 3, -12),

			Entry(T, 2, 2, 2000),

			Entry(T, 3, 1, 1),
			Entry(T, 3, 3, -5),

			Entry(T, 4, 4, 2000),
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

		var pm PM
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 1000,
			Tolerance:    1e-8,
		}, 0, 4, 2, 0, 2, 4)
		if err != nil {
			panic(err)
		}

		err = pm.Next(1)
		if err != nil {
			panic(err)
		}

		// result checking
		xExpect := []float64{1.0, 1.0 / 3.0}
		if math.Abs(pm.E[0].𝑿[1]-xExpect[0]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝑿[1], xExpect[0])
		}
		if math.Abs(pm.E[0].𝑿[3]-xExpect[1]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝑿[3], xExpect[1])
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

		var pm PM
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 1000,
			Tolerance:    1e-8,
		})
		if err != nil {
			panic(err)
		}

		err = pm.Next(2)
		if err != nil {
			panic(err)
		}

		// result checking
		xExpect := []float64{1.0 / 2.0, 1.0 / 2.0, 1.0}
		for i := range pm.E[0].𝑿 {
			if math.Abs(pm.E[0].𝑿[i]-xExpect[i]) > eps {
				t.Errorf("Not correct %d : %e != %e", i, pm.E[0].𝑿[i], xExpect[i])
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

		var pm PM
		err = pm.Factorize(A, &PmConfig{IterationMax: 1})
		if err != nil {
			panic(err)
		}

		err = pm.Next(1)
		if err == nil {
			t.Errorf("cannot check max iter")
		}
		t.Log(err)
		t.Log(err.Error())
	})
	t.Run("2x2: Iteration error", func(t *testing.T) {
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
		var pm PM
		_ = pm.Factorize(A, nil)
	})
	t.Run("oneMax: error", func(t *testing.T) {
		defer func() {
			if r := recover(); r != nil {
				t.Fatalf("Not valid")
			}
		}()
		oneMax(nil)
	})
	t.Run("A is nil", func(t *testing.T) {
		var pm PM
		err := pm.Factorize(nil, nil)
		if err == nil {
			t.Fatalf("not check")
		}
		t.Log(err)
	})
	t.Run("Triplet", func(t *testing.T) {
		var stdin bytes.Buffer
		stdin.WriteString("0 0 1\n 0 1 2\n 1 0 3\n 1 1 4")
		T, err := Load(&stdin)
		if err != nil {
			panic(err)
		}
		var pm PM
		err = pm.Factorize((*Matrix)(T), nil)
		if err == nil {
			t.Fatalf("not check")
		}
		t.Log(err)
	})
	t.Run("Small", func(t *testing.T) {
		var stdin bytes.Buffer
		stdin.WriteString("")
		T, err := Load(&stdin)
		if err != nil {
			panic(err)
		}
		A, err := Compress(T)
		if err != nil {
			panic(err)
		}
		var pm PM
		if err = pm.Factorize(A, nil, -1, 2, 100, 9); err == nil {
			t.Fatalf("not check")
		}
		if err = pm.Next(-1); err == nil {
			t.Fatalf("not check")
		}
		t.Log(err)
	})
}
