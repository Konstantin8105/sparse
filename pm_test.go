package sparse

import (
	"bytes"
	"fmt"
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

	for _, p := range []int{1, 3} {
		b.Run(fmt.Sprintf("%d", p), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				var pm PM
				err = pm.Factorize(A, &PmConfig{
					IterationMax: 1000000000,
					Tolerance:    1e-6,
				})
				if err != nil {
					panic(err)
				}

				err = pm.Next(p)
				if err != nil {
					panic(err)
				}
			}
		})
	}
}

func TestPM(t *testing.T) {
	// tolerance
	eps := 1e-7

	t.Run("specific2x2", func(t *testing.T) {
		T, err := NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			Entry(T, 0, 0, 13),
			Entry(T, 0, 1, 5),

			Entry(T, 1, 0, 2),
			Entry(T, 1, 1, 4),
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

		for i := range pm.E {
			t.Logf("pm.E = %d", i)
			t.Logf("𝑿 = %v", pm.E[i].𝑿)
			t.Logf("𝜦 = %v", pm.E[i].𝜦)
		}

		// result checking

		// lambda = 14
		xExpect := []float64{1.0, 0.2}
		if math.Abs(pm.E[0].𝑿[0]-xExpect[0]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝑿[0], xExpect[0])
		}
		if math.Abs(pm.E[0].𝑿[1]-xExpect[1]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝑿[1], xExpect[1])
		}

		// lambda = 3
		xExpect = []float64{0.55, 1.0}
		if math.Abs(pm.E[1].𝑿[0]-xExpect[0]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[1].𝑿[0], xExpect[0])
		}
		if math.Abs(pm.E[1].𝑿[1]-xExpect[1]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[1].𝑿[1], xExpect[1])
		}

		𝛌Expect := []float64{14, 3}
		if math.Abs(math.Abs(pm.E[0].𝜦)-math.Abs(𝛌Expect[0])) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝜦, 𝛌Expect[0])
		}
		if math.Abs(math.Abs(pm.E[1].𝜦)-math.Abs(𝛌Expect[1])) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[1].𝜦, 𝛌Expect[1])
		}
	})

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

		err = pm.Next(2)
		if err != nil {
			panic(err)
		}

		for i := range pm.E {
			t.Logf("pm.E = %d", i)
			t.Logf("𝑿 = %v", pm.E[i].𝑿)
			t.Logf("𝜦 = %v", pm.E[i].𝜦)
		}

		// result checking

		// lambda = -2
		xExpect := []float64{1.0, 1.0 / 3.0}
		if math.Abs(pm.E[0].𝑿[1]-xExpect[0]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝑿[1], xExpect[0])
		}
		if math.Abs(pm.E[0].𝑿[3]-xExpect[1]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝑿[3], xExpect[1])
		}

		// lambda = -1
		xExpect = []float64{1.0, 1.0 / 4.0}
		if math.Abs(pm.E[1].𝑿[1]-xExpect[0]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[1].𝑿[1], xExpect[0])
		}
		if math.Abs(pm.E[1].𝑿[3]-xExpect[1]) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[1].𝑿[3], xExpect[1])
		}

		𝛌Expect := []float64{-2, -1}
		if math.Abs(math.Abs(pm.E[0].𝜦)-math.Abs(𝛌Expect[0])) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[0].𝜦, 𝛌Expect[0])
		}
		if math.Abs(math.Abs(pm.E[1].𝜦)-math.Abs(𝛌Expect[1])) > eps {
			t.Errorf("Not correct : %e != %e", pm.E[1].𝜦, 𝛌Expect[1])
		}
	})
	t.Run("3x3", func(t *testing.T) {
		T, err := NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			Entry(T, 0, 0, 2),
			Entry(T, 0, 1, 3),

			Entry(T, 1, 0, 3),
			Entry(T, 1, 2, 5),

			Entry(T, 2, 1, 5),
			Entry(T, 2, 2, 2),
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
			IterationMax: 10000000,
			Tolerance:    1e-8,
		})

		if err != nil {
			panic(err)
		}

		err = pm.Next(3)
		for i := range pm.E {
			t.Logf("pm.E = %d", i)
			t.Logf("𝑿 = %v", pm.E[i].𝑿)
			t.Logf("𝜦 = %v", pm.E[i].𝜦)
		}
		if err != nil {
			t.Fatal(err)
		}

		// result checking
		𝛌Expect := []float64{6.9160798, -4.9160798, 2.0}
		for i := range pm.E {
			if math.Abs(math.Abs(pm.E[i].𝜦)-math.Abs(𝛌Expect[i])) > eps {
				t.Errorf("Not correct 𝜦%d: %e != %e", i, pm.E[i].𝜦, 𝛌Expect[i])
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
