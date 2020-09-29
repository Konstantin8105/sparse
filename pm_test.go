package sparse

import (
	"bytes"
	"math"
	"testing"
)

// func ExamplePm() {
// 	T, err := NewTriplet()
// 	if err != nil {
// 		panic(err)
// 	}
// 	// storage
// 	errs := []error{
// 		Entry(T, 0, 0, 1.0000),
// 		Entry(T, 0, 1, 2.0000),
//
// 		Entry(T, 1, 0, -2.000),
// 		Entry(T, 1, 1, 1.0000),
// 		Entry(T, 1, 2, 2.0000),
//
// 		Entry(T, 2, 0, 1.0000),
// 		Entry(T, 2, 1, 3.0000),
// 		Entry(T, 2, 2, 1.0000),
// 	}
// 	for i := range errs {
// 		if errs[i] != nil {
// 			panic(errs[i])
// 		}
// 	}
//
// 	// compress
// 	A, err := Compress(T)
// 	if err != nil {
// 		panic(err)
// 	}
//
// 	var pm Pm
// 	err = pm.Factorize(A, &PmConfig{
// 		IterationMax: 1000000000,
// 		Tolerance:    1e-6,
// 	})
// 	if err != nil {
// 		panic(err)
// 	}
//
// 	err = pm.Eigen()
// 	if err != nil {
// 		panic(err)
// 	}
//
// 	// Output:
// }

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

	for i := 0; i < b.N; i++ {
		var pm Pm
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 1000000000,
			Tolerance:    1e-6,
		})
		if err != nil {
			panic(err)
		}

		err = pm.Eigen()
		if err != nil {
			panic(err)
		}
	}
}

func TestPm(t *testing.T) {
	// tolerance
	eps := 1e-7

	compare := func(t *testing.T, v1, v2 float64, name string) {
		v1 = math.Abs(v1)
		v2 = math.Abs(v2)
		if v1 != 0.0 {
			if e := math.Abs(v1-v2) / v1; e > eps {
				t.Errorf("Not correct `%s`: %e != %e. Diff = %e", name, v1, v2, e)
			}
		} else if v2 != 0.0 {
			if e := math.Abs(v1-v2) / v2; e > eps {
				t.Errorf("Not correct `%s`: %e != %e. Diff = %e", name, v1, v2, e)
			}
		} else {
			if e := math.Abs(v1 - v2); e > eps {
				t.Errorf("Not correct `%s`: %e != %e. Diff = %e", name, v1, v2, e)
			}
		}
	}
	compareSlice := func(t *testing.T, v1, v2 []float64, name string) {
		if len(v1) != len(v2) {
			t.Fatalf("Not same len of slices: %d != %d", len(v1), len(v2))
		}
		for i := range v1 {
			compare(t, v1[i], v2[i], name)
		}
	}

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

		var pm Pm
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 1000,
			Tolerance:    1e-8,
		})
		if err != nil {
			panic(err)
		}

		err = pm.Eigen()
		if err != nil {
			panic(err)
		}

		t.Logf("𝑿 = %v", pm.𝑿)
		t.Logf("𝜦 = %v", pm.𝜦)

		// result checking
		𝛌Expect := 14.0
		xExpect := []float64{1.0, 0.2}
		compare(t, pm.𝜦, 𝛌Expect, "𝛌")
		compareSlice(t, pm.𝑿, xExpect, "𝑿")
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

		var pm Pm
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 1000,
			Tolerance:    1e-8,
		}, 0, 4, 2, 0, 2, 4)
		if err != nil {
			t.Fatal(err)
		}

		err = pm.Eigen()
		if err != nil {
			t.Fatal(err)
		}

		t.Logf("𝑿 = %v", pm.𝑿)
		t.Logf("𝜦 = %v", pm.𝜦)

		// result checking
		𝛌Expect := -2.0
		xExpect := []float64{0.0, 1.0, 0.0, 1.0 / 3.0, 0.0}
		compare(t, pm.𝜦, 𝛌Expect, "𝛌")
		compareSlice(t, pm.𝑿, xExpect, "𝑿")
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

		var pm Pm
		err = pm.Factorize(A, &PmConfig{
			IterationMax: 10000000,
			Tolerance:    1e-8,
		})

		if err != nil {
			panic(err)
		}

		err = pm.Eigen()
		if err != nil {
			t.Fatal(err)
		}

		t.Logf("𝑿 = %v", pm.𝑿)
		t.Logf("𝜦 = %v", pm.𝜦)

		// result checking
		𝛌Expect := 6.916079682949667
		compare(t, pm.𝜦, 𝛌Expect, "𝛌")
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

		var pm Pm
		err = pm.Factorize(A, &PmConfig{IterationMax: 1})
		if err != nil {
			panic(err)
		}

		err = pm.Eigen()
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
		var pm Pm
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
		var pm Pm
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
		var pm Pm
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
		var pm Pm
		if err = pm.Factorize(A, nil, -1, 2, 100, 9); err == nil {
			t.Fatalf("not check")
		}
		t.Log(err)
	})
}
