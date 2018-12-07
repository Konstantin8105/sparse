package sparse

import (
	"bytes"
	"fmt"
	"testing"
)

func BenchmarkIsSym(b *testing.B) {
	for _, size := range []int{30, 100, 300, 1000, 3000} {
		// triplet
		T, err := NewTriplet()
		if err != nil {
			b.Fatal(err)
		}
		// storage
		val := 1.0
		for j := 0; j < size; j++ {
			for k := 0; k < size; k++ {
				// d - distance between diagonal and entry
				d := j - k
				if k > j {
					d = k - j
				}
				if d > 5 { // spacing
					continue
				}

				val = float64(d)
				err = Entry(T, j, k, val)
				if err != nil {
					b.Fatal(err)
				}
			}
		}

		// compress
		A, err := Compress(T)
		if err != nil {
			b.Fatal(err)
		}

		if ok, err := IsSym(A); ok == false || err != nil {
			b.Fatalf("Not sym : %v %v", ok, err)
		}

		b.Run(fmt.Sprintf("%d", size), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_, _ = IsSym(A)
			}
		})
	}
}

func TestIsSym(t *testing.T) {
	tcs := []struct {
		mat     string
		isSym   bool
		isError bool
	}{
		{
			mat:     "0 0 1\n1 0 1\n2 0 1\n2 1 1\n2 2 1\n 1 1 1",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 1\n1 0 1\n2 0 1\n0 2 1\n 1 2 1\n 2 2 1",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 1\n1 0 1\n0 1 2",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 1\n1 0 1\n0 1 1",
			isSym:   true,
			isError: false,
		},
		{
			mat:     "0 0 1\n1 1 1\n2 2 1\n0 0 1\n1 0 1\n0 1 1\n2 0 1\n0 2 1\n3 0 1\n0 3 1",
			isSym:   true,
			isError: false,
		},
		{
			mat:     "0 0 1\n1 1 1\n2 2 1\n0 3 10\n0 0 1\n1 0 1\n0 1 1\n2 0 1\n0 2 1\n3 0 1\n0 3 1\n3 0 20",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 1\n1 1 1\n2 2 1\n3 0 10\n0 0 1\n1 0 1\n0 1 1\n2 0 1\n0 2 1\n3 0 1\n0 3 1\n0 3 20",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 0\n3 3 0\n1 0 5\n0 1 5\n0 2 -1\n2 0 10",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 0\n3 3 0\n1 0 5\n0 1 5\n0 2 -1\n2 0 10",
			isSym:   false,
			isError: true,
		},
		{
			mat:     "0 0 1\n0 1 2\n1 0 2\n2 0 3\n0 2 3\n3 0 4\n0 3 4\n1 1 1\n2 2 1\n0 0 1",
			isSym:   true,
			isError: false,
		},
	}

	for i := range tcs {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			var stdin bytes.Buffer
			stdin.WriteString(tcs[i].mat)
			T, err := Load(&stdin)
			if err != nil {
				t.Fatal(err)
			}
			A, err := Compress(T)
			if err != nil {
				t.Fatal(err)
			}

			ok, err := IsSym(A)

			for i := range A.i {
				if A.i[i] < 0 {
					t.Fatalf("negative index")
				}
			}

			isSame := true

			if ok != tcs[i].isSym {
				isSame = false
			}
			if (err == nil) == tcs[i].isError {
				isSame = false
			}
			if !isSame {
				t.Logf("Is not identical sym:\n%v\n%v", ok, tcs[i].isSym)
				t.Logf("Is not identical error:\n%v\n%v", err, tcs[i].isError)
				t.Fail()
			}
		})
	}
}
