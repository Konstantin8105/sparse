package sparse_test

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
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
	min, max := math.MaxFloat64, 0.0
	_, err = sparse.Fkeep(A, func(i, j int, x float64) bool {
		if i == j { // diagonal
			if math.Abs(x) > max {
				max = math.Abs(x)
			}
			if math.Abs(x) < min {
				min = math.Abs(x)
			}
		}
		// keep entry
		return true
	})
	if err != nil {
		panic(err)
	}
	if min == 0 {
		panic("singular: zero entry on diagonal")
	}
	if max/min > 1e18 {
		panic(fmt.Sprintf("singular: max/min diagonal entry: %v", max/min))
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

func BenchmarkLU(b *testing.B) {
	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		b.Fatal(err)
	}

	for i := range matrixes {
		if testing.Short() {
			if !strings.Contains(matrixes[i], "bcsstk01") {
				continue
			}
		}

		b.Run(matrixes[i], func(b *testing.B) {
			// triplet
			T, err := sparse.NewTriplet()
			if err != nil {
				b.Fatal(err)
			}
			// compress
			A, err := sparse.Compress(T)
			if err != nil {
				b.Fatal(err)
			}

			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				// solving
				lu := new(sparse.LU)

				// order
				lu.Order(sparse.AmdNatural)

				// factorization
				err = lu.Factorize(A)
				if err != nil {
					b.Fatal(err)
				}

				size, _ := A.Size()

				br := make([]float64, size)
				for j := 0; j < size; j++ {
					br[j] = float64(j + 1)
				}

				x, err := lu.Solve(br)
				if err != nil {
					b.Fatal(err)
				}
				_ = x
			}
		})
	}
}

func TestLU(t *testing.T) {
	sizes := []int{1, 2, 3, 4, 256, 512, 1024, 2048}
	tol := 1e-3
	for i := range sizes {
		t.Run(fmt.Sprintf("Size%d", sizes[i]), func(t *testing.T) {
			size := sizes[i]
			// triplet
			T, err := sparse.NewTriplet()
			if err != nil {
				t.Fatal(err)
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

					val += float64(j*size + k)
					err = sparse.Entry(T, j, k, val)
					if err != nil {
						t.Fatal(err)
					}
				}
			}

			// compress
			A, err := sparse.Compress(T)
			if err != nil {
				t.Fatal(err)
			}

			// singinal check
			min, max := math.MaxFloat64, 0.0
			_, err = sparse.Fkeep(A, func(i, j int, x float64) bool {
				if i == j { // diagonal
					if math.Abs(x) > max {
						max = math.Abs(x)
					}
					if math.Abs(x) < min {
						min = math.Abs(x)
					}
				}
				// keep entry
				return true
			})
			if err != nil {
				panic(err)
			}
			if min == 0 {
				panic("singular: zero entry on diagonal")
			}
			if max/min > 1e18 {
				panic(fmt.Sprintf("singular: max/min diagonal entry: %v", max/min))
			}

			// solving
			lu := new(sparse.LU)

			// order
			lu.Order(sparse.AmdNatural)

			// factorization
			err = lu.Factorize(A)
			if err != nil {
				t.Fatal(err)
			}

			b := make([]float64, size)
			for j := 0; j < size; j++ {
				b[j] = float64(j + 1)
			}

			x, err := lu.Solve(b)
			if err != nil {
				t.Fatal(err)
			}

			for j := 0; j < size; j++ {
				b[j] = -b[j]
			}

			err = sparse.Gaxpy(A, x, b)
			if err != nil {
				t.Fatal(err)
			}
			max = 0.0
			for j := range b {
				if math.Abs(b[j]) > max {
					max = math.Abs(b[j])
				}
			}
			if max > tol {
				t.Fatalf("Tolerance problem : %.5e", max)
			}
		})
	}

	t.Run("Wrong", func(t *testing.T) {
		// solving
		lu := new(sparse.LU)

		// factorization
		err := lu.Factorize(nil)
		if err == nil {
			t.Fatal("nil factorization")
		}

		_, err = lu.Solve(nil)
		if err == nil {
			t.Fatal("nil solve")
		}
	})

	t.Run("Wrong-Triplet", func(t *testing.T) {
		// solving
		lu := new(sparse.LU)

		// factorization
		T, err := sparse.NewTriplet()
		if err != nil {
			t.Fatal(err)
		}

		err = lu.Factorize((*sparse.Matrix)(T))
		if err == nil {
			t.Fatal("nil factorization")
		}
	})
}
