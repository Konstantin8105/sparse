package sparse_test

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/Konstantin8105/sparse"
	"gonum.org/v1/gonum/mat"
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
	for _, size := range []int{30, 100, 300, 1000, 3000} {
		// triplet
		T, err := sparse.NewTriplet()
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

				val += float64(j*size + k)
				err = sparse.Entry(T, j, k, val)
				if err != nil {
					b.Fatal(err)
				}
			}
		}

		// compress
		A, err := sparse.Compress(T)
		if err != nil {
			b.Fatal(err)
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

		nonZeros := 0
		sparse.Fkeep(A, func(i, j int, x float64) bool {
			nonZeros++
			return true
		})

		b.ResetTimer()
		b.Run(fmt.Sprintf("Sparse:%4d:%5d", size, nonZeros), func(b *testing.B) {
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

				size, _ := A.Dims()

				br := make([]float64, size)
				for j := 0; j < size; j++ {
					br[j] = float64(j + 1)
				}

				x, err := lu.Solve(br)
				if err != nil {
					panic(err)
				}
				_ = x
			}
		})

		// compare with dense matrix
		a := mat.NewDense(size, size, nil)
		val = 1.0
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
				a.Set(j, k, val)
			}
		}
		b.ResetTimer()
		b.Run(fmt.Sprintf("Dense :%4d:%5d", size, nonZeros), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				// solving
				var lu mat.LU

				// factorization
				lu.Factorize(a)

				size, _ := a.Dims()

				br := mat.NewDense(size, 1, nil)
				for j := 0; j < size; j++ {
					br.Set(j, 0, float64(j+1))
				}

				var x mat.Dense
				err := lu.Solve(&x, false, br)
				if err != nil {
					panic(err)
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

	t.Run("Vector b", func(t *testing.T) {
		// Solve next:
		// [ 0 0 0 0 0 ] [ ?  ]   [  ? ]
		// [ 0 1 0 5 0 ] [ x0 ] = [ 11 ]
		// [ 0 0 0 0 0 ] [ ?  ]   [  ? ]
		// [ 0 2 0 3 0 ] [ x1 ]   [  8 ]
		// [ 0 0 0 0 0 ] [ ?  ]   [  ? ]
		T, err := sparse.NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			sparse.Entry(T, 1, 1, 1),
			sparse.Entry(T, 1, 3, 5),
			sparse.Entry(T, 3, 1, 2),
			sparse.Entry(T, 3, 3, 3),
			sparse.Entry(T, 0, 0, 3000),
			sparse.Entry(T, 2, 2, 3000),
			sparse.Entry(T, 4, 4, 3000),
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
		lu := new(sparse.LU)
		err = lu.Factorize(A, 4, 2, 2, 0)
		if err != nil {
			t.Fatal(err)
		}
		b := []float64{0, 11, 0, 8, 0}

		x, err := lu.Solve(b)
		if err != nil {
			t.Fatal(err)
		}

		if math.Abs(x[1]-1) > 1e-8 {
			t.Errorf("x0 = %e", x[1])
		}
		if math.Abs(x[3]-2) > 1e-8 {
			t.Errorf("x1 = %e", x[3])
		}
	})

	t.Run("Add-Ignore-Begin", func(t *testing.T) {
		// Solve next:
		// [ 0 0 0 ] [ ?  ]   [  ? ]
		// [ 0 1 5 ] [ x0 ] = [ 11 ]
		// [ 0 2 3 ] [ x1 ]   [  8 ]
		T, err := sparse.NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			sparse.Entry(T, 1, 1, 1),
			sparse.Entry(T, 1, 2, 5),
			sparse.Entry(T, 2, 1, 2),
			sparse.Entry(T, 2, 2, 3),
			sparse.Entry(T, 0, 0, 3000),
			sparse.Entry(T, 0, 1, 3000),
			sparse.Entry(T, 1, 0, 3000),
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
		lu := new(sparse.LU)
		err = lu.Factorize(A, 0, 0, 0, 0, 0, 0, 0)
		if err != nil {
			t.Fatal(err)
		}
		b := []float64{11, 8}

		x, err := lu.Solve(b)
		if err != nil {
			t.Fatal(err)
		}

		if math.Abs(x[0]-1) > 1e-8 {
			t.Fatalf("x0 = %e", x[0])
		}
		if math.Abs(x[1]-2) > 1e-8 {
			t.Fatalf("x1 = %e", x[1])
		}
	})

	t.Run("Add-Ignore-Begin-End", func(t *testing.T) {
		// Solve next:
		// [ 0 0 0 0 ] [ ?  ]   [  ? ]
		// [ 0 1 5 0 ] [ x0 ] = [ 11 ]
		// [ 0 2 3 0 ] [ x1 ]   [  8 ]
		// [ 0 0 0 0 ] [ ?  ]   [  ? ]
		T, err := sparse.NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			sparse.Entry(T, 1, 1, 1),
			sparse.Entry(T, 1, 2, 5),
			sparse.Entry(T, 2, 1, 2),
			sparse.Entry(T, 2, 2, 3),
			sparse.Entry(T, 0, 0, 3000),
			sparse.Entry(T, 0, 1, 3000),
			sparse.Entry(T, 1, 0, 3000),
			sparse.Entry(T, 3, 1, 3000),
			sparse.Entry(T, 3, 3, 3000),
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
		lu := new(sparse.LU)
		err = lu.Factorize(A, 3, 0)
		if err != nil {
			t.Fatal(err)
		}
		b := []float64{11, 8}

		x, err := lu.Solve(b)
		if err != nil {
			t.Fatal(err)
		}

		if math.Abs(x[0]-1) > 1e-8 {
			t.Fatalf("x0 = %e", x[0])
		}
		if math.Abs(x[1]-2) > 1e-8 {
			t.Fatalf("x1 = %e", x[1])
		}
	})

	t.Run("Add-Ignore-Middle", func(t *testing.T) {
		// Solve next:
		// [ 1 0 5 ] [ x0 ] = [ 11 ]
		// [ 0 0 0 ] [ ?  ]   [  ? ]
		// [ 2 0 3 ] [ x1 ]   [  8 ]
		T, err := sparse.NewTriplet()
		if err != nil {
			panic(err)
		}
		// storage
		errs := []error{
			sparse.Entry(T, 0, 0, 1),
			sparse.Entry(T, 0, 2, 5),
			sparse.Entry(T, 1, 1, 3000),
			sparse.Entry(T, 2, 0, 2),
			sparse.Entry(T, 2, 2, 3),
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
		lu := new(sparse.LU)
		err = lu.Factorize(A, 1)
		if err != nil {
			t.Fatal(err)
		}
		b := []float64{11, 8}

		x, err := lu.Solve(b)
		if err != nil {
			t.Fatal(err)
		}

		if math.Abs(x[0]-1) > 1e-8 {
			t.Fatalf("x0 = %e", x[0])
		}
		if math.Abs(x[1]-2) > 1e-8 {
			t.Fatalf("x1 = %e", x[1])
		}
	})

	t.Run("Add-Ignore-End", func(t *testing.T) {
		// Solve next:
		// [ 1 5 0 ] [ x0 ] = [ 11 ]
		// [ 2 3 0 ] [ x1 ]   [  8 ]
		// [ 0 0 0 ] [ ?  ]   [  ? ]
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
			sparse.Entry(T, 2, 2, 3000),
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
		lu := new(sparse.LU)
		err = lu.Factorize(A, 2, 2, 2, 2, 2, 2, 2)
		if err != nil {
			t.Fatal(err)
		}
		b := []float64{11, 8}

		x, err := lu.Solve(b)
		if err != nil {
			t.Fatal(err)
		}

		if math.Abs(x[0]-1) > 1e-8 {
			t.Fatalf("x0 = %e", x[0])
		}
		if math.Abs(x[1]-2) > 1e-8 {
			t.Fatalf("x1 = %e", x[1])
		}
	})

	t.Run("Wrong-Ignore", func(t *testing.T) {
		// Solve next:
		// [ 1 5 0 ] [ x0 ] = [ 11 ]
		// [ 2 3 0 ] [ x1 ]   [  8 ]
		// [ 0 0 0 ] [ ?  ]   [  ? ]
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
			sparse.Entry(T, 2, 2, 3),
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
		lu := new(sparse.LU)
		err = lu.Factorize(A, 0, 5, 1, -1, 2)
		if err == nil {
			t.Fatalf("Not correct ignore list")
		}
	})
}

func TestMatLU(t *testing.T) {
	matrixes, err := filepath.Glob("./testdata/matrix/" + "*.lu")
	if err != nil {
		t.Fatal(err)
	}
	for i := range matrixes {
		if testing.Short() && i < len(matrixes)-1 {
			continue
		}
		t.Run(matrixes[i], func(t *testing.T) {
			o, err := ioutil.ReadFile(matrixes[i])
			if err != nil {
				t.Fatal(err)
			}
			var stdin bytes.Buffer
			stdin.Write(o)
			T, err := sparse.Load(&stdin)
			if err != nil {
				t.Fatal(err)
			}
			A, err := sparse.Compress(T)
			if err != nil {
				t.Fatal(err)
			}
			r, c := A.Dims()
			t.Logf("size : %4d %4d", r, c)
			if ok, err := sparse.IsSym(A); !ok || err != nil {
				t.Fatalf("matrix is not symmetrical: %v", err)
			}
			t.Log("symmetrical")
			lu := new(sparse.LU)
			err = lu.Factorize(A)
			if err != nil {
				t.Fatal(err)
			}
		})
	}
}
