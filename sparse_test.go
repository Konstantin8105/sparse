package sparse

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func buildC(t *testing.T, filename string) {
	// CSparse C source
	csFiles, err := filepath.Glob("CSparse/Source/" + "*.c")
	if err != nil {
		t.Fatal(err)
	}

	// build testdata application :
	//
	// clang -ICSparse/Include/   \
	//  ./CSparse/Source/*.c      \
	//  ./testdata/csparse_test.c \
	//  -lm                       \
	//  -o                        \
	//  ./testdata/csparse_test
	var args []string
	args = append(args, "-ICSparse/Include/")
	args = append(args, csFiles...)
	args = append(args, filename)
	args = append(args, "-lm")
	args = append(args, "-o")
	args = append(args, "testdata/csparse_test")

	cmd := exec.Command("clang", args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	err = cmd.Run()
	if err != nil {
		t.Fatalf("cmd.Run() failed with %s.\n%s\n%s\n",
			err,
			stderr.String(),
			stdout.String(),
		)
	}
}

func getCresult(t *testing.T, matrix string) (in []byte, out string) {
	cmd := exec.Command(
		"./testdata/csparse_test",
	)

	var stdin, stdout, stderr bytes.Buffer
	b, err := ioutil.ReadFile(matrix)
	if err != nil {
		t.Fatal(err)
	}
	stdin.Write(b)
	cmd.Stdin = &stdin
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	err = cmd.Run()
	if err != nil {
		t.Fatalf("cmd.Run() failed with %s.\n%s\n%s\n",
			err,
			stderr.String(),
			stdout.String(),
		)
	}
	return b, stdout.String()
}

func TestDemo1(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
		buildC(t, "testdata/csparse_demo1_test.c")
	})

	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		t.Fatal(err)
	}

	for i := range matrixes {

		// TODO : remove
		if !strings.Contains(matrixes[i], "t1") {
			continue
		}

		t.Run("Demo1: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c := getCresult(t, matrixes[i])
			fmt.Println(c)

			fmt.Println("-------")

			var stdin bytes.Buffer
			stdin.Write(b)
			T := cs_load(&stdin)
			// cs_print(T, false)

			A := cs_compress(T)
			// cs_print(A, false)

			AT := cs_transpose(A, true)
			// cs_print(AT, false)

			var m int
			if A != nil {
				m = A.m
			}
			T = cs_spalloc(m, m, m, true, true)
			for i := 0; i < m; i++ {
				cs_entry(T, i, i, 1.0)
			}
			Eye := cs_compress(T)
			// cs_print(Eye, false)

			C := cs_multiply(A, AT)
			// cs_print(C, false)

			// D = C + Eye*norm(C,1)
			D := cs_add(C, Eye, 1, cs_norm(C))
			cs_print(D, false)
		})
	}
}

func TestDemo2(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
		buildC(t, "testdata/csparse_demo2_test.c")
	})

	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		t.Fatal(err)
	}

	for i := range matrixes {

		// TODO : remove
		// if !strings.Contains(matrixes[i], "t1") {
		// 	continue
		// }

		t.Run("Demo2: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c := getCresult(t, matrixes[i])
			fmt.Println(c)

			fmt.Println("-------")

			var stdin bytes.Buffer
			stdin.Write(b)
			_ = get_problem(&stdin, 1e-14)
		})
	}
}

type problem struct {
	A     *cs
	C     *cs
	sym   int
	x     []float64
	b     []float64
	resid []float64
}

// is_sym - 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
func is_sym(A *cs) int {
	var is_upper bool
	var is_lower bool
	n := A.n
	m := A.m
	var Ap []int = A.p
	var Ai []int = A.i
	if m != n {
		return 0
	}
	is_upper = true
	is_lower = true
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			if Ai[p] > j {
				is_upper = false
			}
			if Ai[p] < j {
				is_lower = false
			}
		}
	}
	if is_upper {
		return 1
	}
	if is_lower {
		return -1
	}
	return 0
}

// dropdiag - true for off-diagonal entries
func dropdiag(i int, j int, aij float64, other interface{}) bool {
	return (i != j)
}

// make_sym - C = A + triu(A,1)'
func make_sym(A *cs) *cs {
	var AT *cs
	var C *cs
	// AT = A'
	AT = cs_transpose(A, true)
	// drop diagonal entries from AT
	cs_fkeep(AT, dropdiag, nil)
	// C = A+AT
	C = cs_add(A, AT, 1, 1)
	cs_spfree(AT)
	return (C)
}

// get_problem - read a problem from a file; use %g for integers to avoid csi conflicts */
func get_problem(f io.Reader, tol float64) *problem {
	var C *cs
	var nz1 int
	var nz2 int
	Prob := new(problem)
	if Prob == nil {
		return nil
	}
	// load triplet matrix T from a file */
	T := cs_load(f)
	A := cs_compress(T)
	// A = compressed-column form of T */
	Prob.A = A
	// clear T */
	cs_spfree(T)
	if !cs_dupl(A) {
		// sum up duplicates */
		return nil
	}
	sym := is_sym(A)
	// determine if A is symmetric */
	Prob.sym = sym
	m := A.m
	n := A.n
	mn := n
	if m > n {
		mn = m
	}
	nz1 = A.p[n]
	// drop zero entries */
	cs_dropzeros(A)
	nz2 = A.p[n]
	if tol > 0 {
		// drop tiny entries (just to test) */
		cs_droptol(A, tol)
	}
	C = func() *cs {
		if sym == 1 {
			return make_sym(A)
		}
		return A
	}()
	// C = A + triu(A,1)', or C=A */
	Prob.C = C
	if C == nil {
		return nil
	}
	fmt.Printf("\n--- Matrix: %g-by-%g, nnz: %g (sym: %g: nnz %g), norm: %8.2e\n",
		float64((m)),
		float64((n)),
		float64((A.p[n])),
		float64((sym)),
		float64((func() int {
			if sym == 1 {
				return C.p[n]
			}
			return 0
		}())), cs_norm(C))
	if nz1 != nz2 {
		fmt.Printf("zero entries dropped: %g\n", float64(nz1-nz2))
	}
	if nz2 != A.p[n] {
		fmt.Printf("tiny entries dropped: %g\n", float64((int32(nz2 - A.p[n]))))
	}
	Prob.b = make([]float64, mn)
	Prob.x = make([]float64, mn)
	Prob.resid = make([]float64, mn)
	return (func() *problem {
		if Prob.b == nil || Prob.x == nil || Prob.resid == nil {
			return nil
		}
		return Prob
	}())
}
