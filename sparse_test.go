package sparse

import (
	"bytes"
	"fmt"
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
		if !strings.Contains(matrixes[i], "t1") {
			continue
		}
		t.Run("Demo2: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c := getCresult(t, matrixes[i])
			fmt.Println(c)

			fmt.Println("-------")

			_ = b
		})
	}
}

type problem struct {
	A     *cs
	C     *cs
	sym   bool
	x     *float64
	b     *float64
	resid *float64
}

// is_sym - 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
func is_sym(A []cs) noarch.PtrdiffT {
	var is_upper noarch.PtrdiffT
	var is_lower noarch.PtrdiffT
	var j noarch.PtrdiffT
	var p noarch.PtrdiffT
	var n noarch.PtrdiffT = noarch.PtrdiffT(A[0].n)
	var m noarch.PtrdiffT = noarch.PtrdiffT(A[0].m)
	var Ap []noarch.PtrdiffT = A[0].p
	var Ai []noarch.PtrdiffT = A[0].i
	if m != n {
		return noarch.PtrdiffT((0))
	}
	is_upper = 1
	is_lower = 1
	for j = 0; j < n; j++ {
		for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			if Ai[p] > j {
				is_upper = 0
			}
			if Ai[p] < j {
				is_lower = 0
			}
		}
	}
	return noarch.PtrdiffT((func() int {
		if bool(noarch.PtrdiffT(is_upper)) {
			return 1
		}
		return (func() int {
			if bool(noarch.PtrdiffT(is_lower)) {
				return -1
			}
			return 0
		}())
	}()))
}

// dropdiag - true for off-diagonal entries
func dropdiag(i noarch.PtrdiffT, j noarch.PtrdiffT, aij float64, other interface{}) noarch.PtrdiffT {
	return noarch.PtrdiffT((i != j))
}

// make_sym - C = A + triu(A,1)'
func make_sym(A []cs) []cs {
	var AT []cs
	var C []cs
	// AT = A' */
	AT = cs_transpose(A, 1)
	// drop diagonal entries from AT */
	cs_fkeep(AT, dropdiag, nil)
	// C = A+AT */
	C = cs_add(A, AT, 1, 1)
	cs_spfree(AT)
	return (C)
}

// free_problem - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/testdata/csparse_demo2_test.c:47
// free a problem */
func free_problem(Prob *problem) *problem {
	return nil
}

// get_problem - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/testdata/csparse_demo2_test.c:59
// read a problem from a file; use %g for integers to avoid csi conflicts */
func get_problem(f *noarch.File, tol float64) []problem {
	var T []cs
	var A []cs
	var C []cs
	var sym noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var mn noarch.PtrdiffT
	var nz1 noarch.PtrdiffT
	var nz2 noarch.PtrdiffT
	var Prob []problem
	// Warning (*ast.UnaryExprOrTypeTraitExpr):  $GOPATH/src/github.com/Konstantin8105/sparse/testdata/csparse_demo2_test.c:64 :Cannot determine sizeof : |problem|. err = Cannot canculate `struct` sizeof for `string`. Cannot determine sizeof : |ptrdiff_t|. err = error in array size
	Prob = cs_calloc(1, uint(0)).([]problem)
	if Prob == nil {
		return nil
	}
	// load triplet matrix T from a file */
	T = cs_load(f)
	A = cs_compress(T)
	// A = compressed-column form of T */
	Prob[0].A = A
	// clear T */
	cs_spfree(T)
	if bool(noarch.NotNoarch.PtrdiffT(cs_dupl(A))) {
		// sum up duplicates */
		return (free_problem(Prob))
	}
	sym = is_sym(A)
	// determine if A is symmetric */
	Prob[0].sym = sym
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	mn = noarch.PtrdiffT(func() int32 {
		if m > n {
			return int32(noarch.PtrdiffT((m)))
		}
		return int32(noarch.PtrdiffT((n)))
	}() / 8)
	nz1 = A[0].p[n]
	// drop zero entries */
	cs_dropzeros(A)
	nz2 = A[0].p[n]
	if tol > 0 {
		// drop tiny entries (just to test) */
		cs_droptol(A, tol)
	}
	C = func() []cs {
		if bool(noarch.PtrdiffT(sym)) {
			return make_sym(A)
		}
		return A
	}()
	// C = A + triu(A,1)', or C=A */
	Prob[0].C = C
	if C == nil {
		return (free_problem(Prob))
	}
	noarch.Printf([]byte("\n--- Matrix: %g-by-%g, nnz: %g (sym: %g: nnz %g), norm: %8.2e\n\x00"), float64(noarch.PtrdiffT(m)), float64(noarch.PtrdiffT(n)), float64(noarch.PtrdiffT((A[0].p[n]))), float64(noarch.PtrdiffT(sym)), float64((func() int32 {
		if bool(noarch.PtrdiffT(sym)) {
			return int32(noarch.PtrdiffT(C[0].p[n]))
		}
		return 0
	}())), cs_norm(C))
	if nz1 != nz2 {
		noarch.Printf([]byte("zero entries dropped: %g\n\x00"), float64((int32(nz1 - nz2))))
	}
	if nz2 != A[0].p[n] {
		noarch.Printf([]byte("tiny entries dropped: %g\n\x00"), float64((int32(nz2 - A[0].p[n]))))
	}
	Prob[0].b = cs_malloc(noarch.PtrdiffT(mn), uint(8)).([]float64)
	Prob[0].x = cs_malloc(noarch.PtrdiffT(mn), uint(8)).([]float64)
	Prob[0].resid = cs_malloc(noarch.PtrdiffT(mn), uint(8)).([]float64)
	return (func() []problem {
		if Prob[0].b == nil || Prob[0].x == nil || Prob[0].resid == nil {
			return free_problem(Prob)
		}
		return Prob
	}())
}

// main - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/testdata/csparse_demo2_test.c:92
// cs_demo2: read a matrix and solve a linear system
func main() {
	defer noarch.AtexitRun()
	var Prob []problem = get_problem(noarch.Stdin, 1e-14)
	// demo2 (Prob) ;
	free_problem(Prob)
	os.Exit((0))
}
