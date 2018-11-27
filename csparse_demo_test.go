package sparse

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os/exec"
	"path/filepath"
	"testing"
	"time"
)

type problem struct {
	A     *Cs
	C     *Cs
	sym   int
	x     []float64
	b     []float64
	resid []float64
}

func print_problem(P *problem) {
	fmt.Fprintf(osStdout, "Matrix A:\n")
	Print(P.A, false)
	fmt.Fprintf(osStdout, "Matrix C:\n")
	Print(P.C, false)
	fmt.Fprintf(osStdout, "sym = %v\n", P.sym)

	fmt.Fprintf(osStdout, "Vector x\n")
	for i := 0; i < P.A.n; i++ {
		fmt.Fprintf(osStdout, "x[%d] = %f\n", i, P.x[i])
	}
	for i := 0; i < P.A.n; i++ {
		fmt.Fprintf(osStdout, "b[%d] = %f\n", i, P.b[i])
	}
	for i := 0; i < P.A.n; i++ {
		fmt.Fprintf(osStdout, "resid[%d] = %f\n", i, P.resid[i])
	}
}

// norm - infinity-norm of x
func norm(x []float64, n int) float64 {
	var normx float64
	for i := 0; i < n; i++ {
		normx = func() float64 {
			if normx > math.Abs(x[i]) {
				return normx
			}
			return math.Abs(x[i])
		}()
	}
	return normx
}

// rhs - create a right-hand side
func rhs(x []float64, b []float64, m int) {
	for i := 0; i < m; i++ {
		b[i] = 1 + float64(i)/float64(m)
	}
	for i := 0; i < m; i++ {
		x[i] = b[i]
	}
}

// is_sym - 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise
func is_sym(A *Cs) int {
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
func make_sym(A *Cs) *Cs {
	// AT = A'
	AT, err := Transpose(A)
	if err != nil {
		return nil
	}
	// drop diagonal entries from AT
	cs_fkeep(AT, dropdiag, nil)
	// C = A+AT
	C, err := Add(A, AT, 1, 1)
	if err != nil {
		panic(err)
	}
	cs_free(AT)
	return (C)
}

// print_resid - compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf))
func print_resid(ok bool, A *Cs, x []float64, b []float64, resid []float64) {
	if !ok {
		fmt.Fprintf(osStdout, "    (failed)\n")
		return
	}
	m := A.m
	n := A.n

	// resid = -b
	for i := 0; i < m; i++ {
		resid[i] = -b[i]
	}

	// resid = resid + A*x
	Gaxpy(A, x, resid)
	// fmt.Fprintf(osStdout,"resid: %8.2e\n", norm(resid, m)/func() float64 {
	// 	if n == 0 {
	// 		return 1
	// 	}
	// 	return cs_norm(A)*norm(x, n) + norm(b, m)
	// }())
	_ = n
	fmt.Fprintf(osStdout, "\n")
}

// print_order -
func print_order(order int, output bool) {
	if !output {
		return
	}
	switch order {
	case 0:
		fmt.Fprintf(osStdout, "natural    ")
	case 1:
		fmt.Fprintf(osStdout, "amd(A+A')  ")
	case 2:
		fmt.Fprintf(osStdout, "amd(S'*S)  ")
	case 3:
		fmt.Fprintf(osStdout, "amd(A'*A)  ")
	}

}

// get_problem - read a problem from a file; use %g for integers to avoid csi conflicts */
func get_problem(f io.Reader, tol float64, output bool) *problem {
	var nz1 int
	var nz2 int
	Prob := new(problem)
	if Prob == nil {
		return nil
	}
	// load triplet matrix T from a file */
	T := Load(f)
	A, err := Compress(T)
	if err != nil {
		panic(err)
	}
	// A = compressed-column form of T */
	Prob.A = A
	// clear T */
	cs_free(T)
	if !cs_dupl(A) {
		// sum up duplicates */
		return nil
	}
	// determine if A is symmetric */
	sym := is_sym(A)
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

	if output {
		fmt.Fprintf(osStdout, "n   = %d\n", n)
		fmt.Fprintf(osStdout, "nz1 = %d\n", nz1)
		fmt.Fprintf(osStdout, "nz2 = %d\n", nz2)
		fmt.Fprintf(osStdout, "A->p[%d] = %d\n", n, A.p[n])

		fmt.Fprintf(osStdout, "A before drop\n")
		Print(A, false)

		fmt.Fprintf(osStdout, "tol = %.5e\n", tol)
	}
	if tol > 0 {
		// drop tiny entries (just to test) */
		ok := cs_droptol(A, tol)
		if output {
			fmt.Fprintf(osStdout, "droptol = %d\n", ok)
		}
	}

	if output {
		fmt.Fprintf(osStdout, "A before make_sym\n")
		Print(A, false)
	}

	// C = A + triu(A,1)', or C=A */
	C := func() *Cs {
		if sym != 0 {
			return make_sym(A)
		}
		return A
	}()
	Prob.C = C
	if C == nil {
		return nil
	}
	if output {
		fmt.Fprintf(osStdout, "\n--- Matrix: %g-by-%g, nnz: %g (sym: %g: nnz %g), norm: %8.2e\n",
			float64((m)),
			float64((n)),
			float64((A.p[n])),
			float64((sym)),
			float64((func() int {
				if sym != 0 {
					return C.p[n]
				}
				return 0
			}())),
			0.0)
		// cs_norm(C))
		if nz1 != nz2 {
			fmt.Fprintf(osStdout, "zero entries dropped: %g\n", float64(nz1-nz2))
		}

		Print(A, false)

		fmt.Fprintf(osStdout, "nz2 = %d\n", nz2)
		fmt.Fprintf(osStdout, "A->p[%d] = %d\n", n, A.p[n])
		if nz2 != A.p[n] {
			fmt.Fprintf(osStdout, "tiny entries dropped: %g\n", float64((int32(nz2 - A.p[n]))))
		}
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

// tic -
func tic() float64 {
	return 0
}

// toc -
func toc(t float64) float64 {
	return 0
}

func buildC(t *testing.T, filename string, output bool) {
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
	if output {
		args = append(args, "-DPRINT")
	}
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

func getCresult(t *testing.T, matrix string) (in []byte, out string, dur time.Duration) {
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
	start := time.Now()
	err = cmd.Run()
	end := time.Now()
	if err != nil {
		t.Fatalf("cmd.Run() failed with %s.\n%s\n%s\n",
			err,
			stderr.String(),
			stdout.String(),
		)
	}
	return b, stdout.String(), end.Sub(start)
}
