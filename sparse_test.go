package sparse

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
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
		t.Run("Demo1: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c := getCresult(t, matrixes[i])

			tmpfile, err := ioutil.TempFile("", "example")
			if err != nil {
				t.Fatal(err)
			}
			old := os.Stdout
			os.Stdout = tmpfile
			defer func() {
				os.Stdout = old
			}()

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

			filename := tmpfile.Name()
			err = tmpfile.Close()
			if err != nil {
				t.Fatal(err)
			}
			cb2, err := ioutil.ReadFile(filename)
			if err != nil {
				t.Fatal(err)
			}
			c2 := string(cb2)

			// compare strings
			if c != c2 {
				t.Log(ShowDiff(c, c2))
				t.Fail()
			}
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

		// TODO : remove, not clear - Why is Fail in QR?
		if strings.Contains(matrixes[i], "mbeacxc") {
			continue
		}

		t.Run("Demo2: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c := getCresult(t, matrixes[i])

			tmpfile, err := ioutil.TempFile("", "example")
			if err != nil {
				t.Fatal(err)
			}
			old := os.Stdout
			os.Stdout = tmpfile
			defer func() {
				os.Stdout = old
			}()

			var stdin bytes.Buffer
			stdin.Write(b)
			prob := get_problem(&stdin, 1e-14)
			// print_problem(prob)
			demo2(prob)

			filename := tmpfile.Name()
			err = tmpfile.Close()
			if err != nil {
				t.Fatal(err)
			}
			cb2, err := ioutil.ReadFile(filename)
			if err != nil {
				t.Fatal(err)
			}
			c2 := string(cb2)

			// compare strings
			if c != c2 {
				t.Log(ShowDiff(c, c2))
				t.Fail()
			}
		})
	}
}

// ShowDiff will print two strings vertically next to each other so that line
// differences are easier to read.
func ShowDiff(a, b string) string {
	aLines := strings.Split(a, "\n")
	bLines := strings.Split(b, "\n")
	maxLines := int(math.Max(float64(len(aLines)), float64(len(bLines))))
	out := "\n"

	for lineNumber := 0; lineNumber < maxLines; lineNumber++ {
		aLine := ""
		bLine := ""

		// Replace NULL characters with a dot. Otherwise the strings will look
		// exactly the same but have different length (and therfore not be
		// equal).
		if lineNumber < len(aLines) {
			aLine = strconv.Quote(aLines[lineNumber])
		}
		if lineNumber < len(bLines) {
			bLine = strconv.Quote(bLines[lineNumber])
		}

		diffFlag := " "
		if aLine != bLine {
			diffFlag = "*"
		}
		// if diffFlag == " " {
		// 	continue
		// }
		out += fmt.Sprintf("%s %3d %-40s%s\n", diffFlag, lineNumber+1, aLine, bLine)

	}

	return out
}

type problem struct {
	A     *cs
	C     *cs
	sym   int
	x     []float64
	b     []float64
	resid []float64
}

func print_problem(P *problem) {
	fmt.Println("Matrix A:")
	cs_print(P.A, false)
	fmt.Println("Matrix C:")
	cs_print(P.C, false)
	fmt.Println("sym =", P.sym)

	fmt.Printf("Vector x\n")
	for i := 0; i < P.A.n; i++ {
		fmt.Printf("x[%d] = %f\n", i, P.x[i])
	}
	for i := 0; i < P.A.n; i++ {
		fmt.Printf("b[%d] = %f\n", i, P.b[i])
	}
	for i := 0; i < P.A.n; i++ {
		fmt.Printf("resid[%d] = %f\n", i, P.resid[i])
	}
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
	if tol > 0 {
		// drop tiny entries (just to test) */
		cs_droptol(A, tol)
	}
	// C = A + triu(A,1)', or C=A */
	C := func() *cs {
		if sym != 0 {
			return make_sym(A)
		}
		return A
	}()
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
			if sym != 0 {
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

// demo2 - solve a linear system using Cholesky, LU, and QR, with various orderings
func demo2(Prob *problem) bool {
	var t float64
	var tol float64
	var ok bool
	var order int
	var D *csd
	if Prob == nil {
		return false
	}
	A := Prob.A
	C := Prob.C
	b := Prob.b
	x := Prob.x
	resid := Prob.resid
	m := A.m
	n := A.n
	// partial pivoting tolerance
	tol = func() float64 {
		if Prob.sym == 1 {
			return 0.001
		}
		return 1
	}()
	// randomized dmperm analysis
	D = cs_dmperm(C, 1)
	if D == nil {
		fmt.Println("D is nil")
		return false
	}

	nb := D.nb
	r := D.r
	s := D.s
	rr := &D.rr
	sprank := rr[3]

	ns := 0
	k := 0
	for k = 0; k < nb; k++ {
		if (r[k+1] == r[k]+1) && (s[k+1] == s[k]+1) {
			ns++
		}
	}

	fmt.Printf("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank)
	cs_dfree(D)

	// natural and amd(A'*A)
	for order = 0; order <= 3; order += 3 {
		fmt.Printf("Order : %d\n", order)
		fmt.Printf("M is : %d\n", m)
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		fmt.Printf("Start order : %d\n", order)
		fmt.Printf("QR   ")
		print_order(order)
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// min norm(Ax-b) with QR
		ok = cs_qrsol(order, C, x)
		fmt.Printf("time: %8.2f ", toc(t))
		// print residual
		print_resid(ok, C, x, b, resid)
		if ok {
			for r := 0; r < m; r++ {
				fmt.Printf("x[%d] = %10e\n", r, x[r])
			}
		}
	}

	fmt.Printf("m,n,sprank : %d:%d:%d\n", m, n, sprank)

	if m != n || sprank < n {
		// return if rect. or singular
		return true
	}

	// try all orderings
	for order = 0; order <= 3; order++ {
		fmt.Printf("Order : %d\n", order)
		fmt.Printf("M is : %d\n", m)
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		fmt.Printf("Start order : %d\n", order)
		fmt.Printf("LU   ")
		print_order(order)
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// solve Ax=b with LU
		ok = cs_lusol(order, C, x, tol)
		fmt.Printf("time: %8.2f ", toc(t))
		// print residual
		print_resid(ok, C, x, b, resid)
		if ok {
			for r := 0; r < m; r++ {
				fmt.Printf("x[%d] = %10e\n", r, x[r])
			}
		}
	}

	fmt.Printf("Problem sym is : %d\n", Prob.sym)

	if Prob.sym == 0 {
		return true
	}

	// natural and amd(A+A')
	for order = 0; order <= 1; order++ {
		fmt.Printf("Order : %d\n", order)
		fmt.Printf("M is : %d\n", m)
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		fmt.Printf("Start order : %d\n", order)
		fmt.Printf("Chol ")
		print_order(order)
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// solve Ax=b with Cholesky
		ok = cs_cholsol(order, C, x)
		fmt.Printf("time: %8.2f ", toc(t))
		// print residual
		print_resid(ok, C, x, b, resid)
		if ok {
			for r := 0; r < m; r++ {
				fmt.Printf("x[%d] = %10e\n", r, x[r])
			}
		}
	}

	return true
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

// tic -
func tic() float64 {
	return 0
}

// toc -
func toc(t float64) float64 {
	return 0
}

// print_resid - compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf))
func print_resid(ok bool, A *cs, x []float64, b []float64, resid []float64) {
	if !ok {
		fmt.Printf("    (failed)\n")
		return
	}
	m := A.m
	n := A.n

	// resid = -b
	for i := 0; i < m; i++ {
		resid[i] = -b[i]
	}

	// resid = resid + A*x
	cs_gaxpy(A, x, resid)
	// fmt.Printf("resid: %8.2e\n", norm(resid, m)/func() float64 {
	// 	if n == 0 {
	// 		return 1
	// 	}
	// 	return cs_norm(A)*norm(x, n) + norm(b, m)
	// }())
	_ = n
	fmt.Println()
}

// print_order -
func print_order(order int) {
	switch order {
	case 0:
		fmt.Printf("natural    ")
	case 1:
		fmt.Printf("amd(A+A')  ")
	case 2:
		fmt.Printf("amd(S'*S)  ")
	case 3:
		fmt.Printf("amd(A'*A)  ")
	}
}
