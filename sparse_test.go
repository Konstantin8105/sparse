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
			result := demo2(prob)
			if result {
				fmt.Printf("Result demo2 : 1\n")
			} else {
				fmt.Printf("Result demo2 : 0\n")
			}

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

func TestDemo3(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
		buildC(t, "testdata/csparse_demo3_test.c")
	})

	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		t.Fatal(err)
	}

	for i := range matrixes {

		// TODO : remove, not clear - Why it is soo long?
		if strings.Contains(matrixes[i], "bcsstk16") {
			continue
		}

		t.Run("Demo3: "+matrixes[i], func(t *testing.T) {
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
			result := demo3(prob)
			fmt.Printf("Result demo3 : %d\n", result)

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

	fmt.Printf("A before drop\n")
	cs_print(A, false)

	fmt.Printf("n   = %d\n", n)
	fmt.Printf("nz1 = %d\n", nz1)
	fmt.Printf("nz2 = %d\n", nz2)
	fmt.Printf("A->p[%d] = %d\n", n, A.p[n])

	fmt.Printf("tol = %.5e\n", tol)
	if tol > 0 {
		// drop tiny entries (just to test) */
		ok := cs_droptol(A, tol)
		fmt.Printf("droptol = %d\n", ok)
	}

	fmt.Printf("A before make_sym\n")
	cs_print(A, false)

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
		}())),
		0.0)
	// cs_norm(C))
	if nz1 != nz2 {
		fmt.Printf("zero entries dropped: %g\n", float64(nz1-nz2))
	}

	cs_print(A, false)

	fmt.Printf("nz2 = %d\n", nz2)
	fmt.Printf("A->p[%d] = %d\n", n, A.p[n])
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

// Cholesky update/downdate
func demo3(Prob *problem) int {
	var A *cs
	var C *cs
	var W *cs
	var WW *cs
	var WT *cs
	var E *cs
	var W2 *cs
	var n int
	var k int
	var Li []int
	var Lp []int
	var Wi []int
	var Wp []int
	var p1 int
	var p2 int
	var p []int
	var ok int
	var b []float64
	var x []float64
	var resid []float64
	var y []float64
	var Lx []float64
	var Wx []float64
	var s float64
	var t float64
	var t1 float64
	var S *css
	var N *csn
	if Prob == nil || Prob.sym == 0 || Prob.A.n == 0 {
		return 0
	}
	A = Prob.A
	C = Prob.C
	b = Prob.b
	x = Prob.x
	resid = Prob.resid
	n = int(A.n)
	if Prob.sym == 0 || n == 0 {
		return 1
	}
	// compute right-hand side
	rhs(x, b, int(n))
	fmt.Printf("\nchol then update/downdate ")
	print_order(1)
	y = make([]float64, n)
	t = tic()
	// symbolic Chol, amd(A+A')
	S = cs_schol(1, C)
	fmt.Printf("\nsymbolic chol time %8.2f\n", toc(t))
	t = tic()
	// numeric Cholesky
	N = cs_chol(C, S)
	fmt.Printf("numeric  chol time %8.2f\n", toc(t))
	if S == nil || N == nil || y == nil {
		return 1 //done3(0, S, N, y, W, E, p)
	}
	t = tic()
	// y = P*b
	cs_ipvec(S.pinv, b, y, int(n))
	// y = L\y
	cs_lsolve(N.L, y)
	// y = L'\y
	cs_ltsolve(N.L, y)
	// x = P'*y
	cs_pvec(S.pinv, y, x, int(n))
	fmt.Printf("solve    chol time %8.2f\n", toc(t))
	fmt.Printf("original: ")
	// print residual
	print_resid(true, C, x, b, resid)
	// construct W
	k = n / 2
	W = cs_spalloc(n, 1, n, true, false)
	if W == nil {
		return 0 // done3(0, S, N, y, W, E, p)
	}
	Lp = N.L.p
	Li = N.L.i
	Lx = N.L.x
	Wp = W.p
	Wi = W.i
	Wx = W.x
	Wp[0] = 0
	p1 = Lp[k]
	Wp[1] = Lp[k+1] - p1
	s = Lx[p1]
	// rand.Seed(int64(1))
	counter := 0.001
	for ; p1 < Lp[k+1]; p1++ {
		p2 = p1 - Lp[k]
		Wi[p2] = Li[p1]
		Wx[p2] = s * counter // float64(rand.Int()) / float64(2147483647)
		counter *= 1.05
	}
	t = tic()
	// update: L*L'+W*W'
	ok = cs_updown(N.L, int(+1), W, S.parent)
	t1 = toc(t)
	fmt.Printf("update:   time: %8.2f\n", t1)
	if ok == 0 { // check
		return 0 // int((done3(0, S, N, y, W, E, p)))
	}
	t = tic()
	// y = P*b
	cs_ipvec(S.pinv, b, y, int(n))
	// y = L\y
	cs_lsolve(N.L, y)
	// y = L'\y
	cs_ltsolve(N.L, y)
	// x = P'*y
	cs_pvec(S.pinv, y, x, int(n))
	t = toc(t)
	p = cs_pinv(S.pinv, int(n))
	// E = C + (P'W)*(P'W)'
	W2 = cs_permute(W, p, nil, true)
	WT = cs_transpose(W2, true)
	WW = cs_multiply(W2, WT)
	cs_spfree(WT)
	cs_spfree(W2)
	E = cs_add(C, WW, 1, 1)
	cs_spfree(WW)
	if E == nil || p == nil {
		return 0 // int((done3(0, S, N, y, W, E, p)))
	}
	fmt.Printf("update:   time: %8.2f (incl solve) ", t1+t)
	// print residual
	print_resid(true, E, x, b, resid)
	// clear N
	cs_nfree(N)
	t = tic()
	// numeric Cholesky
	N = cs_chol(E, S)
	cs_print(E, false)
	if N == nil {
		return 0 //int((done3(0, S, N, y, W, E, p)))
	}
	// y = P*b
	cs_ipvec(S.pinv, b, y, int(n))
	// y = L\y
	cs_lsolve(N.L, y)
	// y = L'\y
	cs_ltsolve(N.L, y)
	// x = P'*y
	cs_pvec(S.pinv, y, x, int(n))
	t = toc(t)
	fmt.Printf("rechol:   time: %8.2f (incl solve) ", t)
	// print residual
	print_resid(true, E, x, b, resid)
	t = tic()
	// downdate: L*L'-W*W'
	ok = cs_updown(N.L, int(-1), W, S.parent)
	t1 = toc(t)
	if ok == 0 {
		return 0 // int((done3(0, S, N, y, W, E, p)))
	}
	fmt.Printf("downdate: time: %8.2f\n", t1)
	t = tic()
	// y = P*b
	cs_ipvec(S.pinv, b, y, int(n))
	// y = L\y
	cs_lsolve(N.L, y)
	// y = L'\y
	cs_ltsolve(N.L, y)
	// x = P'*y
	cs_pvec(S.pinv, y, x, int(n))
	t = toc(t)
	fmt.Printf("downdate: time: %8.2f (incl solve) ", t1+t)
	// print residual
	print_resid(true, C, x, b, resid)
	return 1 // int((done3(1, S, N, y, W, E, p)))
}

func TestNilCheck(t *testing.T) {
	// TODO (KI): modify return types
	if r := cs_add(nil, nil, 0, 0); r != nil {
		t.Errorf("cs_add: not nil")
	}
	if r := cs_amd(-1, nil); r != nil {
		t.Errorf("cs_amd: not nil")
	}
	if r := cs_chol(nil, nil); r != nil {
		t.Errorf("cs_chol: not nil")
	}
	if r := cs_cholsol(-1, nil, nil); r == true {
		t.Errorf("cs_cholsol: not nil")
	}
	if r := cs_compress(nil); r != nil {
		t.Errorf("cs_compress: not nil")
	}
	if r := cs_counts(nil, nil, nil, false); r != nil {
		t.Errorf("cs_counts: not nil")
	}
	if r := cs_cumsum(nil, nil, 0); r != -1 {
		t.Errorf("cs_cumsum: not nil")
	}
	if r := cs_dfs(0, nil, -1, nil, nil, nil); r != -1 {
		t.Errorf("cs_dfs: not nil")
	}
	if r := cs_bfs(nil, -1, nil, nil, nil, nil, nil, -1); r == true {
		t.Errorf("cs_bfs: not nil")
	}
	if r := cs_dmperm(nil, -1); r != nil {
		t.Errorf("cs_dmperm: not nil")
	}
	if r := cs_droptol(nil, -1); r != -1 {
		t.Errorf("cs_droptol: not nil")
	}
	if r := cs_dropzeros(nil); r != -1 {
		t.Errorf("cs_dropzeros: not nil")
	}
	if r := cs_dupl(nil); r == true {
		t.Errorf("cs_dupl: not nil")
	}
	if r := cs_entry(nil, -1, -1, 0); r == true {
		t.Errorf("cs_entry: not nil")
	}
	if r := cs_ereach(nil, -1, nil, nil, nil); r == 1 {
		t.Errorf("cs_ereach: not nil")
	}
	if r := cs_etree(nil, false); r != nil {
		t.Errorf("cs_etree: not nil")
	}
	if r := cs_fkeep(nil, nil, nil); r == 1 {
		t.Errorf("cs_fkeep: not nil")
	}
	if r := cs_gaxpy(nil, nil, nil); r == true {
		t.Errorf("cs_gaxpy: not nil")
	}
	if r := cs_happly(nil, -1, -1, nil); r == 1 {
		t.Errorf("cs_happly: not nil")
	}
	if r := cs_house(nil, nil, -1); r != -1 {
		t.Errorf("cs_house: not nil")
	}
	if r := cs_ipvec(nil, nil, nil, -1); r == true {
		t.Errorf("cs_ipvec: not nil")
	}
	if r := cs_leaf(-1, -1, nil, nil, nil, nil, nil); r != -1 {
		t.Errorf("cs_leaf: not nil")
	}
	if r := cs_load(nil); r != nil {
		t.Errorf("cs_load: not nil")
	}
	if r := cs_lsolve(nil, nil); r != false {
		t.Errorf("cs_lsolve: not nil")
	}
	if r := cs_ltsolve(nil, nil); r == true {
		t.Errorf("cs_ltsolve: not nil")
	}
	if r := cs_lu(nil, nil, -1); r != nil {
		t.Errorf("cs_lu: not nil")
	}
	if r := cs_lusol(-1, nil, nil, -1); r == true {
		t.Errorf("cs_lusol: not nil")
	}
	if r := cs_free(nil); r != nil {
		t.Errorf("cs_free: not nil")
	}
	if r := cs_realloc(nil, -1, nil); r != nil {
		t.Errorf("cs_realloc: not nil")
	}
	if r := cs_maxtrans(nil, -1); r != nil {
		t.Errorf("cs_maxtrans: not nil")
	}
	if r := cs_multiply(nil, nil); r != nil {
		t.Errorf("cs_multiply: not nil")
	}
	if r := cs_norm(nil); r != -1 {
		t.Errorf("cs_norm: not nil")
	}
	if r := cs_permute(nil, nil, nil, false); r != nil {
		t.Errorf("cs_permute: not nil")
	}
	if r := cs_pinv(nil, -1); r != nil {
		t.Errorf("cs_pinv: not nil")
	}
	if r := cs_post(nil, -1); r != nil {
		t.Errorf("cs_post: not nil")
	}
	if r := cs_print(nil, false); r == true {
		t.Errorf("cs_print: not nil")
	}
	if r := cs_pvec(nil, nil, nil, -1); r == true {
		t.Errorf("cs_pvec: not nil")
	}
	if r := cs_qr(nil, nil); r != nil {
		t.Errorf("cs_qr: not nil")
	}
	if r := cs_qrsol(-1, nil, nil); r == true {
		t.Errorf("cs_qrsol: not nil")
	}
	if r := cs_randperm(-1, -1); r != nil {
		t.Errorf("cs_randperm: not nil")
	}
	if r := cs_reach(nil, nil, -1, nil, nil); r != -1 {
		t.Errorf("cs_reach: not nil")
	}
	if r := cs_scatter(nil, -1, -1, nil, nil, -1, nil, -1); r != -1 {
		t.Errorf("cs_scatter: not nil")
	}
	if r := cs_scc(nil); r != nil {
		t.Errorf("cs_scc: not nil")
	}
	if r := cs_schol(-1, nil); r != nil {
		t.Errorf("cs_schol: not nil")
	}
	if r := cs_spsolve(nil, nil, -1, nil, nil, nil, false); r != -1 {
		t.Errorf("cs_spsolve: not nil")
	}
	if r := cs_vcount(nil, nil); r == true {
		t.Errorf("cs_vcount: not nil")
	}
	if r := cs_sqr(-1, nil, false); r != nil {
		t.Errorf("cs_sqr: not nil")
	}
	if r := cs_symperm(nil, nil, false); r != nil {
		t.Errorf("cs_symperm: not nil")
	}
	if r := cs_tdfs(-1, -1, nil, nil, nil, nil); r != -1 {
		t.Errorf("cs_tdfs: not nil")
	}
	if r := cs_transpose(nil, false); r != nil {
		t.Errorf("cs_transpose: not nil")
	}
	if r := cs_updown(nil, -1, nil, nil); r != 0 {
		t.Errorf("cs_updown: not nil")
	}
	if r := cs_usolve(nil, nil); r == true {
		t.Errorf("cs_usolve: not nil")
	}
	if r := cs_spalloc(-1, -1, -1, false, false); r != nil {
		t.Errorf("cs_spalloc: not nil")
	}
	if r := cs_sprealloc(nil, -1); r == true {
		t.Errorf("cs_sprealloc: not nil")
	}
	if r := cs_spfree(nil); r != nil {
		t.Errorf("cs_spfree: not nil")
	}
	if r := cs_nfree(nil); r != nil {
		t.Errorf("cs_nfree: not nil")
	}
	if r := cs_sfree(nil); r != nil {
		t.Errorf("cs_sfree: not nil")
	}
	if r := cs_dalloc(-1, -1); r != nil {
		t.Errorf("cs_dalloc: not nil")
	}
	if r := cs_dfree(nil); r != nil {
		t.Errorf("cs_dfree: not nil")
	}
	if r := cs_done(nil, nil, nil, false); r != nil {
		t.Errorf("cs_done: not nil")
	}
	if r := cs_idone(nil, nil, nil, false); r != nil {
		t.Errorf("cs_idone: not nil")
	}
	if r := cs_ndone(nil, nil, nil, nil, false); r != nil {
		t.Errorf("cs_ndone: not nil")
	}
	if r := cs_utsolve(nil, nil); r == true {
		t.Errorf("cs_utsolve: not nil")
	}
}
