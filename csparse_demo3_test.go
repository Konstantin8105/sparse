package sparse

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func TestDemo3(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
		buildC(t, "testdata/csparse_demo3_test.c", true)
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

		if testing.Short() {
			if !strings.Contains(matrixes[i], "bcsstk01") {
				continue
			}
		}

		t.Run("Demo3: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c, cDur := getCresult(t, matrixes[i])
			t.Log("CSparse is ok")

			start := time.Now() // start timer

			tmpfile, err := ioutil.TempFile("", "example")
			if err != nil {
				t.Fatal(err)
			}
			old := osStdout
			osStdout = tmpfile
			defer func() {
				osStdout = old
			}()

			var stdin bytes.Buffer
			stdin.Write(b)
			prob := get_problem(&stdin, 1e-14, true)
			// print_problem(prob)
			result := demo3(prob, true)
			fmt.Fprintf(osStdout, "Result demo3 : %d\n", result)

			end := time.Now() // end timer

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

			t.Logf("Compare bench:\nC program  : %15d\nGo program : %15d\n",
				cDur, end.Sub(start))
		})
	}
}

// Cholesky update/downdate
func demo3(Prob *problem, output bool) int {
	var A *Matrix
	var C *Matrix
	var W *Matrix
	var WW *Matrix
	var WT *Matrix
	var W2 *Matrix
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
	if output {
		fmt.Fprintf(osStdout, "\nchol then update/downdate ")
		print_order(1, output)
	}
	y = make([]float64, n)
	t = tic()
	// symbolic Chol, amd(A+A')
	S = cs_schol(1, C)
	if output {
		fmt.Fprintf(osStdout, "\nsymbolic chol time %8.2f\n", toc(t))
	}
	t = tic()
	// numeric Cholesky
	N = cs_chol(C, S)
	if output {
		fmt.Fprintf(osStdout, "numeric  chol time %8.2f\n", toc(t))
	}
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
	if output {
		fmt.Fprintf(osStdout, "solve    chol time %8.2f\n", toc(t))
		fmt.Fprintf(osStdout, "original: ")
		// print residual
		print_resid(true, C, x, b, resid)
	}
	// construct W
	k = n / 2
	W, err := cs_spalloc(n, 1, n, true, cscFormat)
	if err != nil {
		panic(err)
	}
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
	if output {
		fmt.Fprintf(osStdout, "update:   time: %8.2f\n", t1)
	}
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
	WT, err = Transpose(W2)
	if err != nil {
		panic(err)
	}
	WW, err = Multiply(W2, WT)
	if err != nil {
		panic(err)
	}
	cs_free(WT)
	cs_free(W2)
	E, err := Add(C, WW, 1, 1)
	if err != nil {
		panic(err)
	}
	cs_free(WW)
	if E == nil || p == nil {
		return 0 // int((done3(0, S, N, y, W, E, p)))
	}
	if output {
		fmt.Fprintf(osStdout, "update:   time: %8.2f (incl solve) ", t1+t)
		// print residual
		print_resid(true, E, x, b, resid)
	}
	// clear N
	cs_free(N)
	t = tic()
	// numeric Cholesky
	N = cs_chol(E, S)
	if output {
		Print(E, false)
	}
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
	if output {
		fmt.Fprintf(osStdout, "rechol:   time: %8.2f (incl solve) ", t)
		// print residual
		print_resid(true, E, x, b, resid)
	}
	t = tic()
	// downdate: L*L'-W*W'
	ok = cs_updown(N.L, int(-1), W, S.parent)
	t1 = toc(t)
	if ok == 0 {
		return 0 // int((done3(0, S, N, y, W, E, p)))
	}
	if output {
		fmt.Fprintf(osStdout, "downdate: time: %8.2f\n", t1)
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
	if output {
		fmt.Fprintf(osStdout, "downdate: time: %8.2f (incl solve) ", t1+t)
		// print residual
		print_resid(true, C, x, b, resid)
	}
	return 1 // int((done3(1, S, N, y, W, E, p)))
}
