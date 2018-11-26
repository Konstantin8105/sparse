package sparse

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"math"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"
	"time"

	codestyle "github.com/Konstantin8105/cs"
)

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

func Benchmark(b *testing.B) {
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
			o, err := ioutil.ReadFile(matrixes[i])
			if err != nil {
				b.Fatal(err)
			}
			var stdin bytes.Buffer

			b.Run("cs_load", func(b *testing.B) {
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					stdin.Write(o)
					_ = Load(&stdin)
				}
			})

			b.Run("cs_compress", func(b *testing.B) {
				stdin.Write(o)
				T := Load(&stdin)
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = Compress(T)
				}
			})

			b.Run("cs_transpose", func(b *testing.B) {
				stdin.Write(o)
				T := Load(&stdin)
				A, err := Compress(T)
				if err != nil {
					b.Fatal(err)
				}
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = Transpose(A)
				}
			})

			b.Run("cs_add", func(b *testing.B) {
				stdin.Write(o)
				T := Load(&stdin)
				A, err := Compress(T)
				if err != nil {
					b.Fatal(err)
				}
				_, err = Add(A, A, 1, 2)
				if err != nil {
					b.Fatal(err)
				}
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = Add(A, A, 1, 2)
				}
			})

			b.Run("cs_multiply", func(b *testing.B) {
				stdin.Write(o)
				T := Load(&stdin)
				A, err := Compress(T)
				if err != nil {
					b.Fatal(err)
				}
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_ = Multiply(A, A)
				}
			})

			b.Run("cs_gaxpy", func(b *testing.B) {
				stdin.Write(o)
				T := Load(&stdin)
				A, err := Compress(T)
				if err != nil {
					b.Fatal(err)
				}
				x := make([]float64, A.n)
				y := make([]float64, A.m)
				err = Gaxpy(A, x, y)
				if err != nil {
					b.Fatal(err)
				}
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_ = Gaxpy(A, x, y)
				}
			})

			b.Run("demo2", func(b *testing.B) {
				tmpfile, err := ioutil.TempFile("", "example")
				if err != nil {
					b.Fatal(err)
				}
				old := osStdout
				osStdout = tmpfile
				defer func() {
					osStdout = old
				}()

				var stdin bytes.Buffer
				stdin.Write(o)
				prob := get_problem(&stdin, 1e-14, false)
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_ = demo2(prob, false)
				}
			})
		})
	}
}

func TestDemo2(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
		buildC(t, "testdata/csparse_demo2_test.c", true)
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

		if testing.Short() {
			if !strings.Contains(matrixes[i], "bcsstk01") {
				continue
			}
		}

		t.Run("Demo2: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c, cDur := getCresult(t, matrixes[i])
			t.Log("CSparse is ok")

			tmpfile, err := ioutil.TempFile("", "example")
			if err != nil {
				t.Fatal(err)
			}
			old := osStdout
			osStdout = tmpfile
			defer func() {
				osStdout = old
			}()

			start := time.Now() // start timer

			var stdin bytes.Buffer
			stdin.Write(b)
			prob := get_problem(&stdin, 1e-14, true)
			// print_problem(prob)
			result := demo2(prob, true)
			if result {
				fmt.Fprintf(osStdout, "Result demo2 : 1\n")
			} else {
				fmt.Fprintf(osStdout, "Result demo2 : 0\n")
			}

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

// demo2 - solve a linear system using Cholesky, LU, and QR, with various orderings
func demo2(Prob *problem, output bool) bool {
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
		fmt.Fprintf(osStdout, "D is nil\n")
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

	if output {
		fmt.Fprintf(osStdout, "blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank)
	}
	cs_free(D)

	// natural and amd(A'*A)
	for order = 0; order <= 3; order += 3 {
		if output {
			fmt.Fprintf(osStdout, "Order : %d\n", order)
			fmt.Fprintf(osStdout, "M is : %d\n", m)
		}
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		if output {
			fmt.Fprintf(osStdout, "Start order : %d\n", order)
			fmt.Fprintf(osStdout, "QR   ")
			print_order(order, output)
		}
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// min norm(Ax-b) with QR
		ok = cs_qrsol(order, C, x)
		if output {
			fmt.Fprintf(osStdout, "time: %8.2f ", toc(t))
			// print residual
			print_resid(ok, C, x, b, resid)
			if ok {
				for r := 0; r < m; r++ {
					fmt.Fprintf(osStdout, "x[%d] = %10e\n", r, x[r])
				}
			}
		}
	}

	if output {
		fmt.Fprintf(osStdout, "m,n,sprank : %d:%d:%d\n", m, n, sprank)
	}

	if m != n || sprank < n {
		// return if rect. or singular
		return true
	}

	// try all orderings
	for order = 0; order <= 3; order++ {
		if output {
			fmt.Fprintf(osStdout, "Order : %d\n", order)
			fmt.Fprintf(osStdout, "M is : %d\n", m)
		}
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		if output {
			fmt.Fprintf(osStdout, "Start order : %d\n", order)
			fmt.Fprintf(osStdout, "LU   ")
			print_order(order, output)
		}
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// solve Ax=b with LU
		ok = cs_lusol(order, C, x, tol)
		if output {
			fmt.Fprintf(osStdout, "time: %8.2f ", toc(t))
			// print residual
			print_resid(ok, C, x, b, resid)
			if ok {
				for r := 0; r < m; r++ {
					fmt.Fprintf(osStdout, "x[%d] = %10e\n", r, x[r])
				}
			}
		}
	}

	if output {
		fmt.Fprintf(osStdout, "Problem sym is : %d\n", Prob.sym)
	}

	if Prob.sym == 0 {
		return true
	}

	// natural and amd(A+A')
	for order = 0; order <= 1; order++ {
		if output {
			fmt.Fprintf(osStdout, "Order : %d\n", order)
			fmt.Fprintf(osStdout, "M is : %d\n", m)
		}
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		if output {
			fmt.Fprintf(osStdout, "Start order : %d\n", order)
			fmt.Fprintf(osStdout, "Chol ")
			print_order(order, output)
		}
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// solve Ax=b with Cholesky
		ok = cs_cholsol(order, C, x)
		if output {
			fmt.Fprintf(osStdout, "time: %8.2f ", toc(t))
			// print residual
			print_resid(ok, C, x, b, resid)
			if ok {
				for r := 0; r < m; r++ {
					fmt.Fprintf(osStdout, "x[%d] = %10e\n", r, x[r])
				}
			}
		}
	}

	return true
}

// Cholesky update/downdate
func demo3(Prob *problem, output bool) int {
	var A *Cs
	var C *Cs
	var W *Cs
	var WW *Cs
	var WT *Cs
	var W2 *Cs
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
	WW = Multiply(W2, WT)
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

func TestNilCheck(t *testing.T) {
	tcs := []struct {
		name string
		fs   []error
	}{
		{
			name: "Add",
			fs: []error{
				func() error {
					_, err := Add(nil, nil, 0, 0)
					return err
				}(),
				func() error {
					_, err := Add(nil, nil, math.NaN(), math.NaN())
					return err
				}(),
				func() error {
					_, err := Add(nil, nil, math.Inf(0), math.Inf(0))
					return err
				}(),
				func() error {
					var stdin bytes.Buffer
					stdin.WriteString(`0 0 1
0 1 2
1 0 3
1 1 4`)
					T := Load(&stdin)
					_, err := Add(T, T, 0, 0)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T := Load(&s)
					A, err := Compress(T)
					if err != nil {
						panic(err)
					}

					var s2 bytes.Buffer
					s2.WriteString("0 0 1")
					T2 := Load(&s2)
					A2, err := Compress(T2)
					if err != nil {
						panic(err)
					}

					_, err = Add(A, A2, 0, 0)
					return err
				}(),
			},
		},
		{
			name: "Gaxpy",
			fs: []error{
				func() error {
					return Gaxpy(nil, nil, nil)
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T := Load(&s)
					return Gaxpy(T, nil, nil)
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T := Load(&s)
					A, err := Compress(T)
					if err != nil {
						panic(err)
					}
					x := make([]float64, 100)
					y := make([]float64, 80)
					return Gaxpy(A, x, y)
				}(),
			},
		},
		{
			name: "Transpose",
			fs: []error{
				func() error {
					_, err := Transpose(nil)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T := Load(&s)
					_, err := Transpose(T)
					return err
				}(),
			},
		},
	}

	for i := range tcs {
		t.Run(tcs[i].name, func(t *testing.T) {
			for j := range tcs[i].fs {
				if tcs[i].fs[j] == nil {
					t.Fatalf("Error is nil in case %d", j)
				}
				t.Log(tcs[i].fs[j])
			}
		})
	}

	// TODO (KI): modify return types
	if _, err := Add(nil, nil, 0, 0); err == nil {
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
	if _, err := Compress(nil); err == nil {
		t.Errorf("cs_compress: not nil")
	}
	if r := cs_counts(nil, nil, nil, false); r != nil {
		t.Errorf("cs_counts: not nil")
	}
	if _, err := cs_cumsum(nil, nil); err == nil {
		t.Errorf("cs_cumsum: not nil")
	}
	if r := cs_dfs(0, nil, -1, nil, nil, nil); r != -1 {
		t.Errorf("cs_dfs: not nil")
	}
	// TODO
	// if r := cs_bfs(nil, -1, nil, nil, nil, nil, nil, -1); r == true {
	// 	t.Errorf("cs_bfs: not nil")
	// }
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
	if r := Entry(nil, -1, -1, 0); r == true {
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
	// if r := Gaxpy(nil, nil, nil); r == true {
	// 	t.Errorf("cs_gaxpy: not nil")
	// }
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
	if r := Load(nil); r != nil {
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
	// TODO (KI) : no need
	// if r := cs_free(nil); r != nil {
	// 	t.Errorf("cs_free: not nil")
	// }
	if r := cs_realloc(nil, -1, nil); r != nil {
		t.Errorf("cs_realloc: not nil")
	}
	if r := cs_maxtrans(nil, -1); r != nil {
		t.Errorf("cs_maxtrans: not nil")
	}
	if r := Multiply(nil, nil); r != nil {
		t.Errorf("cs_multiply: not nil")
	}
	if r := Norm(nil); r != -1 {
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
	if r := Print(nil, false); r == true {
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
	// TODO(KI)
	// if r := cs_randperm(-1, -1); r != nil {
	// 	t.Errorf("cs_randperm: not nil")
	// }
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
	// TODO
	// if r := cs_vcount(nil, nil); r == true {
	// 	t.Errorf("cs_vcount: not nil")
	// }
	if r := cs_sqr(-1, nil, false); r != nil {
		t.Errorf("cs_sqr: not nil")
	}
	if r := cs_symperm(nil, nil, false); r != nil {
		t.Errorf("cs_symperm: not nil")
	}
	if r := cs_tdfs(-1, -1, nil, nil, nil, nil); r != -1 {
		t.Errorf("cs_tdfs: not nil")
	}
	if _, err := Transpose(nil); err == nil {
		t.Errorf("cs_transpose: not nil")
	}
	if r := cs_updown(nil, -1, nil, nil); r != 0 {
		t.Errorf("cs_updown: not nil")
	}
	if r := cs_usolve(nil, nil); r == true {
		t.Errorf("cs_usolve: not nil")
	}
	// TODO
	// if r := cs_spalloc(-1, -1, -1, false, false); r != nil {
	// 	t.Errorf("cs_spalloc: not nil")
	// }
	if r := cs_sprealloc(nil, -1); r == true {
		t.Errorf("cs_sprealloc: not nil")
	}
	// if r := cs_spfree(nil); r != nil {
	// 	t.Errorf("cs_spfree: not nil")
	// }
	// if r := cs_nfree(nil); r != nil {
	// 	t.Errorf("cs_nfree: not nil")
	// }
	// if r := cs_sfree(nil); r != nil {
	// 	t.Errorf("cs_sfree: not nil")
	// }
	// TODO
	// if r := cs_dalloc(-1, -1); r != nil {
	// 	t.Errorf("cs_dalloc: not nil")
	// }
	// if r := cs_dfree(nil); r != nil {
	// 	t.Errorf("cs_dfree: not nil")
	// }
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

func TestCsCompress(t *testing.T) {
	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		t.Fatal(err)
	}

	f := func(A *Cs) string {
		tmpfile, err := ioutil.TempFile("", "cs_compress")
		if err != nil {
			panic(err)
		}
		old := osStdout
		osStdout = tmpfile
		defer func() {
			osStdout = old
		}()

		// sort
		var ss []string
		for p := 0; p < A.nz; p++ {
			s := fmt.Sprintf("%8d %8d %10e", A.i[p], A.p[p], A.x[p])
			ss = append(ss, s)
		}
		sort.Strings(ss)

		// print
		for i := range ss {
			fmt.Fprintf(osStdout, ss[i]+"\n")
		}

		filename := tmpfile.Name()
		err = tmpfile.Close()
		if err != nil {
			panic(err)
		}
		cb2, err := ioutil.ReadFile(filename)
		if err != nil {
			panic(err)
		}
		return string(cb2)
	}

	for i := range matrixes {

		// TODO : remove, not clear - Why is Fail in QR?
		if strings.Contains(matrixes[i], "bcsstk16") {
			continue
		}

		t.Run(matrixes[i], func(t *testing.T) {
			b, err := ioutil.ReadFile(matrixes[i])
			if err != nil {
				t.Fatal(err)
			}

			var stdin bytes.Buffer
			stdin.Write(b)
			T := Load(&stdin)
			A, err := Compress(T)
			if err != nil {
				t.Fatal(err)
			}

			t.Run("invert", func(t *testing.T) {
				// invert file
				lines := bytes.Split(b, []byte("\n"))
				var buf bytes.Buffer
				for i := range lines {
					line := lines[len(lines)-1-i]
					if len(line) == 0 {
						continue
					}
					buf.Write(line)
					if i == len(lines)-1 {
						continue
					}
					buf.Write([]byte("\n"))
				}
				T2 := Load(&buf)
				if T2 == nil {
					t.Fatalf("T2 is nil")
				}
				A2, err := Compress(T2)
				if err != nil {
					t.Fatal(err)
				}

				if f(A) != f(A2) {
					t.Log(ShowDiff(f(A), f(A2)))
					t.Errorf("matrix is not same")
				}
			})

			t.Run("pair", func(t *testing.T) {
				// pair swapping
				lines := bytes.Split(b, []byte("\n"))
				var buf bytes.Buffer
				for i := 0; i < len(lines)/2; i++ {
					lines[i], lines[len(lines)-1-i] = lines[len(lines)-1-i], lines[i] // swap
				}
				for i := range lines {
					line := lines[len(lines)-1-i]
					if len(line) == 0 {
						continue
					}
					buf.Write(line)
					if i == len(lines)-1 {
						continue
					}
					buf.Write([]byte("\n"))
				}
				T2 := Load(&buf)
				if T2 == nil {
					t.Fatalf("T2 is nil")
				}
				A2, err := Compress(T2)
				if err != nil {
					t.Fatal(err)
				}

				if f(A) != f(A2) {
					t.Log(ShowDiff(f(A), f(A2)))
					t.Errorf("matrix is not same")
				}
			})
		})
	}
}

func TestCodeStyle(t *testing.T) {
	codestyle.All(t)
}

func snapshot(filename string, t *testing.T, f func()) {
	b, err := ioutil.ReadFile(filename)
	if err != nil {
		t.Fatal(err)
	}

	tmpfile, err := ioutil.TempFile("", "example")
	if err != nil {
		t.Fatal(err)
	}
	old := osStdout
	osStdout = tmpfile
	defer func() {
		osStdout = old
	}()

	f()

	file := tmpfile.Name()
	err = tmpfile.Close()
	if err != nil {
		t.Fatal(err)
	}
	b2, err := ioutil.ReadFile(file)
	if err != nil {
		t.Fatal(err)
	}

	if !bytes.Equal(b, b2) {
		t.Fatalf("Results is not same:\n`%s`\n`%s`", string(b), string(b2))
	}
	t.Logf("`%s`", string(b))
}

func TestGaxpy(t *testing.T) {
	snapshot("./testdata/.snapshot.gaxpy", t, func() {
		var s bytes.Buffer
		s.WriteString("0 0 1\n1 0 3\n2 0 5\n0 1 2\n1 1 4\n2 1 6")
		T := Load(&s)
		A, err := Compress(T)
		if err != nil {
			t.Fatal(err)
		}
		x := []float64{7, 8}
		y := []float64{9, 10, 11}
		err = Gaxpy(A, x, y)
		fmt.Fprintf(osStdout, "%v\n", y)
		if err != nil {
			t.Fatal(err)
		}
	})
}

func TestAdd(t *testing.T) {
	snapshot("./testdata/.snapshot.add", t, func() {
		var stdin bytes.Buffer
		stdin.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
		T := Load(&stdin)
		A, err := Compress(T)
		if err != nil {
			t.Fatal(err)
		}
		AT, err := Transpose(A)
		if err != nil {
			t.Fatal(err)
		}
		R, err := Add(A, AT, 1, 2)
		if err != nil {
			t.Fatal(err)
		}
		Print(R, false)
	})
}

func TestCumsum(t *testing.T) {
	t.Run("simple check", func(t *testing.T) {
		p := []int{0, 0, 0, 0, 0}
		c := []int{8, 8, 8, 6}
		nz, err := cs_cumsum(p, c)
		if err != nil {
			t.Fatal(err)
		}
		if nz != 30 {
			t.Fatalf("Not correct result: %d", nz)
		}
	})
	t.Run("overflow check", func(t *testing.T) {
		p := []int{0, 0, 0}
		c := []int{math.MaxInt64, math.MaxInt64}
		_, err := cs_cumsum(p, c)
		if err == nil {
			t.Fatalf("Error for overflow is not happen")
		}
		t.Log(err)
	})
	t.Run("overflow check 2", func(t *testing.T) {
		p := []int{0, 0, 0, 0}
		c := []int{math.MaxInt64, 1, 5}
		_, err := cs_cumsum(p, c)
		if err == nil {
			t.Fatalf("Error for overflow is not happen")
		}
		t.Log(err)
	})
}
