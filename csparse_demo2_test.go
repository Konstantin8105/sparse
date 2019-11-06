package sparse

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

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

			tmpfile, err := ioutil.TempFile("", "demo2")
			if err != nil {
				t.Fatal(err)
			}

			start := time.Now() // start timer

			var stdin bytes.Buffer
			stdin.Write(b)
			prob := get_problem(&stdin, 1e-14, true, tmpfile)
			// print_problem(prob)
			result := demo2(prob, true, tmpfile)
			if result {
				fmt.Fprintf(tmpfile, "Result demo2 : 1\n")
			} else {
				fmt.Fprintf(tmpfile, "Result demo2 : 0\n")
			}

			end := time.Now() // end timer

			filename := tmpfile.Name()
			defer func() { _ = os.Remove(filename) }()
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

// demo2 - solve a linear system using Cholesky, LU, and QR, with various orderings
func demo2(Prob *problem, output bool, tmpfile io.Writer) bool {
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
		fmt.Fprintf(tmpfile, "D is nil\n")
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
		fmt.Fprintf(tmpfile, "blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank)
	}
	cs_free(D)

	// natural and amd(A'*A)
	for order = 0; order <= 3; order += 3 {
		if output {
			fmt.Fprintf(tmpfile, "Order : %d\n", order)
			fmt.Fprintf(tmpfile, "M is : %d\n", m)
		}
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		if output {
			fmt.Fprintf(tmpfile, "Start order : %d\n", order)
			fmt.Fprintf(tmpfile, "QR   ")
			print_order(order, output, tmpfile)
		}
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// min norm(Ax-b) with QR
		ok = cs_qrsol(Order(order), C, x)
		if output {
			fmt.Fprintf(tmpfile, "time: %8.2f ", toc(t))
			// print residual
			print_resid(ok, C, x, b, resid, tmpfile)
			if ok {
				for r := 0; r < m; r++ {
					fmt.Fprintf(tmpfile, "x[%d] = %10e\n", r, x[r])
				}
			}
		}
	}

	if output {
		fmt.Fprintf(tmpfile, "m,n,sprank : %d:%d:%d\n", m, n, sprank)
	}

	if m != n || sprank < n {
		// return if rect. or singular
		return true
	}

	// try all orderings
	for order = 0; order <= 3; order++ {
		if output {
			fmt.Fprintf(tmpfile, "Order : %d\n", order)
			fmt.Fprintf(tmpfile, "M is : %d\n", m)
		}
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		if output {
			fmt.Fprintf(tmpfile, "Start order : %d\n", order)
			fmt.Fprintf(tmpfile, "LU   ")
			print_order(order, output, tmpfile)
		}
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// solve Ax=b with LU
		ok = cs_lusol(Order(order), C, x, tol)
		if output {
			fmt.Fprintf(tmpfile, "time: %8.2f ", toc(t))
			// print residual
			print_resid(ok, C, x, b, resid, tmpfile)
			if ok {
				for r := 0; r < m; r++ {
					fmt.Fprintf(tmpfile, "x[%d] = %10e\n", r, x[r])
				}
			}
		}
	}

	if output {
		fmt.Fprintf(tmpfile, "Problem sym is : %d\n", Prob.sym)
	}

	if Prob.sym == 0 {
		return true
	}

	// natural and amd(A+A')
	for order = 0; order <= 1; order++ {
		if output {
			fmt.Fprintf(tmpfile, "Order : %d\n", order)
			fmt.Fprintf(tmpfile, "M is : %d\n", m)
		}
		// if order != 0 && m > 1000 {
		// 	continue
		// }
		if output {
			fmt.Fprintf(tmpfile, "Start order : %d\n", order)
			fmt.Fprintf(tmpfile, "Chol ")
			print_order(order, output, tmpfile)
		}
		// compute right-hand side
		rhs(x, b, m)
		t = tic()
		// solve Ax=b with Cholesky
		ok = cs_cholsol(Order(order), C, x)
		if output {
			fmt.Fprintf(tmpfile, "time: %8.2f ", toc(t))
			// print residual
			print_resid(ok, C, x, b, resid, tmpfile)
			if ok {
				for r := 0; r < m; r++ {
					fmt.Fprintf(tmpfile, "x[%d] = %10e\n", r, x[r])
				}
			}
		}
	}

	return true
}
