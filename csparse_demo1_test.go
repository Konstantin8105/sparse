package sparse

import (
	"bytes"
	"io/ioutil"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

func TestDemo1(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
		buildC(t, "testdata/csparse_demo1_test.c", true)
	})

	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		t.Fatal(err)
	}

	for i := range matrixes {

		if testing.Short() {
			if !strings.Contains(matrixes[i], "bcsstk01") {
				continue
			}
		}

		t.Run("Demo1: "+matrixes[i], func(t *testing.T) {
			// data checking
			b, c, cDur := getCresult(t, matrixes[i])

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

			// load data
			T := Load(&stdin)

			// A compressed-column form of T
			A := Compress(T)

			// Transpose
			AT, err := Transpose(A)
			if err != nil {
				t.Fatal(err)
			}

			// m = # of rows of A
			var m int
			if A != nil {
				m = A.m
			}

			// triplet identify
			T, err = cs_spalloc(m, m, m, true, tripletFormat)
			if err != nil {
				t.Fatal(err)
			}

			for i := 0; i < m; i++ {
				Entry(T, i, i, 1.0)
			}
			Eye := Compress(T)

			// C = A*A'
			C := Multiply(A, AT)

			// D = C + Eye*norm(C,1)
			D, err := Add(C, Eye, 1, Norm(C))
			if err != nil {
				t.Fatal(err)
			}
			Print(D, false)

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
