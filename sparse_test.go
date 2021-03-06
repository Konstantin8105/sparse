package sparse

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	codestyle "github.com/Konstantin8105/cs"
)

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
					_, _ = Load(&stdin)
				}
			})

			b.Run("cs_compress", func(b *testing.B) {
				stdin.Write(o)
				T, err := Load(&stdin)
				if err != nil {
					b.Fatal(err)
				}
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = Compress(T)
				}
			})

			b.Run("cs_transpose", func(b *testing.B) {
				stdin.Write(o)
				T, err := Load(&stdin)
				if err != nil {
					b.Fatal(err)
				}
				A, err := Compress(T)
				if err != nil {
					b.Fatal(err)
				}
				_, err = Transpose(A)
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
				T, err := Load(&stdin)
				if err != nil {
					b.Fatal(err)
				}
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
				T, err := Load(&stdin)
				if err != nil {
					b.Fatal(err)
				}
				A, err := Compress(T)
				if err != nil {
					b.Fatal(err)
				}
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_, _ = Multiply(A, A)
				}
			})

			b.Run("cs_gaxpy", func(b *testing.B) {
				stdin.Write(o)
				T, err := Load(&stdin)
				if err != nil {
					b.Fatal(err)
				}
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
				tmpfile, err := ioutil.TempFile("", "bench_demo2")
				if err != nil {
					b.Fatal(err)
				}

				var stdin bytes.Buffer
				stdin.Write(o)
				prob := get_problem(&stdin, 1e-14, false, tmpfile)
				b.ResetTimer()
				for i := 0; i < b.N; i++ {
					_ = demo2(prob, false, tmpfile)
				}
			})
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
					stdin.WriteString("0 0 1\n 0 1 2\n 1 0 3\n 1 1 4")
					T, err := Load(&stdin)
					if err != nil {
						panic(err)
					}
					_, err = Add((*Matrix)(T), (*Matrix)(T), 0, 0)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					A, err := Compress(T)
					if err != nil {
						panic(err)
					}

					var s2 bytes.Buffer
					s2.WriteString("0 0 1")
					T2, err := Load(&s2)
					if err != nil {
						panic(err)
					}
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
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					return Gaxpy((*Matrix)(T), nil, nil)
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
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
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					_, err = Transpose((*Matrix)(T))
					return err
				}(),
			},
		},
		{
			name: "Dupl",
			fs: []error{
				func() error {
					err := Dupl(nil)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					err = Dupl((*Matrix)(T))
					return err
				}(),
			},
		},
		{
			name: "Entry",
			fs: []error{
				func() error {
					err := Entry(nil, -1, -1, math.NaN())
					return err
				}(),
				func() error {
					err := Entry(nil, -1, -1, math.Inf(0))
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					A, err := Compress(T)
					if err != nil {
						panic(err)
					}
					err = Entry((*Triplet)(A), 1, 1, 12)
					return err
				}(),
			},
		},
		{
			name: "Multiply",
			fs: []error{
				func() error {
					_, err := Multiply(nil, nil)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					_, err = Multiply((*Matrix)(T), (*Matrix)(T))
					return err
				}(),
				func() error {
					var s bytes.Buffer
					var A, B *Matrix
					var err error
					{
						s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
						T, err := Load(&s)
						if err != nil {
							panic(err)
						}
						A, err = Compress(T)
						if err != nil {
							panic(err)
						}
					}
					{
						s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4\n12 12 1")
						T, err := Load(&s)
						if err != nil {
							panic(err)
						}
						B, err = Compress(T)
						if err != nil {
							panic(err)
						}
					}
					_, err = Multiply(A, B)
					return err
				}(),
			},
		},
		{
			name: "IsSym",
			fs: []error{
				func() error {
					_, err := IsSym(nil)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					// triplet in input
					_, err = IsSym((*Matrix)(T))
					return err
				}(),
				func() error {
					var s bytes.Buffer
					// empty matrix
					s.WriteString("")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					A, err := Compress(T)
					if err != nil {
						panic(err)
					}
					_, err = IsSym(A)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					// rectangle matrix
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4\n3 7 2")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					A, err := Compress(T)
					if err != nil {
						panic(err)
					}
					_, err = IsSym(A)
					return err
				}(),
			},
		},
		{
			name: "Fkeep",
			fs: []error{
				func() error {
					_, err := Fkeep(nil, nil)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					// triplet in input
					_, err = Fkeep((*Matrix)(T), nil)
					return err
				}(),
			},
		},
		{
			name: "Load",
			fs: []error{
				func() error {
					_, err := Load(nil)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n1 1 1 1\n1 1 4")
					_, err := Load(&s)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n1  1\n1 1 4")
					_, err := Load(&s)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n1 a 1\n1 1 4")
					_, err := Load(&s)
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n-1 -1 1\n1 1 4")
					_, err := Load(&s)
					return err
				}(),
			},
		},
		{
			name: "Copy",
			fs: []error{
				func() error {
					A := new(Matrix)
					_, err := A.Copy()
					return err
				}(),
				func() error {
					var s bytes.Buffer
					s.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
					T, err := Load(&s)
					if err != nil {
						panic(err)
					}
					// triplet in input
					_, err = ((*Matrix)(T)).Copy()
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
	if r := cs_amd(100, nil); r != nil {
		t.Errorf("cs_amd: not nil")
	}
	if r := cs_chol(nil, nil); r != nil {
		t.Errorf("cs_chol: not nil")
	}
	if r := cs_cholsol(100, nil, nil); r == true {
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
	if _, err := cs_droptol(nil, -1); err == nil {
		t.Errorf("cs_droptol: not nil")
	}
	if _, err := cs_dropzeros(nil); err == nil {
		t.Errorf("cs_dropzeros: not nil")
	}
	// if err := Dupl(nil); err == nil {
	// 	t.Errorf("cs_dupl: not nil")
	// }
	// if r := Entry(nil, -1, -1, 0); r == true {
	// 	t.Errorf("cs_entry: not nil")
	// }
	if r := cs_ereach(nil, -1, nil, nil, nil); r == 1 {
		t.Errorf("cs_ereach: not nil")
	}
	if r := cs_etree(nil, false); r != nil {
		t.Errorf("cs_etree: not nil")
	}
	if _, err := Fkeep(nil, nil); err == nil {
		t.Errorf("cs_fkeep: not nil")
	}
	// if r := Gaxpy(nil, nil, nil); r == true {
	// 	t.Errorf("cs_gaxpy: not nil")
	// }
	if r := cs_happly(nil, -1, -1, nil); r == nil {
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
	// if err := Load(nil); r != nil {
	// 	t.Errorf("cs_load: not nil")
	// }
	if r := cs_lsolve(nil, nil); r != false {
		t.Errorf("cs_lsolve: not nil")
	}
	if r := cs_ltsolve(nil, nil); r == true {
		t.Errorf("cs_ltsolve: not nil")
	}
	if r := cs_lu(nil, nil, -1); r != nil {
		t.Errorf("cs_lu: not nil")
	}
	if r := cs_lusol(100, nil, nil, -1); r == true {
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
	// if _, r := Multiply(nil, nil); r == nil {
	// 	t.Errorf("cs_multiply: not nil")
	// }
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
	a := new(Matrix)
	if err := a.Print(nil, false); err == nil {
		t.Errorf("cs_print: not nil")
	}
	b := new(Matrix)
	if err := b.Print(nil, false); err == nil {
		t.Errorf("cs_print: not nil")
	}
	if r := cs_pvec(nil, nil, nil, -1); r == true {
		t.Errorf("cs_pvec: not nil")
	}
	if r := cs_qr(nil, nil); r != nil {
		t.Errorf("cs_qr: not nil")
	}
	if r := cs_qrsol(100, nil, nil); r == true {
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
	if r := cs_schol(100, nil); r != nil {
		t.Errorf("cs_schol: not nil")
	}
	if r := cs_spsolve(nil, nil, -1, nil, nil, nil, false); r != -1 {
		t.Errorf("cs_spsolve: not nil")
	}
	// TODO
	// if r := cs_vcount(nil, nil); r == true {
	// 	t.Errorf("cs_vcount: not nil")
	// }
	if r := cs_sqr(100, nil, false); r != nil {
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

const (
	MinUint uint = 0 // binary: all zeroes

	// Perform a bitwise NOT to change every bit from 0 to 1
	MaxUint = ^MinUint // binary: all ones

	// Shift the binary number to the right (i.e. divide by two)
	// to change the high bit to 0
	MaxInt = int(MaxUint >> 1) // binary: all ones except high bit

	// Perform another bitwise NOT to change the high bit to 1 and
	// all other bits to 0
	MinInt = ^MaxInt // binary: all zeroes except high bit
)

func TestCsCompress(t *testing.T) {

	t.Run("BigMatrix", func(t *testing.T) {
		T, err := NewTriplet()
		if err != nil {
			t.Fatal(err)
		}
		err = Entry(T, 0, 0, 1)
		if err != nil {
			t.Fatal(err)
		}
		err = Entry(T, MaxInt-2, MaxInt-2, 1)
		if err != nil {
			t.Fatal(err)
		}
		_, err = Compress(T)
		if err == nil {
			t.Fatal("Haven`t overflow error")
		}
		t.Log(err)
	})

	matrixes, err := filepath.Glob("CSparse/Matrix/" + "*")
	if err != nil {
		t.Fatal(err)
	}

	f := func(A *Matrix) string {
		tmpfile, err := ioutil.TempFile("", "cs_compress")
		if err != nil {
			panic(err)
		}

		// sort
		var ss []string
		for p := 0; p < A.nz; p++ {
			s := fmt.Sprintf("%8d %8d %10e", A.i[p], A.p[p], A.x[p])
			ss = append(ss, s)
		}
		sort.Strings(ss)

		// print
		for i := range ss {
			fmt.Fprintf(tmpfile, ss[i]+"\n")
		}

		filename := tmpfile.Name()
		defer func() { _ = os.Remove(filename) }()
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
			T, err := Load(&stdin)
			if err != nil {
				t.Fatal(err)
			}
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
				T2, err := Load(&buf)
				if err != nil {
					t.Fatalf("T2 is nil: %v", err)
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
				T2, err := Load(&buf)
				if err != nil {
					t.Fatalf("T2 is nil: %v", err)
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

func snapshot(filename string, t *testing.T, f func(tmpfile io.Writer)) {
	b, err := ioutil.ReadFile(filename)
	if err != nil {
		t.Fatal(err)
	}

	tmpfile, err := ioutil.TempFile("", "example")
	if err != nil {
		t.Fatal(err)
	}

	f(tmpfile)

	file := tmpfile.Name()
	defer func() { _ = os.Remove(file) }()
	err = tmpfile.Close()
	if err != nil {
		t.Fatal(err)
	}
	b2, err := ioutil.ReadFile(file)
	if err != nil {
		t.Fatal(err)
	}

	if !bytes.Equal(b, b2) {
		t.Logf("%s", ShowDiff(string(b), string(b2)))
		t.Fatalf("Results is not same:\n`%s`\n`%s`", string(b), string(b2))
	}
	t.Logf("`%s`", string(b))
}

func TestGaxpy(t *testing.T) {
	snapshot("./testdata/.snapshot.gaxpy", t, func(out io.Writer) {
		var s bytes.Buffer
		s.WriteString("0 0 1\n1 0 3\n2 0 5\n0 1 2\n1 1 4\n2 1 6")
		T, err := Load(&s)
		if err != nil {
			t.Fatal(err)
		}
		A, err := Compress(T)
		if err != nil {
			t.Fatal(err)
		}
		x := []float64{7, 8}
		y := []float64{9, 10, 11}
		err = Gaxpy(A, x, y)
		fmt.Fprintf(out, "%v\n", y)
		if err != nil {
			t.Fatal(err)
		}
	})
}

func ExampleGaxpy() {
	var s bytes.Buffer
	s.WriteString("0 0 1\n1 0 3\n2 0 5\n0 1 2\n1 1 4\n2 1 6")
	T, err := Load(&s)
	if err != nil {
		panic(err)
	}
	A, err := Compress(T)
	if err != nil {
		panic(err)
	}
	x := []float64{7, 8}
	y := []float64{9, 10, 11}
	fmt.Fprintln(os.Stdout, "Vector `y` before:")
	fmt.Fprintln(os.Stdout, y)
	err = Gaxpy(A, x, y)
	if err != nil {
		panic(err)
	}
	fmt.Fprintln(os.Stdout, "Vector `y` before:")
	fmt.Fprintln(os.Stdout, y)

	// Output:
	// Vector `y` before:
	// [9 10 11]
	// Vector `y` before:
	// [32 63 94]
}

func TestAdd(t *testing.T) {
	snapshot("./testdata/.snapshot.add", t, func(out io.Writer) {
		var stdin bytes.Buffer
		stdin.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
		T, err := Load(&stdin)
		if err != nil {
			t.Fatal(err)
		}
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
		R.Print(out, false)
	})
}

func ExampleAdd() {
	var stdin bytes.Buffer
	stdin.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4")
	T, err := Load(&stdin)
	if err != nil {
		panic(err)
	}
	A, err := Compress(T)
	if err != nil {
		panic(err)
	}
	AT, err := Transpose(A)
	if err != nil {
		panic(err)
	}
	R, err := Add(A, AT, 1, 2)
	if err != nil {
		panic(err)
	}
	R.Print(os.Stdout, false)

	// Output:
	// Sparse
	// 2-by-2, nzmax: 4 nnz: 4, 1-norm: 2.000000e+01
	//     col 0 : locations 0 to 1
	//       0 : 3.000000e+00
	//       1 : 7.000000e+00
	//     col 1 : locations 2 to 3
	//       0 : 8.000000e+00
	//       1 : 1.200000e+01
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
		c := []int{MaxInt, MaxInt}
		_, err := cs_cumsum(p, c)
		if err == nil {
			t.Fatalf("Error for overflow is not happen")
		}
		t.Log(err)
	})
	t.Run("overflow check 2", func(t *testing.T) {
		p := []int{0, 0, 0, 0}
		c := []int{MaxInt, 1, 5}
		_, err := cs_cumsum(p, c)
		if err == nil {
			t.Fatalf("Error for overflow is not happen")
		}
		t.Log(err)
	})
}

func TestDupl(t *testing.T) {
	snapshot("./testdata/.snapshot.dupl", t, func(out io.Writer) {
		var stdin bytes.Buffer
		stdin.WriteString("0 0 1\n0 1 2\n1 0 3\n1 1 4\n 0 0 1\n 1 0 10")
		T, err := Load(&stdin)
		if err != nil {
			t.Fatal(err)
		}
		A, err := Compress(T)
		if err != nil {
			t.Fatal(err)
		}
		A.Print(out, false)
		err = Dupl(A)
		if err != nil {
			t.Fatal(err)
		}
		A.Print(out, false)
	})
	// invert data
	snapshot("./testdata/.snapshot.dupl.invert", t, func(out io.Writer) {
		var stdin bytes.Buffer
		stdin.WriteString(" 1 0 10\n 0 0 1\n 1 1 4\n 1 0 3\n 0 1 2\n 0 0 1 ")
		T, err := Load(&stdin)
		if err != nil {
			t.Fatal(err)
		}
		A, err := Compress(T)
		if err != nil {
			t.Fatal(err)
		}
		A.Print(out, false)
		err = Dupl(A)
		if err != nil {
			t.Fatal(err)
		}
		A.Print(out, false)
	})
}

func ExampleDupl() {
	var stdin bytes.Buffer
	stdin.WriteString(" 1 0 10\n 0 0 1\n 1 1 4\n 1 0 3\n 0 1 2\n 0 0 1 ")
	T, err := Load(&stdin)
	if err != nil {
		panic(err)
	}
	A, err := Compress(T)
	if err != nil {
		panic(err)
	}
	fmt.Fprintln(os.Stdout, "Before:")
	A.Print(os.Stdout, false)
	err = Dupl(A)
	if err != nil {
		panic(err)
	}
	fmt.Fprintln(os.Stdout, "After:")
	A.Print(os.Stdout, false)

	// Output:
	// Before:
	// Sparse
	// 2-by-2, nzmax: 6 nnz: 6, 1-norm: 1.500000e+01
	//     col 0 : locations 0 to 3
	//       1 : 1.000000e+01
	//       0 : 1.000000e+00
	//       1 : 3.000000e+00
	//       0 : 1.000000e+00
	//     col 1 : locations 4 to 5
	//       1 : 4.000000e+00
	//       0 : 2.000000e+00
	// After:
	// Sparse
	// 2-by-2, nzmax: 4 nnz: 4, 1-norm: 1.500000e+01
	//     col 0 : locations 0 to 1
	//       1 : 1.300000e+01
	//       0 : 2.000000e+00
	//     col 1 : locations 2 to 3
	//       1 : 4.000000e+00
	//       0 : 2.000000e+00
}

func ExamplePrint() {
	T, err := NewTriplet()
	if err != nil {
		panic(err)
	}
	for i := 0; i < 25; i++ {
		err = Entry(T, i, i, 10+float64(i))
		if err != nil {
			panic(err)
		}
	}

	fmt.Fprintln(os.Stdout, "Full print of triplets:")
	T.Print(os.Stdout, false)

	fmt.Fprintln(os.Stdout, "Short print of triplets:")
	T.Print(os.Stdout, true)

	A, err := Compress(T)
	if err != nil {
		panic(err)
	}

	fmt.Fprintln(os.Stdout, "Full print of CSC matrix:")
	A.Print(os.Stdout, false)

	fmt.Fprintln(os.Stdout, "Short print of CSC matrix:")
	A.Print(os.Stdout, true)

	// Output:
	// Full print of triplets:
	// Sparse
	// triplet: 25-by-25, nzmax: 32 nnz: 25
	//     0 0 : 1.000000e+01
	//     1 1 : 1.100000e+01
	//     2 2 : 1.200000e+01
	//     3 3 : 1.300000e+01
	//     4 4 : 1.400000e+01
	//     5 5 : 1.500000e+01
	//     6 6 : 1.600000e+01
	//     7 7 : 1.700000e+01
	//     8 8 : 1.800000e+01
	//     9 9 : 1.900000e+01
	//     10 10 : 2.000000e+01
	//     11 11 : 2.100000e+01
	//     12 12 : 2.200000e+01
	//     13 13 : 2.300000e+01
	//     14 14 : 2.400000e+01
	//     15 15 : 2.500000e+01
	//     16 16 : 2.600000e+01
	//     17 17 : 2.700000e+01
	//     18 18 : 2.800000e+01
	//     19 19 : 2.900000e+01
	//     20 20 : 3.000000e+01
	//     21 21 : 3.100000e+01
	//     22 22 : 3.200000e+01
	//     23 23 : 3.300000e+01
	//     24 24 : 3.400000e+01
	// Short print of triplets:
	// Sparse
	// triplet: 25-by-25, nzmax: 32 nnz: 25
	//     0 0 : 1.000000e+01
	//     1 1 : 1.100000e+01
	//     2 2 : 1.200000e+01
	//     3 3 : 1.300000e+01
	//     4 4 : 1.400000e+01
	//     5 5 : 1.500000e+01
	//     6 6 : 1.600000e+01
	//     7 7 : 1.700000e+01
	//     8 8 : 1.800000e+01
	//     9 9 : 1.900000e+01
	//     10 10 : 2.000000e+01
	//     11 11 : 2.100000e+01
	//     12 12 : 2.200000e+01
	//     13 13 : 2.300000e+01
	//     14 14 : 2.400000e+01
	//     15 15 : 2.500000e+01
	//     16 16 : 2.600000e+01
	//     17 17 : 2.700000e+01
	//     18 18 : 2.800000e+01
	//     19 19 : 2.900000e+01
	//     20 20 : 3.000000e+01
	//     21 21 : 3.100000e+01
	//   ...
	// Full print of CSC matrix:
	// Sparse
	// 25-by-25, nzmax: 25 nnz: 25, 1-norm: 3.400000e+01
	//     col 0 : locations 0 to 0
	//       0 : 1.000000e+01
	//     col 1 : locations 1 to 1
	//       1 : 1.100000e+01
	//     col 2 : locations 2 to 2
	//       2 : 1.200000e+01
	//     col 3 : locations 3 to 3
	//       3 : 1.300000e+01
	//     col 4 : locations 4 to 4
	//       4 : 1.400000e+01
	//     col 5 : locations 5 to 5
	//       5 : 1.500000e+01
	//     col 6 : locations 6 to 6
	//       6 : 1.600000e+01
	//     col 7 : locations 7 to 7
	//       7 : 1.700000e+01
	//     col 8 : locations 8 to 8
	//       8 : 1.800000e+01
	//     col 9 : locations 9 to 9
	//       9 : 1.900000e+01
	//     col 10 : locations 10 to 10
	//       10 : 2.000000e+01
	//     col 11 : locations 11 to 11
	//       11 : 2.100000e+01
	//     col 12 : locations 12 to 12
	//       12 : 2.200000e+01
	//     col 13 : locations 13 to 13
	//       13 : 2.300000e+01
	//     col 14 : locations 14 to 14
	//       14 : 2.400000e+01
	//     col 15 : locations 15 to 15
	//       15 : 2.500000e+01
	//     col 16 : locations 16 to 16
	//       16 : 2.600000e+01
	//     col 17 : locations 17 to 17
	//       17 : 2.700000e+01
	//     col 18 : locations 18 to 18
	//       18 : 2.800000e+01
	//     col 19 : locations 19 to 19
	//       19 : 2.900000e+01
	//     col 20 : locations 20 to 20
	//       20 : 3.000000e+01
	//     col 21 : locations 21 to 21
	//       21 : 3.100000e+01
	//     col 22 : locations 22 to 22
	//       22 : 3.200000e+01
	//     col 23 : locations 23 to 23
	//       23 : 3.300000e+01
	//     col 24 : locations 24 to 24
	//       24 : 3.400000e+01
	// Short print of CSC matrix:
	// Sparse
	// 25-by-25, nzmax: 25 nnz: 25, 1-norm: 3.400000e+01
	//     col 0 : locations 0 to 0
	//       0 : 1.000000e+01
	//     col 1 : locations 1 to 1
	//       1 : 1.100000e+01
	//     col 2 : locations 2 to 2
	//       2 : 1.200000e+01
	//     col 3 : locations 3 to 3
	//       3 : 1.300000e+01
	//     col 4 : locations 4 to 4
	//       4 : 1.400000e+01
	//     col 5 : locations 5 to 5
	//       5 : 1.500000e+01
	//     col 6 : locations 6 to 6
	//       6 : 1.600000e+01
	//     col 7 : locations 7 to 7
	//       7 : 1.700000e+01
	//     col 8 : locations 8 to 8
	//       8 : 1.800000e+01
	//     col 9 : locations 9 to 9
	//       9 : 1.900000e+01
	//     col 10 : locations 10 to 10
	//       10 : 2.000000e+01
	//     col 11 : locations 11 to 11
	//       11 : 2.100000e+01
	//     col 12 : locations 12 to 12
	//       12 : 2.200000e+01
	//     col 13 : locations 13 to 13
	//       13 : 2.300000e+01
	//     col 14 : locations 14 to 14
	//       14 : 2.400000e+01
	//     col 15 : locations 15 to 15
	//       15 : 2.500000e+01
	//     col 16 : locations 16 to 16
	//       16 : 2.600000e+01
	//     col 17 : locations 17 to 17
	//       17 : 2.700000e+01
	//     col 18 : locations 18 to 18
	//       18 : 2.800000e+01
	//     col 19 : locations 19 to 19
	//       19 : 2.900000e+01
	//     col 20 : locations 20 to 20
	//       20 : 3.000000e+01
	//     col 21 : locations 21 to 21
	//       21 : 3.100000e+01
	//   ...
}

func ExampleMultiply() {
	var stdin bytes.Buffer
	stdin.WriteString(" 1 0 10\n 0 0 1\n 1 1 4\n 1 0 3\n 0 1 2\n 0 0 1 ")
	T, err := Load(&stdin)
	if err != nil {
		panic(err)
	}
	A, err := Compress(T)
	if err != nil {
		panic(err)
	}
	AT, err := Transpose(A)
	if err != nil {
		panic(err)
	}
	M, err := Multiply(A, AT)
	if err != nil {
		panic(err)
	}
	M.Print(os.Stdout, false)

	// Output:
	// Sparse
	// 2-by-2, nzmax: 4 nnz: 4, 1-norm: 2.190000e+02
	//     col 0 : locations 0 to 1
	//       1 : 3.400000e+01
	//       0 : 8.000000e+00
	//     col 1 : locations 2 to 3
	//       1 : 1.850000e+02
	//       0 : 3.400000e+01
}

func TestAmd(t *testing.T) {
	for _, a := range []Order{AmdNatural, AmdChol, AmdLU, AmdQR} {
		t.Logf("%s", a)
	}
}
