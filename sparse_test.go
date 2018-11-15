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

func Test(t *testing.T) {

	t.Run("Build test", func(t *testing.T) {
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
		args = append(args, "testdata/csparse_test.c")
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

		t.Run(matrixes[i], func(t *testing.T) {
			// data checking
			cmd := exec.Command(
				"./testdata/csparse_test",
			)

			var stdin, stdout, stderr bytes.Buffer
			b, err := ioutil.ReadFile(matrixes[i])
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
			fmt.Println(stdout.String())

			fmt.Println("-------")

			stdin.Write(b)
			T := cs_load(&stdin)
			T.cs_print(false)

			AC := T.cs_compress()
			AC.cs_print(false)
		})
	}
}
