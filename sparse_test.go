package sparse

import (
	"fmt"
	"os/exec"
	"testing"
)

func Test(t *testing.T) {

	// build testdata application
	// clang -ICSparse/Include/ ./CSparse/Source/*.c ./testdata/csparse_test.c -lm -o ./testdata/csparse_test
	cmd := exec.Command("clang",
		"-ICSparse/Include/",
		"./CSparse/Source/*.c",
		"./testdata/csparse_test.c",
		"-lm",
		"-o",
		"./testdata/csparse_test",
	)
	_, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("cmd.Run() failed with %s\n", err)
	}

	// data checking
	cmd = exec.Command(
		"./testdata/csparse_test",
		"<",
		"./CSparse/Matrix/t1",
	)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("cmd.Run() failed with %s\n", err)
	}

	fmt.Println(string(out))
}
