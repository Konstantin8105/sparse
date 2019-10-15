package mms

import (
	"fmt"
	"runtime"
)

// Debug is flag for switching to debug mode
var Debug bool = false

func called() (info string) {
	for i := 0; ; i++ {
		function, file, line, ok := runtime.Caller(i)
		if !ok {
			break
		}
		info += fmt.Sprintf(
			"\tFile     : %s\n"+
				"\tFunction : %s\n"+
				"\tLine     : %d\n"+
				"\n",
			file,
			runtime.FuncForPC(function).Name(),
			line,
		)
	}
	return
}
