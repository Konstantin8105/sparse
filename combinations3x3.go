// +build ignore

package main

import (
	"fmt"
	"strings"
)

// generate test matrixes
func main() {
	ind := node([]int{}, 3, 2)
	// fmt.Println("matrix positions: ", len(ind))
	vals := node([]int{}, 3, 3*3)
	// fmt.Println("amount combinations : ", len(vals))
	// fmt.Println("amount in each combination : ", len(vals[0]))
	str := []string{"-1", "0", "1"}

	for i := range vals {
		if len(vals[0]) != len(ind) {
			panic("sizes")
		}
		var out string
		for j := 0; j < len(ind); j++ {
			out += fmt.Sprintf("%1d %2s", ind[j], str[vals[i][j]])
			if j != len(ind)-1 {
				out += "\\n"
			}
		}
		out = strings.Replace(out, "[", " ", -1)
		out = strings.Replace(out, "]", " ", -1)
		fmt.Printf("%s\n", out)
	}
}

func node(last []int, vals int, deep int) (result [][]int) {
	if deep == 0 {
		result = append(result, last)
		return
	}

	deep--

	for i := 0; i < vals; i++ {
		cp := make([]int, len(last)+1)
		copy(cp, last)
		cp[len(last)] = i
		result = append(result, node(cp, vals, deep)...)
	}

	return
}
