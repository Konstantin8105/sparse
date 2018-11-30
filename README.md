# sparse

[![Coverage Status](https://coveralls.io/repos/github/Konstantin8105/sparse/badge.svg?branch=master)](https://coveralls.io/github/Konstantin8105/sparse?branch=master)
[![Build Status](https://travis-ci.org/Konstantin8105/sparse.svg?branch=master)](https://travis-ci.org/Konstantin8105/sparse)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/sparse)](https://goreportcard.com/report/github.com/Konstantin8105/sparse)
[![GitHub license](https://img.shields.io/badge/license-LGPL%20v2.1-blue.svg)](https://github.com/Konstantin8105/sparse/blob/master/LICENSE)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/sparse?status.svg)](https://godoc.org/github.com/Konstantin8105/sparse)

**This package based on program CSparse from [SuiteSparse 5.3.0](http://faculty.cse.tamu.edu/davis/SuiteSparse/)**

### Example of comparing CSparse and CXSparse

```cmd
rm -rf /tmp/CXSparse
cp -R ./CXSparse/ /tmp/
sed -i 's/CS_ENTRY/csi/g' /tmp/CXSparse/Source/*.c
sed -i 's/CS_INT/csi/g'   /tmp/CXSparse/Source/*.c
meld  /tmp/CXSparse/ ./CSparse/
```

### How updating package

* Check new version of `CSparse` is exist on [page](http://faculty.cse.tamu.edu/davis/SuiteSparse/).
* Download new version.
* Compare file-by-file with file in folder `CSparse` of that package.
* Inject changes into Go files of package.

> Note:
> CSparse haven`t any updates at the last few years, so
> that package is actual at the future.
>

### Just for information

**This package transpiled CSparse from C to Go by [c4go](https://github.com/Konstantin8105/c4go).**

### Profiling

```
go test -v -cpuprofile cpu.prof -memprofile mem.prof -run=Benchmark -bench=Benchmark -benchmem
go tool pprof cpu.prof
go tool pprof mem.prof
```

### Questions for CSparse

* Variables `css.lnz` and `css.unz` is not float type `double`, better to use integer type like `int`.
* Not clear size of variable `w` in function `cs_dupl`. By default it is take - amount of rows, but we can take less memory, if we will use next calculation. But we have to reindex :
```go
	// find maximal amount of non-zero values in row
	m := 0
	for j := 0; j < n; j++ {      // iteration by columns in CSC
		if mz < Ap[j+1]-Ap[j] {   // comparing
			mz = Ap[j+1] - Ap[j]  // save new size
		}
	}
	fmt.Printf("Max. non-zero values in row: %d \n", mz)
	for i := 0; i < mz; i++ {
		w[i] = -1
	}
	// TODO: need reindex
	...
	i := Ai[p]
	if w[i] >= g{ // need reindex w[i]
	...
```
* 
