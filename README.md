# sparse

[![Coverage Status](https://coveralls.io/repos/github/Konstantin8105/sparse/badge.svg?branch=master)](https://coveralls.io/github/Konstantin8105/sparse?branch=master)
[![Build Status](https://travis-ci.org/Konstantin8105/sparse.svg?branch=master)](https://travis-ci.org/Konstantin8105/sparse)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/sparse)](https://goreportcard.com/report/github.com/Konstantin8105/sparse)
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Konstantin8105/sparse/blob/master/LICENSE)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/sparse?status.svg)](https://godoc.org/github.com/Konstantin8105/sparse)

**This package based on program CSparse from [SuiteSparse 5.3.0](http://faculty.cse.tamu.edu/davis/SuiteSparse/)**



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
