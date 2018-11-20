# sparse

[![Coverage Status](https://coveralls.io/repos/github/Konstantin8105/sparse/badge.svg?branch=master)](https://coveralls.io/github/Konstantin8105/sparse?branch=master)
[![Build Status](https://travis-ci.org/Konstantin8105/sparse.svg?branch=master)](https://travis-ci.org/Konstantin8105/sparse)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/sparse)](https://goreportcard.com/report/github.com/Konstantin8105/sparse)
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Konstantin8105/sparse/blob/master/LICENSE)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/sparse?status.svg)](https://godoc.org/github.com/Konstantin8105/sparse)


At the base on transpilated CSparse from C to Go by [c4go](https://github.com/Konstantin8105/c4go).



Result of:
```
└─▪ cat sparse.go | grep -v "func()" | grep -v "function" | grep "func" | wc -l
75

└─▪ cat sparse.go | grep -v "func()" | grep -v "function" | grep "func" |  grep -v "//" | wc -l
71
```

#### **CSparse base on SuiteSparse 5.3.0**

http://faculty.cse.tamu.edu/davis/SuiteSparse/
