# mms

[![Coverage Status](https://coveralls.io/repos/github/Konstantin8105/mms/badge.svg?branch=master)](https://coveralls.io/github/Konstantin8105/mms?branch=master)
[![Build Status](https://travis-ci.org/Konstantin8105/mms.svg?branch=master)](https://travis-ci.org/Konstantin8105/mms)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/mms)](https://goreportcard.com/report/github.com/Konstantin8105/mms)
[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Konstantin8105/mms/blob/master/LICENSE)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/mms?status.svg)](https://godoc.org/github.com/Konstantin8105/mms)

Memory managenet for slice


```
go test -v -bench=. -benchmem -count=10 > bench.prof
benchstat bench.prof 

name       time/op
/Direct-5   192µs ±18%
/Cache-5    144µs ± 9%

name       alloc/op
/Direct-5  23.1kB ± 0%
/Cache-5   4.06kB ± 0%

name       allocs/op
/Direct-5     123 ± 0%
/Cache-5      123 ± 0%
```

```
go test -v -bench=. -benchmem

=== RUN   Test
--- PASS: Test (0.00s)
goos: linux
goarch: amd64
pkg: github.com/Konstantin8105/mms
Benchmark/Direct-5         	    8082	    200886 ns/op	   23077 B/op	     123 allocs/op
Benchmark/Cache-5          	   10000	    142438 ns/op	    4065 B/op	     123 allocs/op
PASS
ok  	github.com/Konstantin8105/mms	3.307s
```


```
go test -v -bench=. -benchmem -memprofile=mem.prof
go tool pprof mem.prof

(pprof) top10 -cum
Showing nodes accounting for 123.56MB, 99.60% of 124.06MB total
Dropped 5 nodes (cum <= 0.62MB)
      flat  flat%   sum%        cum   cum%
         0     0%     0%   121.06MB 97.58%  github.com/Konstantin8105/mms.Benchmark.func1.1
  105.06MB 84.68% 84.68%   105.06MB 84.68%  github.com/Konstantin8105/mms.(*Direct).Get
      16MB 12.90% 97.58%       16MB 12.90%  github.com/Konstantin8105/mms.(*Cache).Put
         0     0% 97.58%     2.50MB  2.02%  github.com/Konstantin8105/mms.Benchmark.func1
    2.50MB  2.02% 99.60%     2.50MB  2.02%  github.com/Konstantin8105/mms.getChan
         0     0% 99.60%     2.50MB  2.02%  testing.(*B).launch
         0     0% 99.60%     2.50MB  2.02%  testing.(*B).runN
(pprof) exit
```

```
go test -v -bench=. -benchmem -cpuprofile=cpu.prof
go tool pprof cpu.prof
(pprof) top10 -cum
Showing nodes accounting for 1.89s, 40.04% of 4.72s total
Dropped 71 nodes (cum <= 0.02s)
Showing top 10 nodes out of 84
      flat  flat%   sum%        cum   cum%
     0.10s  2.94%  2.94%      2.78s 81.76%  github.com/Konstantin8105/mms.Benchmark.func1.1
     0.05s  1.47%  4.41%      2.10s 61.76%  math/rand.Float64
     0.05s  1.47%  5.88%      2.05s 60.29%  math/rand.(*Rand).Float64
     0.01s  0.29%  6.18%         2s 58.82%  math/rand.(*Rand).Int63
     0.17s  5.00% 11.18%      1.99s 58.53%  math/rand.(*lockedSource).Int63
     0.37s 10.88% 22.06%      0.98s 28.82%  sync.(*Mutex).Lock
     0.44s 12.94% 35.00%      0.68s 20.00%  sync.(*Mutex).Unlock
     0.40s 11.76% 46.76%      0.61s 17.94%  sync.(*Mutex).lockSlow
         0     0% 46.76%      0.35s 10.29%  runtime.mcall
     0.03s  0.88% 47.65%      0.27s  7.94%  runtime.schedule
(pprof) exit
```

