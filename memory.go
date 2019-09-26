package sparse

import (
	"sync"
)

// TODO: check data race

// TODO: change is actual. memory less ~5 times
//
// FROM:
//
// goos: linux
// goarch: amd64
// pkg: github.com/Konstantin8105/sparse
// BenchmarkLU/Sparse:__30:__300-5         	   41133	     30024 ns/op	   57120 B/op	      33 allocs/op
// BenchmarkLU/Sparse:_100:_1070-5         	   10000	    112641 ns/op	  221328 B/op	      33 allocs/op
// BenchmarkLU/Sparse:_300:_3270-5         	    3720	    313935 ns/op	  623760 B/op	      33 allocs/op
// BenchmarkLU/Sparse:1000:10970-5         	    1047	   1086014 ns/op	 1983378 B/op	      33 allocs/op
// BenchmarkLU/Sparse:3000:32970-5         	     378	   3110375 ns/op	 5948316 B/op	      33 allocs/op
//
// TO:
//
// goos: linux
// goarch: amd64
// pkg: github.com/Konstantin8105/sparse
// BenchmarkLU/Sparse:__30:__300-5         	   38524	     30243 ns/op	   10705 B/op	      37 allocs/op
// BenchmarkLU/Sparse:_100:_1070-5         	   12150	     98506 ns/op	   37533 B/op	      37 allocs/op
// BenchmarkLU/Sparse:_300:_3270-5         	    4090	    289935 ns/op	  108049 B/op	      37 allocs/op
// BenchmarkLU/Sparse:1000:10970-5         	    1082	   1089524 ns/op	  380083 B/op	      40 allocs/op
// BenchmarkLU/Sparse:3000:32970-5         	     399	   2951034 ns/op	 1102933 B/op	      40 allocs/op

type floatCache struct {
	// TODO : try to use RWMutex
	sync.Mutex
	caps  []int
	pools []*sync.Pool
}

var floats = floatCache{pools: []*sync.Pool{}}

func (f *floatCache) get(cap int) []float64 {
	f.Lock()
	defer func() {
		f.Unlock()
	}()

	// finding index
	index := -1
	for i := range f.caps {
		if f.caps[i] == cap {
			index = i
		}
	}

	if index < 0 {
		// create a new
		f.caps = append(f.caps, cap)
		f.pools = append(f.pools, &sync.Pool{
			New: func() interface{} {
				return make([]float64, cap)
			},
		})
		// TODO: typical list of caps
		// [1 2 4 8 16 32 64 128 256 512 300 31 1230 165 275 1024 2048 1070
		//  101 4380 585 1045 4096 3270 301 13380 1785 3245 8192 16384 10970
		//  1001 44880 5985 10945 32768 65536 32970 3001 134880 17985 32945]
		// TODO: add sorting and finding near size for minimaze Pools
		// TODO: add minimal capacity
		index = len(f.pools) - 1
		return f.pools[index].Get().([]float64)
	}
	// reuse
	arr := f.pools[index].Get().([]float64)
	for i := 0; i < cap; i++ {
		arr[i] = 0.0
	}
	return arr
}

func (f *floatCache) put(fs []float64) {
	f.Lock()
	defer func() {
		f.Unlock()
	}()

	cap := len(fs)

	// finding index
	index := -1
	for i := range f.caps {
		if f.caps[i] == cap {
			index = i
		}
	}

	if index < 0 {
		return
	}

	f.pools[index].Put(fs)
}

// ===========================================================================

// TODO: add golang template
type intCache struct {
	sync.Mutex
	caps  []int
	pools []*sync.Pool
}

var ints = intCache{pools: []*sync.Pool{}}

func (f *intCache) get(cap int) []int {
	f.Lock()
	defer func() {
		f.Unlock()
	}()

	// finding index
	index := -1
	for i := range f.caps {
		if f.caps[i] == cap {
			index = i
		}
	}

	if index < 0 {
		// create a new
		f.caps = append(f.caps, cap)
		f.pools = append(f.pools, &sync.Pool{
			New: func() interface{} {
				return make([]int, cap)
			},
		})
		index = len(f.pools) - 1
		return f.pools[index].Get().([]int)
	}
	// reuse
	arr := f.pools[index].Get().([]int)
	for i := 0; i < cap; i++ {
		arr[i] = 0
	}
	return arr
}

func (f *intCache) put(fs []int) {
	f.Lock()
	defer func() {
		f.Unlock()
	}()

	cap := len(fs)

	// finding index
	index := -1
	for i := range f.caps {
		if f.caps[i] == cap {
			index = i
		}
	}

	if index < 0 {
		return
	}

	f.pools[index].Put(fs)
}
