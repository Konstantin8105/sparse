// Code generated by mmsgen; DO NOT EDIT.

package mms

import (
	"fmt"
	"sort"
	"sync"
)

// Int32sCache of slices
type Int32sCache struct {
	mutex  sync.RWMutex
	ps     []poolInt32sCache
	putarr []debugInt32sCache
}

type poolInt32sCache struct {
	p    *sync.Pool
	size int
}

type debugInt32sCache struct {
	arr  *[]int32
	line string
}

// Get return slice
func (c *Int32sCache) Get(size int) []int32 {

	// lock
	c.mutex.Lock()
	defer func() {
		c.mutex.Unlock()
	}()

	// finding index
	index := c.index(size)

	// creating a new pool
	if index < 0 {
		c.ps = append(c.ps, poolInt32sCache{
			p: &sync.Pool{
				New: func() interface{} {
					return make([]int32, size)
				},
			},
			size: size,
		})
		// sort
		sort.SliceStable(c.ps, func(i, j int) bool {
			return c.ps[i].size < c.ps[j].size
		})

		// return
		index = c.index(size)
		return c.ps[index].p.New().([]int32)
	}

	// pool is found
	arr := c.ps[index].p.Get().([]int32)

	if len(arr) == 0 {
		arr = arr[:size]
	}

	if Debug {
		if len(arr) < size {
			panic(fmt.Errorf("not same sizes: %d != %d", len(arr), size))
		}
		if len(arr) != cap(arr) {
			panic(fmt.Errorf("not valid capacity: %d != %d", len(arr), cap(arr)))
		}
	}

	for i := range arr {
		// initialization of slice
		arr[i] = 0
	}
	return arr
}

// Put slice into pool
func (c *Int32sCache) Put(arr *[]int32) {
	c.mutex.RLock() // lock
	var (
		size  = cap(*arr)
		index = c.index(size) // finding index
	)
	c.mutex.RUnlock() // unlock

	if index < 0 {
		// pool is not exist
		return
	}
	if size == 0 {
		// empty size
		return
	}
	if len(*arr) == 0 {
		// propably it is a dublicate putting
		return
	}

	// lock and add
	c.mutex.Lock()
	defer func() {
		c.mutex.Unlock()
	}()
	if !(index < len(c.ps) && c.ps[index].size == size) {
		return
	}
	*arr = (*arr)[:0]
	if Debug {
		// check if putting same arr
		for i := range c.putarr {
			if c.putarr[i].arr != arr {
				continue
			}
			panic(fmt.Errorf(
				"Dublicate of Put.\n"+
					"Last is called in :\n%v\n"+
					"Present call in   :\n%v\n",
				c.putarr[i].line,
				called(),
			))
		}
		c.putarr = append(c.putarr, debugInt32sCache{
			arr:  arr,
			line: called(),
		})
		return
	}
	c.ps[index].p.Put(*arr)
}

// return index with excepted size
func (c *Int32sCache) index(size int) int {
	index := -1
	for i := range c.ps {
		if c.ps[i].size < size {
			continue
		}
		if c.ps[i].size == size {
			index = i
		}
		break
	}
	return index
}

// Reset internal structure.
// In debug case - better for founding double putting.
// In normal case - for memory management with different memory profile.
//
//	Example of code:
//	w := cache.Get(10)
//	defer func() {
//		if mms.Debug {
//			cache.Reset()
//		}
//	}
//	... // Put memory in cache in next lines of code
//
func (c *Int32sCache) Reset() {
	// lock
	c.mutex.Lock()
	defer func() {
		c.mutex.Unlock()
	}()

	// remove
	c.ps = make([]poolInt32sCache, 0)
	c.putarr = make([]debugInt32sCache, 0)
}