package sparse

import (
	"fmt"
	"io"
	"math"

	"golang.org/x/exp/rand"
)

// matrix in compressed-column or triplet form
type cs struct { // struct cs_sparse
	nzmax int       // maximum number of entries
	m     int       // number of rows
	n     int       // number of columns
	p     []int     // column pointers (size n+1) or col indices (size nzmax)
	i     []int     // row indices, size nzmax
	x     []float64 // numerical values, size nzmax
	nz    int       // # of entries in triplet matrix, -1 for compressed-col
}

// symbolic Cholesky, LU, or QR analysis
type css struct { // struct cs_symbolic
	pinv     *int    // inverse row perm. for QR, fill red. perm for Chol
	q        *int    // fill-reducing column permutation for LU and QR
	parent   *int    // elimination tree for Cholesky and QR
	cp       *int    // column pointers for Cholesky, row counts for QR
	leftmost *int    // TODO
	m2       int     // # of rows for QR, after adding fictitious rows
	lnz      float64 // # entries in L for LU or Cholesky; in V for QR
	unz      float64 // # entries in U for LU; in R for QR
}

// numeric Cholesky, LU, or QR factorization
type csn struct { // struct cs_numeric
	L    *cs      // L for LU and Cholesky, V for QR
	U    *cs      // U for LU, R for QR, not used for Cholesky
	pinv *int     // partial pivoting for LU
	B    *float64 // beta [0..n-1] for QR
}

// cs_dmperm or cs_scc output
type csd struct { // struct cs_dmperm_results
	p  []int  // size m, row permutation
	q  []int  // size n, column permutation
	r  []int  // size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q)
	s  []int  // size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q)
	nb int    // # of blocks in fine dmperm decomposition
	rr [5]int // coarse row decomposition
	cc [5]int // coarse column decomposition
}

// TODO(KI) : remove comments like "transpiled function from ..."

// cs_add - C = alpha*A + beta*B
func cs_add(A *cs, B *cs, alpha float64, beta float64) *cs {
	var x []float64
	if !(A != nil && A.nz == -1) || !(B != nil && B.nz == -1) {
		// check inputs
		return nil
	}
	if (A.m != B.m) || (A.n != B.n) {
		return nil
	}
	var (
		m   = A.m
		anz = A.p[A.n]
		n   = B.n
		Bp  = B.p
		Bx  = B.x
		bnz = Bp[n]
		// get workspace
		w = make([]int, m)
	)
	values := (A.x != nil && Bx != nil)
	// get workspace
	if values {
		x = make([]float64, m)
	}
	// allocate result
	C := cs_spalloc(m, n, anz+bnz, values, false)
	if C == nil || w == nil || values && x == nil {
		return cs_done(C, w, x, false)
	}
	var (
		Cp = C.p
		Ci = C.i
		Cx = C.x
	)
	var nz int
	for j := 0; j < n; j++ {
		// column j of C starts here
		Cp[j] = nz
		// alpha*A(:,j)
		nz = cs_scatter(A, j, alpha, w, x, j+1, C, nz)
		// beta*B(:,j)
		nz = cs_scatter(B, j, beta, w, x, j+1, C, nz)
		if values {
			for p := Cp[j]; p < nz; p++ {
				Cx[p] = x[Ci[p]]
			}
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	// remove extra space from C
	cs_sprealloc(C, 0)
	// success; free workspace, return C
	return cs_done(C, w, x, true)
}

//
// // cs_wclear - clear w
// func cs_wclear(mark, lemax int, w []int, n int) int {
// 	if mark < 2 || mark+lemax < 0 {
// 		for k := 0; k < n; k++ {
// 			if w[k] != 0 {
// 				w[k] = 1
// 			}
// 		}
// 		mark = 2
// 	}
// 	// at this point, w [0..n-1] < mark holds
// 	return mark
// }
//
// // cs_diag - keep off-diagonal entries; drop diagonal entries
// func cs_diag(i, j int, aij float64, other interface{}) bool {
//
// 	// TODO (KI) : remove arguments aij, other
//
// 	return (i != j)
// }

const (
	Natural int = iota
	Chol
	LU
	QR
)

// // cs_amd - p = amd(A+A') if symmetric is true, or amd(A'A) otherwise
// // order 0:natural, 1:Chol, 2:LU, 3:QR
// func cs_amd(order Order, A *cs) *cs {
// 	var C *cs
// 	var A2 *cs
// 	// var AT []cs
// 	// var Cp []int
// 	// var Ci []noarch.PtrdiffT
// 	// var last []noarch.PtrdiffT
// 	// var W []noarch.PtrdiffT
// 	// var len []noarch.PtrdiffT
// 	// var nv []noarch.PtrdiffT
// 	// var next []noarch.PtrdiffT
// 	// var P []noarch.PtrdiffT
// 	// var head []noarch.PtrdiffT
// 	// var elen []noarch.PtrdiffT
// 	// var degree []noarch.PtrdiffT
// 	// var w []noarch.PtrdiffT
// 	// var hhead []noarch.PtrdiffT
// 	var ATp []int
// 	var ATi []int
// 	var d int     // noarch.PtrdiffT
// 	var dk int    // noarch.PtrdiffT
// 	var dext int  // noarch.PtrdiffT
// 	var lemax int // noarch.PtrdiffT
// 	var e int     // noarch.PtrdiffT
// 	// var elenk noarch.PtrdiffT
// 	var eln int // noarch.PtrdiffT
// 	var i int   // noarch.PtrdiffT
// 	var j int
// 	// var k noarch.PtrdiffT
// 	var k1 int // noarch.PtrdiffT
// 	var k2 int // noarch.PtrdiffT
// 	// var k3 noarch.PtrdiffT
// 	var jlast int // noarch.PtrdiffT
// 	var ln int    // noarch.PtrdiffT
// 	// var dense noarch.PtrdiffT
// 	// var nzmax noarch.PtrdiffT
// 	// var mindeg noarch.PtrdiffT
// 	var nvi int // noarch.PtrdiffT
// 	var nvj int // noarch.PtrdiffT
// 	// var nvk noarch.PtrdiffT
// 	var mark int //  noarch.PtrdiffT
// 	var wnvi int // noarch.PtrdiffT
// 	var ok bool  // noarch.PtrdiffT
// 	// var cnz noarch.PtrdiffT
// 	var nel int // noarch.PtrdiffT
// 	var p int
// 	var p1 int // noarch.PtrdiffT
// 	var p2 int
// 	var p3 int  // noarch.PtrdiffT
// 	var p4 int  // noarch.PtrdiffT
// 	var pj int  // noarch.PtrdiffT
// 	var pk int  // noarch.PtrdiffT
// 	var pk1 int // noarch.PtrdiffT
// 	var pk2 int // noarch.PtrdiffT
// 	var pn int  // noarch.PtrdiffT
// 	// var q noarch.PtrdiffT
// 	// var n noarch.PtrdiffT
// 	// var m noarch.PtrdiffT
// 	// var t noarch.PtrdiffT
// 	var h int // noarch.PtrdiffT
// 	if !(A != nil && A.nz == -1) || order <= 0 || order > 3 {
// 		// --- Construct matrix C -----------------------------------------------
// 		// check
// 		return nil
// 	}
// 	// compute A'
// 	AT := cs_transpose(A, false)
// 	if AT == nil {
// 		return nil
// 	}
// 	var (
// 		m = A.m
// 		n = A.n
// 	)
// 	// find dense threshold
// 	dense := 16
// 	if val := int(10.0 * math.Sqrt(float64(n))); dense < val {
// 		dense = val
// 	}
// 	if n-2 < dense {
// 		dense = n - 2
// 	}
//
// 	if order == 1 && n == m {
// 		// C = A+A'
// 		C = cs_add(A, AT, 0, 0)
// 	} else if order == 2 {
// 		// drop dense columns from AT
// 		ATp = AT.p
// 		ATi = AT.i
//
// 		for p2, j = 0, 0; j < m; j++ {
// 			// column j of AT starts here
// 			p = ATp[j]
// 			// new column j starts here
// 			ATp[j] = p2
// 			if ATp[j+1]-p > dense {
// 				// skip dense col j
// 				continue
// 			}
// 			for ; p < ATp[j+1]; p++ {
// 				ATi[p2] = ATi[p]
// 				p2++
// 			}
// 		}
//
// 		// finalize AT
// 		ATp[m] = p2
// 		// A2 = AT'
// 		A2 = cs_transpose(AT, false)
// 		// C=A'*A with no dense rows
// 		if A2 != nil {
// 			C = cs_multiply(AT, A2)
// 		} else {
// 			C = nil
// 		}
// 		cs_spfree(A2) // TODO (KI) : remove
// 	} else {
// 		// C=A'*A
// 		C = cs_multiply(AT, A)
// 	}
// 	cs_spfree(AT) // TODO (KI) : remove
// 	if C == nil {
// 		return nil
// 	}
// 	// drop diagonal entries
// 	cs_fkeep(C, cs_diag, nil)
// 	Cp := C.p
// 	cnz := Cp[n]
// 	// allocate result
// 	P := make([]int, n+1)
// 	// get workspace
// 	var W [8][]int
// 	for i := 0; i < 8; i++ {
// 		W[i] = make([]int, n+1)
// 	}
// 	// add elbow room to C
// 	t := cnz + cnz/5 + 2*n/8
// 	if P == nil || cs_sprealloc(C, t) { //  || W == nil
// 		return cs_idone(P, C, W, false)
// 	}
// 	var (
// 		len    = W[0]
// 		nv     = W[1*(n+1)] //(*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(n+1))*unsafe.Sizeof(W[0]))))[:]
// 		next   = W[2*(n+1)] // (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(2*int32(n+1)))*unsafe.Sizeof(W[0]))))[:]
// 		head   = W[3*(n+1)] // (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(3*int32(n+1)))*unsafe.Sizeof(W[0]))))[:]
// 		elen   = W[4*(n+1)] // (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(4*int32(n+1)))*unsafe.Sizeof(W[0]))))[:]
// 		degree = W[5*(n+1)] // (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(5*int32(n+1)))*unsafe.Sizeof(W[0]))))[:]
// 		w      = W[6*(n+1)] // (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(6*int32(n+1)))*unsafe.Sizeof(W[0]))))[:]
// 		hhead  = W[7*(n+1)] // (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(7*int32(n+1)))*unsafe.Sizeof(W[0]))))[:]
//
// 		// use P as workspace for last
// 		last = P
// 	)
//
// 	// --- Initialize quotient graph ----------------------------------------
// 	for k := 0; k < n; k++ {
// 		len[k] = Cp[k+1] - Cp[k]
// 	}
//
// 	len[n] = 0
// 	nzmax := C.nzmax
// 	Ci := C.i
//
// 	// type graph struct { // TODO (KI) : I think struct is look like that
// 	// 	head   int
// 	// 	last   int
// 	// 	hhead  int
// 	// 	nv     int
// 	// 	elen   int
// 	// 	degree int
// 	// }
//
// 	for i := 0; i <= n; i++ {
// 		// degree list i is empty
// 		head[i] = -1
// 		last[i] = -1
// 		next[i] = -1
// 		// hash list i is empty
// 		hhead[i] = -1
// 		// node i is just one node
// 		nv[i] = 1
// 		// node i is alive
// 		w[i] = 1
// 		// Ek of node i is empty
// 		elen[i] = 0
// 		// degree of node i
// 		degree[i] = len[i]
// 	}
// 	// clear w
// 	mark = cs_wclear(0, 0, w, n)
// 	// n is a dead element
// 	elen[n] = -2
// 	// n is a root of assembly tree
// 	Cp[n] = -1
// 	// n is a dead element
// 	w[n] = 0
//
// 	// --- Initialize degree lists ------------------------------------------
// 	for i := 0; i < n; i++ {
// 		d := degree[i]
// 		if d == 0 {
// 			// node i is empty
// 			// element i is dead
// 			elen[i] = -2
// 			nel++
// 			// i is a root of assembly tree
// 			Cp[i] = -1
// 			w[i] = 0
// 		} else if d > dense {
// 			// node i is dense
// 			// absorb i into element n
// 			nv[i] = 0
// 			// node i is dead
// 			elen[i] = -1
// 			nel++
// 			Cp[i] = -n - 2
// 			nv[n]++
// 		} else {
// 			if head[d] != -1 {
// 				last[head[d]] = i
// 			}
// 			// put node i in degree list d
// 			next[i] = head[d]
// 			head[d] = i
// 		}
// 	}
//
// 	var (
// 		mindeg int
// 		k      int
// 		elenk  int
// 		nvk    int
// 	)
//
// 	// while (selecting pivots) do
// 	for nel < n {
//
// 		// --- Select node of minimum approximate degree --------------------
// 		for k = -1; mindeg < n && (func() int {
// 			k = head[mindeg]
// 			return k
// 		}()) == -1; mindeg++ {
// 		}
//
// 		if next[k] != -1 {
// 			last[next[k]] = -1
// 		}
// 		// remove k from degree list
// 		head[mindeg] = next[k]
// 		// elenk = |Ek|
// 		elenk = elen[k]
// 		// # of nodes k represents
// 		nvk = nv[k]
// 		// nv[k] nodes of A eliminated
// 		nel += nvk
// 		if elenk > 0 && cnz+mindeg >= nzmax {
//
// 			// --- Garbage collection -------------------------------------------
// 			for j = 0; j < n; j++ {
// 				if (func() int {
// 					p = Cp[j]
// 					return p
// 				}()) >= 0 {
// 					// j is a live node or element
// 					// save first entry of object
// 					Cp[j] = Ci[p]
// 					// first entry is now CS_FLIP(j)
// 					Ci[p] = -j - 2
// 				}
// 			}
//
// 			// scan all of memory
// 			var (
// 				q, p int
// 			)
// 			for p = 0; p < cnz; {
// 				if (func() int {
// 					j = -(Ci[func() int {
// 						defer func() {
// 							p++
// 						}()
// 						return p
// 					}()]) - 2
// 					return j
// 				}()) >= 0 {
// 					// found object j
// 					// restore first entry of object
// 					Ci[q] = Cp[j]
// 					// new pointer to object j
// 					Cp[j] = q
// 					q++
// 					for k3 := 0; k3 < len[j]-1; k3++ {
// 						Ci[func() int {
// 							defer func() {
// 								q++
// 							}()
// 							return q
// 						}()] = Ci[func() int {
// 							defer func() {
// 								p++
// 							}()
// 							return p
// 						}()]
// 					}
// 				}
// 			}
//
// 			// Ci [cnz...nzmax-1] now free
// 			cnz = q
// 		}
// 		// --- Construct new element ----------------------------------------
// 		dk = 0
// 		// flag k as in Lk
// 		nv[k] = -nvk
// 		p = Cp[k]
// 		// do in place if elen[k] == 0
// 		pk1 = func() int {
// 			if elenk == 0 {
// 				return p
// 			}
// 			return cnz
// 		}()
// 		pk2 = pk1
// 		for k1 := 1; k1 <= elenk+1; k1++ {
// 			if k1 > elenk {
// 				// search the nodes in k
// 				e = k
// 				// list of nodes starts at Ci[pj]
// 				pj = p
// 				// length of list of nodes in k
// 				ln = len[k] - elenk
// 			} else {
// 				// search the nodes in e
// 				e = Ci[func() int {
// 					defer func() {
// 						p++
// 					}()
// 					return p
// 				}()]
// 				pj = Cp[e]
// 				// length of list of nodes in e
// 				ln = len[e]
// 			}
// 			for k2 := 1; k2 <= ln; k2++ {
// 				i := Ci[func() int {
// 					defer func() {
// 						pj++
// 					}()
// 					return pj
// 				}()]
// 				if (func() int {
// 					nvi = nv[i]
// 					return nvi
// 				}()) <= 0 {
// 					// node i dead, or seen
// 					continue
// 				}
// 				// degree[Lk] += size of node i
// 				dk += nvi
// 				// negate nv[i] to denote i in Lk
// 				nv[i] = -nvi
// 				// place i in Lk
// 				Ci[func() int {
// 					defer func() {
// 						pk2++
// 					}()
// 					return pk2
// 				}()] = i
// 				if next[i] != -1 {
// 					last[next[i]] = last[i]
// 				}
// 				if last[i] != -1 {
// 					// remove i from degree list
// 					next[last[i]] = next[i]
// 				} else {
// 					head[degree[i]] = next[i]
// 				}
// 			}
// 			if e != k {
// 				// absorb e into k
// 				Cp[e] = -k - 2
// 				// e is now a dead element
// 				w[e] = 0
// 			}
// 		}
// 		if elenk != 0 {
// 			// Ci [cnz...nzmax] is free
// 			cnz = pk2
// 		}
// 		// external degree of k - |Lk\i|
// 		degree[k] = dk
// 		// element k is in Ci[pk1..pk2-1]
// 		Cp[k] = pk1
// 		len[k] = pk2 - pk1
// 		// k is now an element
// 		elen[k] = -2
// 		// --- Find set differences -----------------------------------------
// 		// clear w if necessary
// 		mark = cs_wclear(mark, lemax, w, n)
//
// 		// scan 1: find |Le\Lk|
// 		for pk = pk1; pk < pk2; pk++ {
// 			i = Ci[pk]
// 			if (func() int {
// 				eln = elen[i]
// 				return eln
// 			}()) <= 0 {
// 				// skip if elen[i] empty
// 				continue
// 			}
// 			// nv [i] was negated
// 			nvi = -nv[i]
// 			wnvi = mark - nvi
// 			{
// 				// scan Ei
// 				for p = Cp[i]; p <= Cp[i]+eln-1; p++ {
// 					e = Ci[p]
// 					if w[e] >= mark {
// 						// decrement |Le\Lk|
// 						w[e] -= nvi
// 					} else if w[e] != 0 {
// 						// ensure e is a live element
// 						// 1st time e seen in scan 1
// 						w[e] = degree[e] + wnvi
// 					}
// 				}
// 			}
// 		}
//
// 		// --- Degree update ------------------------------------------------
// 		// scan2: degree update
// 		for pk = pk1; pk < pk2; pk++ {
// 			// consider node i in Lk
// 			i = Ci[pk]
// 			p1 = Cp[i]
// 			p2 = p1 + elen[i] - 1
// 			pn = p1
// 			{
// 				// scan Ei
// 				h = 0
// 				d = 0
// 				p = p1
// 				for p = p1; p <= p2; p++ {
// 					e = Ci[p]
// 					if w[e] != 0 {
// 						// e is an unabsorbed element
// 						// dext = |Le\Lk|
// 						dext = w[e] - mark
// 						if dext > 0 {
// 							// sum up the set differences
// 							d += dext
// 							// keep e in Ei
// 							Ci[func() int {
// 								defer func() {
// 									pn++
// 								}()
// 								return pn
// 							}()] = e
// 							// compute the hash of node i
// 							h += e
// 						} else {
// 							// aggressive absorb. e->k
// 							Cp[e] = -k - 2
// 							// e is a dead element
// 							w[e] = 0
// 						}
// 					}
// 				}
// 			}
// 			// elen[i] = |Ei|
// 			elen[i] = pn - p1 + 1
// 			p3 = pn
// 			p4 = p1 + len[i]
// 			{
// 				// prune edges in Ai
// 				for p = p2 + 1; p < p4; p++ {
// 					j = Ci[p]
// 					if (func() int {
// 						nvj = nv[j]
// 						return nvj
// 					}()) <= 0 {
// 						// node j dead or in Lk
// 						continue
// 					}
// 					// degree(i) += |j|
// 					d += nvj
// 					// place j in node list of i
// 					Ci[func() int {
// 						defer func() {
// 							pn++
// 						}()
// 						return pn
// 					}()] = j
// 					// compute hash for node i
// 					h += j
// 				}
// 			}
// 			if d == 0 {
// 				// check for mass elimination
// 				// absorb i into k
// 				Cp[i] = -k - 2
// 				nvi = -nv[i]
// 				// |Lk| -= |i|
// 				dk -= nvi
// 				// |k| += nv[i]
// 				nvk += nvi
// 				nel += nvi
// 				nv[i] = 0
// 				// node i is dead
// 				elen[i] = -1
// 			} else {
// 				// update degree(i)
// 				degree[i] = func() int {
// 					if degree[i] < d {
// 						return degree[i]
// 					}
// 					return d
// 				}()
// 				// move first node to end
// 				Ci[pn] = Ci[p3]
// 				// move 1st el. to end of Ei
// 				Ci[p3] = Ci[p1]
// 				// add k as 1st element in of Ei
// 				Ci[p1] = k
// 				// new len of adj. list of node i
// 				len[i] = pn - p1 + 1
// 				// finalize hash of i
// 				h = func() int {
// 					if h < 0 {
// 						return -h
// 					}
// 					return h
// 				}() % int(n)
// 				// place i in hash bucket
// 				next[i] = hhead[h]
// 				hhead[h] = i
// 				// save hash of i in last[i]
// 				last[i] = h
// 			}
// 		}
//
// 		// scan2 is done
// 		// finalize |Lk|
// 		degree[k] = dk
// 		lemax = func() int {
// 			if lemax > dk {
// 				return lemax
// 			}
// 			return dk
// 		}()
// 		// clear w
// 		mark = cs_wclear(mark+lemax, lemax, w, n)
//
// 		// --- Supernode detection ------------------------------------------
// 		for pk = pk1; pk < pk2; pk++ {
// 			i = Ci[pk]
// 			if nv[i] >= 0 {
// 				// skip if i is dead
// 				continue
// 			}
// 			// scan hash bucket of node i
// 			h = last[i]
// 			i = hhead[h]
// 			// hash bucket will be empty
// 			hhead[h] = -1
// 			for i != -1 && next[i] != -1 {
// 				ln = len[i]
// 				eln = elen[i]
// 				for p = Cp[i] + 1; p <= Cp[i]+ln-1; p++ {
// 					w[Ci[p]] = mark
// 				}
// 				jlast = i
// 				{
// 					// compare i with all j
// 					for j = next[i]; j != -1; {
// 						ok = (len[j] == ln && elen[j] == eln)
// 						for p = Cp[j] + 1; bool(ok) && p <= Cp[j]+ln-1; p++ {
// 							if w[Ci[p]] != mark {
// 								// compare i and j
// 								ok = false
// 							}
// 						}
// 						if bool(ok) {
// 							// i and j are identical
// 							// absorb j into i
// 							Cp[j] = -i - 2
// 							nv[i] += nv[j]
// 							nv[j] = 0
// 							// node j is dead
// 							elen[j] = -1
// 							// delete j from hash bucket
// 							j = next[j]
// 							next[jlast] = j
// 						} else {
// 							// j and i are different
// 							jlast = j
// 							j = next[j]
// 						}
// 					}
// 				}
// 				i = next[i]
// 				mark++
// 			}
// 		}
//
// 		// --- Finalize new element------------------------------------------
// 		// finalize Lk
// 		p = pk1
// 		pk = pk1
// 		for pk = pk1; pk < pk2; pk++ {
// 			i = Ci[pk]
// 			if (func() int {
// 				nvi = -nv[i]
// 				return nvi
// 			}()) <= 0 {
// 				// skip if i is dead
// 				continue
// 			}
// 			// restore nv[i]
// 			nv[i] = nvi
// 			// compute external degree(i)
// 			d = degree[i] + dk - nvi
// 			d = func() int {
// 				if d < n-nel-nvi {
// 					return d
// 				}
// 				return n - nel - nvi
// 			}()
// 			if head[d] != -1 {
// 				last[head[d]] = i
// 			}
// 			// put i back in degree list
// 			next[i] = head[d]
// 			last[i] = -1
// 			head[d] = i
// 			// find new minimum degree
// 			mindeg = func() int {
// 				if mindeg < d {
// 					return mindeg
// 				}
// 				return d
// 			}()
// 			degree[i] = d
// 			// place i in Lk
// 			Ci[func() int {
// 				defer func() {
// 					p++
// 				}()
// 				return p
// 			}()] = i
// 		}
//
// 		// # nodes absorbed into k
// 		nv[k] = nvk
// 		if (func() int {
// 			len[k] = p - pk1
// 			return len[k]
// 		}()) == 0 {
// 			// length of adj list of element k
// 			// k is a root of the tree
// 			Cp[k] = -1
// 			// k is now a dead element
// 			w[k] = 0
// 		}
// 		if elenk != 0 {
// 			// free unused space in Lk
// 			cnz = p
// 		}
// 	}
//
// 	// --- Postordering -----------------------------------------------------
// 	// fix assembly tree
// 	for i = 0; i < n; i++ {
// 		Cp[i] = -Cp[i] - 2
// 	}
//
// 	for j = 0; j <= n; j++ {
// 		head[j] = -1
// 	}
//
// 	// place unordered nodes in lists
// 	for j = n; j >= 0; j-- {
// 		if nv[j] > 0 {
// 			// skip if j is an element
// 			continue
// 		}
// 		// place j in list of its parent
// 		next[j] = head[Cp[j]]
// 		head[Cp[j]] = j
// 	}
//
// 	// place elements in lists
// 	for e = n; e >= 0; e-- {
// 		if nv[e] <= 0 {
// 			// skip unless e is an element
// 			continue
// 		}
// 		if Cp[e] != -1 {
// 			// place e in list of its parent
// 			next[e] = head[Cp[e]]
// 			head[Cp[e]] = e
// 		}
// 	}
//
// 	// postorder the assembly tree
// 	for i, k := 0, 0; i <= n; i++ {
// 		if Cp[i] == -1 {
// 			k = cs_tdfs(i, k, head, next, P, w)
// 		}
// 	}
//
// 	return cs_idone(P, C, W, true)
// }
//
// // cs_chol - L = chol (A, [pinv parent cp]), pinv is optional
// func cs_chol(A []cs, S []css) []csn {
// 	var d float64
// 	var lki float64
// 	var Lx []float64
// 	var x []float64
// 	var Cx []float64
// 	var top int
// 	var i int
// 	var p int
// 	var k int
// 	var n int
// 	var Li []int
// 	var Lp []int
// 	var cp []int
// 	var pinv []int
// 	var s []int
// 	var c []int
// 	var parent []int
// 	var Cp []int
// 	var Ci []int
// 	var L []cs
// 	var C []cs
// 	var E []cs
// 	var N []csn
// 	if !(A != nil && A.nz == -1) || S == nil || S.cp == nil || S.parent == nil {
// 		return nil
// 	}
// 	n = A.n
// 	// allocate result
// 	N = new(*csn) // cs_calloc(1, uint(32)).([]csn)
// 	// get csi workspace
// 	c = cs_malloc(noarch.PtrdiffT(2*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
// 	// get double workspace
// 	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
// 	cp = S[0].cp
// 	pinv = S[0].pinv
// 	parent = S[0].parent
// 	C = func() []cs {
// 		if pinv != nil {
// 			return cs_symperm(A, pinv, 1)
// 		}
// 		return (A)
// 	}()
// 	// E is alias for A, or a copy E=A(p,p)
// 	E = func() []cs {
// 		if pinv != nil {
// 			return C
// 		}
// 		return nil
// 	}()
// 	if N == nil || c == nil || x == nil || C == nil {
// 		return (cs_ndone(N, E, c, x, 0))
// 	}
// 	s = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&c[0])) + (uintptr)(int(n))*unsafe.Sizeof(c[0]))))[:]
// 	Cp = C[0].p
// 	Ci = C[0].i
// 	Cx = C[0].x
// 	L = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(cp[n]), 1, 0)
// 	// allocate result
// 	N[0].L = L
// 	if L == nil {
// 		return (cs_ndone(N, E, c, x, 0))
// 	}
// 	Lp = L[0].p
// 	Li = L[0].i
// 	Lx = L[0].x
// 	for k = 0; k < n; k++ {
// 		c[k] = cp[k]
// 		Lp[k] = c[k]
// 	}
// 	{
// 		// compute L(k,:) for L*L' = C
// 		for k = 0; k < n; k++ {
// 			// --- Nonzero pattern of L(k,:) ------------------------------------
// 			// find pattern of L(k,:)
// 			top = cs_ereach(C, noarch.PtrdiffT(k), parent, s, c)
// 			// x (0:k) is now zero
// 			x[k] = 0
// 			{
// 				// x = full(triu(C(:,k)))
// 				for p = Cp[k]; p < Cp[k+1]; p++ {
// 					if Ci[p] <= k {
// 						x[Ci[p]] = Cx[p]
// 					}
// 				}
// 			}
// 			// d = C(k,k)
// 			d = x[k]
// 			// clear x for k+1st iteration
// 			x[k] = 0
// 			for ; top < n; top++ {
// 				// --- Triangular solve ---------------------------------------------
// 				// solve L(0:k-1,0:k-1) * x = C(:,k)
// 				// s [top..n-1] is pattern of L(k,:)
// 				i = s[top]
// 				// L(k,i) = x (i) / L(i,i)
// 				lki = x[i] / Lx[Lp[i]]
// 				// clear x for k+1st iteration
// 				x[i] = 0
// 				for p = Lp[i] + 1; p < c[i]; p++ {
// 					x[Li[p]] -= Lx[p] * lki
// 				}
// 				// d = d - L(k,i)*L(k,i)
// 				d -= lki * lki
// 				p = func() noarch.PtrdiffT {
// 					tempVar := &c[i]
// 					defer func() {
// 						*tempVar++
// 					}()
// 					return *tempVar
// 				}()
// 				// store L(k,i) in column i
// 				Li[p] = k
// 				Lx[p] = lki
// 			}
// 			if d <= 0 {
// 				// --- Compute L(k,k) -----------------------------------------------
// 				// not pos def
// 				return (cs_ndone(N, E, c, x, 0))
// 			}
// 			p = func() noarch.PtrdiffT {
// 				tempVar := &c[k]
// 				defer func() {
// 					*tempVar++
// 				}()
// 				return *tempVar
// 			}()
// 			// store L(k,k) = sqrt (d) in column k
// 			Li[p] = k
// 			Lx[p] = math.Sqrt(d)
// 		}
// 	}
// 	// finalize L
// 	Lp[n] = cp[n]
// 	// success: free E,s,x; return N
// 	return (cs_ndone(N, E, c, x, 1))
// }
//
// // cs_cholsol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_cholsol.c:3
// // x=A\b where A is symmetric positive definite; b overwritten with solution
// func cs_cholsol(order noarch.PtrdiffT, A []cs, b []float64) noarch.PtrdiffT {
// 	var x []float64
// 	var S []css
// 	var N []csn
// 	var n noarch.PtrdiffT
// 	var ok noarch.PtrdiffT
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) || b == nil {
// 		// check inputs
// 		return 0
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	// ordering and symbolic analysis
// 	S = cs_schol(noarch.PtrdiffT(order), A)
// 	// numeric Cholesky factorization
// 	N = cs_chol(A, S)
// 	// get workspace
// 	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
// 	ok = noarch.PtrdiffT(S != nil && N != nil && x != nil)
// 	if bool(ok) {
// 		// x = P*b
// 		cs_ipvec(S[0].pinv, b, x, noarch.PtrdiffT(n))
// 		// x = L\x
// 		cs_lsolve(N[0].L, x)
// 		// x = L'\x
// 		cs_ltsolve(N[0].L, x)
// 		// b = P'*x
// 		cs_pvec(S[0].pinv, x, b, noarch.PtrdiffT(n))
// 	}
// 	cs_free(x)
// 	cs_sfree(S)
// 	cs_nfree(N)
// 	return noarch.PtrdiffT((ok))
// }

// cs_compress - C = compressed-column form of a triplet matrix T
func cs_compress(T *cs) *cs {
	var m int
	var n int
	var nz int
	var p int
	var k int
	var Cp []int
	var Ci []int
	var w []int
	var Ti []int
	var Tj []int
	var Cx []float64
	var Tx []float64
	var C *cs
	if !(T != nil && T.nz >= 0) {
		// check inputs
		return nil
	}
	m = T.m
	n = T.n
	Ti = T.i
	Tj = T.p
	Tx = T.x
	nz = T.nz
	// allocate result
	C = cs_spalloc(m, n, nz, Tx != nil, false)
	// get workspace
	w = make([]int, n)
	if C == nil || w == nil {
		// out of memory
		return cs_done(C, w, nil, false)
	}
	Cp = C.p
	Ci = C.i
	Cx = C.x
	{
		// column counts
		for k = 0; k < nz; k++ {
			w[Tj[k]]++
		}
	}
	// column pointers
	cs_cumsum(Cp, w, n)
	for k = 0; k < nz; k++ {
		// A(i,j) is the pth entry in C
		Ci[(func() int {
			p = func() int {
				tempVar := &w[Tj[k]]
				defer func() {
					*tempVar++
				}()
				return *tempVar
			}()
			return p
		}())] = Ti[k]
		if Cx != nil {
			Cx[p] = Tx[k]
		}
	}
	// success; free w and return C
	return cs_done(C, w, nil, true)
}

// // init_ata - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_counts.c:5
// // column counts of LL'=A or LL'=A'A, given parent & post ordering
// func init_ata(AT []cs, post []noarch.PtrdiffT, w []noarch.PtrdiffT, head [][]noarch.PtrdiffT, next [][]noarch.PtrdiffT) {
// 	var i noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var m noarch.PtrdiffT = noarch.PtrdiffT(AT[0].n)
// 	var n noarch.PtrdiffT = noarch.PtrdiffT(AT[0].m)
// 	var ATp []noarch.PtrdiffT = AT[0].p
// 	var ATi []noarch.PtrdiffT = AT[0].i
// 	head[0] = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(4*int32(n)))*unsafe.Sizeof(w[0]))))[:]
// 	next[0] = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(5*int32(n)))*unsafe.Sizeof(w[0]))))[:][1:]
// 	{
// 		// invert post
// 		for k = 0; k < n; k++ {
// 			w[post[k]] = k
// 		}
// 	}
// 	for i = 0; i < m; i++ {
// 		{
// 			k = n
// 			p = ATp[i]
// 			for p = ATp[i]; p < ATp[i+1]; p++ {
// 				k = noarch.PtrdiffT(func() int32 {
// 					if k < w[ATi[p]] {
// 						return int32(noarch.PtrdiffT((k)))
// 					}
// 					return int32(noarch.PtrdiffT((w[ATi[p]])))
// 				}() / 8)
// 			}
// 		}
// 		// place row i in linked list k
// 		(next[0])[i] = (head[0])[k]
// 		(head[0])[k] = i
// 	}
// }
//
// // cs_counts - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_counts.c:17
// func cs_counts(A []cs, parent []noarch.PtrdiffT, post []noarch.PtrdiffT, ata noarch.PtrdiffT) []noarch.PtrdiffT {
// 	var i noarch.PtrdiffT
// 	var j noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var m noarch.PtrdiffT
// 	var J noarch.PtrdiffT
// 	var s noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var q noarch.PtrdiffT
// 	var jleaf noarch.PtrdiffT
// 	var ATp []noarch.PtrdiffT
// 	var ATi []noarch.PtrdiffT
// 	var maxfirst []noarch.PtrdiffT
// 	var prevleaf []noarch.PtrdiffT
// 	var ancestor []noarch.PtrdiffT
// 	var head []noarch.PtrdiffT
// 	var next []noarch.PtrdiffT
// 	var colcount []noarch.PtrdiffT
// 	var w []noarch.PtrdiffT
// 	var first []noarch.PtrdiffT
// 	var delta []noarch.PtrdiffT
// 	var AT []cs
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) || parent == nil || post == nil {
// 		// check inputs
// 		return nil
// 	}
// 	m = noarch.PtrdiffT(A[0].m)
// 	n = noarch.PtrdiffT(A[0].n)
// 	s = noarch.PtrdiffT((4*int32(n) + func() int32 {
// 		if bool(noarch.PtrdiffT(ata)) {
// 			return (int32(n + m + 1))
// 		}
// 		return 0
// 	}()) / 8)
// 	colcount = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
// 	// allocate result
// 	delta = colcount
// 	// get workspace
// 	w = cs_malloc(noarch.PtrdiffT(s), uint(0)).([]noarch.PtrdiffT)
// 	// AT = A'
// 	AT = cs_transpose(A, 0)
// 	if AT == nil || colcount == nil || w == nil {
// 		return (cs_idone(colcount, AT, w, 0))
// 	}
// 	ancestor = w
// 	maxfirst = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
// 	prevleaf = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof(w[0]))))[:]
// 	first = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(3*int32(n)))*unsafe.Sizeof(w[0]))))[:]
// 	{
// 		// clear workspace w [0..s-1]
// 		for k = 0; k < s; k++ {
// 			w[k] = -1
// 		}
// 	}
// 	{
// 		// find first [j]
// 		for k = 0; k < n; k++ {
// 			j = post[k]
// 			// delta[j]=1 if j is a leaf
// 			delta[j] = noarch.PtrdiffT(func() int {
// 				if first[j] == -1 {
// 					return 1
// 				}
// 				return 0
// 			}())
// 			for ; j != -1 && first[j] == -1; j = parent[j] {
// 				first[j] = k
// 			}
// 		}
// 	}
// 	ATp = AT[0].p
// 	ATi = AT[0].i
// 	if bool(noarch.PtrdiffT(ata)) {
// 		init_ata(AT, post, w, (*[100000000][]noarch.PtrdiffT)(unsafe.Pointer(&head))[:], (*[100000000][]noarch.PtrdiffT)(unsafe.Pointer(&next))[:])
// 	}
// 	{
// 		// each node in its own set
// 		for i = 0; i < n; i++ {
// 			ancestor[i] = i
// 		}
// 	}
// 	for k = 0; k < n; k++ {
// 		// j is the kth node in postordered etree
// 		j = post[k]
// 		if parent[j] != -1 {
// 			// j is not a root
// 			delta[parent[j]]--
// 		}
// 		{
// 			// J=j for LL'=A case
// 			for J = noarch.PtrdiffT(func() int32 {
// 				if bool(noarch.PtrdiffT(ata)) {
// 					return int32(noarch.PtrdiffT(head[k]))
// 				}
// 				return int32(noarch.PtrdiffT(j))
// 			}() / 8); J != -1; J = noarch.PtrdiffT(func() int32 {
// 				if bool(noarch.PtrdiffT(ata)) {
// 					return int32(noarch.PtrdiffT(next[J]))
// 				}
// 				return int32(-1)
// 			}() / 8) {
// 				for p = ATp[J]; p < ATp[J+1]; p++ {
// 					i = ATi[p]
// 					q = cs_leaf(noarch.PtrdiffT(i), noarch.PtrdiffT(j), first, maxfirst, prevleaf, ancestor, (*[100000000]noarch.PtrdiffT)(unsafe.Pointer(&jleaf))[:])
// 					if jleaf >= 1 {
// 						// A(i,j) is in skeleton
// 						delta[j]++
// 					}
// 					if jleaf == 2 {
// 						// account for overlap in q
// 						delta[q]--
// 					}
// 				}
// 			}
// 		}
// 		if parent[j] != -1 {
// 			ancestor[j] = parent[j]
// 		}
// 	}
// 	{
// 		// sum up delta's of each child
// 		for j = 0; j < n; j++ {
// 			if parent[j] != -1 {
// 				colcount[parent[j]] += colcount[j]
// 			}
// 		}
// 	}
// 	// success: free workspace
// 	return (cs_idone(colcount, AT, w, 1))
// }

// cs_cumsum - p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
func cs_cumsum(p []int, c []int, n int) int64 {
	var nz int
	var nz2 int64
	if p == nil || c == nil {
		// check inputs
		return -1
	}
	for i := 0; i < n; i++ {
		p[i] = nz
		nz += c[i]
		// also in double to avoid csi overflow
		nz2 += int64(c[i])
		// also copy p[0..n-1] back into c[0..n-1]
		c[i] = p[i]
	}
	p[n] = nz
	// return sum (c [0..n-1])
	return nz2
}

// cs_dfs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dfs.c:3
// depth-first-search of the graph of a matrix, starting at node j
func cs_dfs(j int, G *cs, top int, xi []int, pstack []int, pinv []int) int {
	var i int
	var p int
	var p2 int
	var done bool
	var jnew int
	var head int
	if !(G != nil && G.nz == -1) || xi == nil || pstack == nil {
		// check inputs
		return -1
	}
	Gp := G.p
	Gi := G.i
	// initialize the recursion stack
	xi[0] = j
	for head >= 0 {
		// get j from the top of the recursion stack
		j = xi[head]
		jnew = func() int {
			if pinv != nil {
				return pinv[j]
			}
			return j
		}()
		if !(Gp[j] < 0) {
			{
				// mark node j as visited
				Gp[j] = -Gp[j] - 2
			}
			pstack[head] = func() int {
				if jnew < 0 {
					return 0
				}
				return (func() int {
					if Gp[jnew] < 0 {
						return -Gp[jnew] - 2
					}
					return Gp[jnew]
				}())
			}()
		}
		// node j done if no unvisited neighbors
		done = true
		p2 = func() int {
			if jnew < 0 {
				return 0
			}
			return (func() int {
				if Gp[jnew+1] < 0 {
					return -Gp[jnew+1] - 2
				}
				return Gp[jnew+1]
			}())
		}()

		// examine all neighbors of j
		for p = pstack[head]; p < p2; p++ {
			// consider neighbor node i
			i = Gi[p]
			if Gp[i] < 0 {
				// skip visited node i
				continue
			}
			// pause depth-first search of node j
			pstack[head] = p
			// start dfs at node i
			xi[func() int {
				head++
				return head
			}()] = i
			// node j is not done
			done = false
			// break, to start dfs (i)
			break
		}

		if done {
			// depth-first search at node j is done
			// remove j from the recursion stack
			head--
			// and place in the output stack
			xi[func() int {
				top--
				return top
			}()] = j
		}
	}
	return top
}

// cs_bfs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:3
// breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3)
func cs_bfs(A *cs,
	n int,
	wi []int,
	wj []int,
	queue []int,
	imatch []int,
	jmatch []int,
	mark int) bool {

	var head int
	var tail int
	var j2 int
	var C *cs

	// place all unmatched nodes in queue
	for j := 0; j < n; j++ {
		if imatch[j] >= 0 {
			// skip j if matched
			continue
		}
		// j in set C0 (R0 if transpose)
		wj[j] = 0
		// place unmatched col j in queue
		queue[func() int {
			defer func() {
				tail++
			}()
			return tail
		}()] = j
	}

	if tail == 0 {
		// quick return if no unmatched nodes
		return true
	}
	C = func() *cs {
		if mark == 1 {
			return A
		}
		return cs_transpose(A, false)
	}()
	if C == nil {
		// bfs of C=A' to find R3,C3 from R0
		return false
	}
	Ap := C.p
	Ai := C.i
	for head < tail {
		// while queue is not empty
		// get the head of the queue
		j := queue[func() int {
			defer func() {
				head++
			}()
			return head
		}()]
		for p := Ap[j]; p < Ap[j+1]; p++ {
			i := Ai[p]
			if wi[i] >= 0 {
				// skip if i is marked
				continue
			}
			// i in set R1 (C3 if transpose)
			wi[i] = mark
			// traverse alternating path to j2
			j2 = jmatch[i]
			if wj[j2] >= 0 {
				// skip j2 if it is marked
				continue
			}
			// j2 in set C1 (R3 if transpose)
			wj[j2] = mark
			// add j2 to queue
			queue[func() int {
				defer func() {
					tail++
				}()
				return tail
			}()] = j2
		}
	}
	if mark != 1 {
		// free A' if it was created
		cs_spfree(C) // TODO (KI) : remove
	}
	return true
}

// cs_matched - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:37
// collect matched rows and columns into p and q
func cs_matched(n int,
	wj []int,
	imatch []int,
	p []int,
	q []int,
	cc [5]int,
	rr [5]int,
	set int,
	mark int) {

	kc := cc[set]
	kr := rr[set-1]

	for j := 0; j < n; j++ {
		if wj[j] != mark {
			// skip if j is not in C set
			continue
		}
		p[func() int {
			defer func() {
				kr++
			}()
			return kr
		}()] = imatch[j]
		q[func() int {
			defer func() {
				kc++
			}()
			return kc
		}()] = j
	}
	cc[set+1] = kc
	rr[set] = kr
}

// cs_unmatched - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:53
// collect unmatched rows into the permutation vector p
func cs_unmatched(m int, wi []int, p []int, rr [5]int, set int) {
	kr := rr[set]
	for i := 0; i < m; i++ {
		if wi[i] == 0 {
			p[func() int {
				defer func() {
					kr++
				}()
				return kr
			}()] = i
		}
	}
	rr[set+1] = kr
}

// cs_rprune - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:61
// return 1 if row i is in R2
func cs_rprune(i, j int, aij float64, other interface{}) bool {
	rr := other.([5]int)
	return (i >= rr[1] && i < rr[2])
}

// cs_dmperm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:68
// Given A, compute coarse and then fine dmperm
func cs_dmperm(A *cs, seed int) *csd {
	var cnz int
	// var nc noarch.PtrdiffT
	// var pinv []noarch.PtrdiffT
	// var Cp []noarch.PtrdiffT
	// var Ci []noarch.PtrdiffT
	// var ps []noarch.PtrdiffT
	// var rs []noarch.PtrdiffT
	// var nb1 noarch.PtrdiffT
	// var nb2 noarch.PtrdiffT
	// var C []cs
	// var scc []csd
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	// --- Maximum matching -------------------------------------------------
	m := A.m
	n := A.n
	// allocate result
	D := cs_dalloc(m, n)
	if D == nil {
		return nil
	}
	p := D.p
	q := D.q
	r := D.r
	s := D.s
	cc := D.cc
	rr := D.rr
	// max transversal
	jmatch := cs_maxtrans(A, seed)
	// imatch = inverse of jmatch
	imatch := jmatch[m:] // jmatch+m
	if jmatch == nil {
		return cs_ddone(D, nil, jmatch, false)
	}
	// --- Coarse decomposition ---------------------------------------------
	// use r and s as workspace
	wi := r
	wj := s

	// unmark all cols for bfs
	for j := 0; j < n; j++ {
		wj[j] = -1
	}

	// unmark all rows for bfs
	for i := 0; i < m; i++ {
		wi[i] = -1
	}

	// find C1, R1 from C0
	cs_bfs(A, n, wi, wj, q, imatch, jmatch, 1)
	// find R3, C3 from R0
	ok := cs_bfs(A, m, wj, wi, p, jmatch, imatch, 3)
	if ok {
		return (cs_ddone(D, nil, jmatch, false))
	}
	// unmatched set C0
	cs_unmatched(n, wj, q, cc, 0)
	// set R1 and C1
	cs_matched(n, wj, imatch, p, q, cc, rr, 1, 1)
	// set R2 and C2
	cs_matched(n, wj, imatch, p, q, cc, rr, 2, -1)
	// set R3 and C3
	cs_matched(n, wj, imatch, p, q, cc, rr, 3, 3)
	// unmatched set R0
	cs_unmatched(m, wi, p, rr, 3)
	cs_free(jmatch)
	// --- Fine decomposition -----------------------------------------------
	// pinv=p'
	pinv := cs_pinv(p, m)
	if pinv == nil {
		return (cs_ddone(D, nil, nil, false))
	}
	// C=A(p,q) (it will hold A(R2,C2))
	C := cs_permute(A, pinv, q, false)
	cs_free(pinv)
	if C == nil {
		return (cs_ddone(D, nil, nil, false))
	}
	Cp := C.p
	// delete cols C0, C1, and C3 from C
	nc := cc[3] - cc[2]
	if cc[2] > 0 {
		for j := cc[2]; j <= cc[3]; j++ {
			Cp[j-cc[2]] = Cp[j]
		}
	}
	C.n = nc
	if rr[2]-rr[1] < m {
		// delete rows R0, R1, and R3 from C
		cs_fkeep(C, cs_rprune, rr)
		cnz = Cp[nc]
		Ci := C.i
		if rr[1] > 0 {
			for k := 0; k < cnz; k++ {
				Ci[k] -= rr[1]
			}
		}
	}
	C.m = nc
	// find strongly connected components of C
	scc := cs_scc(C)
	if scc == nil {
		return (cs_ddone(D, C, nil, false))
	}
	// --- Combine coarse and fine decompositions ---------------------------
	// C(ps,ps) is the permuted matrix
	ps := scc.p
	// kth block is rs[k]..rs[k+1]-1
	rs := scc.r
	// # of blocks of A(R2,C2)
	nb1 := scc.nb
	for k := 0; k < nc; k++ {
		wj[k] = q[ps[k]+cc[2]]
	}
	for k := 0; k < nc; k++ {
		q[k+cc[2]] = wj[k]
	}
	for k := 0; k < nc; k++ {
		wi[k] = p[ps[k]+rr[1]]
	}
	for k := 0; k < nc; k++ {
		p[k+rr[1]] = wi[k]
	}
	// create the fine block partitions
	nb2 := 0
	s[0] = 0
	r[0] = s[0]
	if cc[2] > 0 {
		// leading coarse block A (R1, [C0 C1])
		nb2++
	}

	// coarse block A (R2,C2)
	for k := 0; k < nb1; k++ {
		// A (R2,C2) splits into nb1 fine blocks
		r[nb2] = rs[k] + rr[1]
		s[nb2] = rs[k] + cc[2]
		nb2++
	}

	if rr[2] < m {
		// trailing coarse block A ([R3 R0], C3)
		r[nb2] = rr[2]
		s[nb2] = cc[3]
		nb2++
	}
	r[nb2] = m
	s[nb2] = n
	D.nb = nb2
	cs_dfree(scc)
	return (cs_ddone(D, C, nil, true))
}

// cs_tol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_droptol.c:2
func cs_tol(i, j int, aij float64, other interface{}) bool {
	return (math.Abs(aij) > *(other.(*float64)))
}

// cs_droptol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_droptol.c:6
func cs_droptol(A *cs, tol float64) int {
	// keep all large entries
	return cs_fkeep(A, cs_tol, &tol)
}

// cs_nonzero
func cs_nonzero(i, j int, aij float64, other interface{}) bool {
	return aij != 0
}

// cs_dropzeros
func cs_dropzeros(A *cs) int {
	// keep all nonzero entries
	return cs_fkeep(A, cs_nonzero, nil)
}

// cs_dupl - remove duplicate entries from A
func cs_dupl(A *cs) bool {
	var nz int
	if !(A != nil && A.nz == -1) {
		// check inputs
		return false
	}
	m := A.m
	n := A.n
	Ap := A.p
	Ai := A.i
	Ax := A.x
	// get workspace
	w := make([]int, m)
	if w == nil {
		// out of memory
		return false
	}

	// row i not yet seen
	for i := 0; i < m; i++ {
		w[i] = -1
	}

	for j := 0; j < n; j++ {
		// column j will start at q
		q := nz
		for p := Ap[j]; p < Ap[j+1]; p++ {
			// A(i,j) is nonzero
			i := Ai[p]
			if w[i] >= q {
				// A(i,j) is a duplicate
				Ax[w[i]] += Ax[p]
			} else {
				// record where row i occurs
				w[i] = nz
				// keep A(i,j)
				Ai[nz] = i
				Ax[func() int {
					defer func() {
						nz++
					}()
					return nz
				}()] = Ax[p]
			}
		}
		// record start of column j
		Ap[j] = q
	}
	// finalize A
	Ap[n] = nz
	// free workspace
	cs_free(w)
	// remove extra space from A
	return (cs_sprealloc(A, 0))
}

// cs_entry - add an entry to a triplet matrix; return 1 if ok, 0 otherwise
func cs_entry(T *cs, i, j int, x float64) bool {
	if !(T != nil && T.nz >= 0) || i < 0 || j < 0 {
		// check inputs
		return false
	}
	if T.nz >= T.nzmax && !cs_sprealloc(T, 2*T.nzmax) {
		return false
	}
	if T.x != nil {
		T.x[T.nz] = x
	}
	T.i[T.nz] = i
	T.p[T.nz] = j
	T.nz++
	if T.m < i+1 {
		T.m = i + 1
	}
	if T.n < j+1 {
		T.n = j + 1
	}
	return true
}

// // cs_ereach - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_ereach.c:3
// // find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k))
// func cs_ereach(A []cs, k noarch.PtrdiffT, parent []noarch.PtrdiffT, s []noarch.PtrdiffT, w []noarch.PtrdiffT) noarch.PtrdiffT {
// 	var i noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var len noarch.PtrdiffT
// 	var top noarch.PtrdiffT
// 	var Ap []noarch.PtrdiffT
// 	var Ai []noarch.PtrdiffT
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) || parent == nil || s == nil || w == nil {
// 		// check inputs
// 		return -1
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	top = n
// 	Ap = A[0].p
// 	Ai = A[0].i
// 	{
// 		// mark node k as visited
// 		w[k] = -noarch.PtrdiffT((w[k])) - 2
// 	}
// 	for p = Ap[k]; p < Ap[k+1]; p++ {
// 		// A(i,k) is nonzero
// 		i = Ai[p]
// 		if i > k {
// 			// only use upper triangular part of A
// 			continue
// 		}
// 		{
// 			// traverse up etree
// 			for len = 0; !(w[i] < 0); i = parent[i] {
// 				// L(k,i) is nonzero
// 				s[func() noarch.PtrdiffT {
// 					defer func() {
// 						len++
// 					}()
// 					return len
// 				}()] = i
// 				{
// 					// mark i as visited
// 					w[i] = -noarch.PtrdiffT((w[i])) - 2
// 				}
// 			}
// 		}
// 		for len > 0 {
// 			// push path onto stack
// 			s[func() noarch.PtrdiffT {
// 				top--
// 				return top
// 			}()] = s[func() noarch.PtrdiffT {
// 				len--
// 				return len
// 			}()]
// 		}
// 	}
// 	{
// 		// unmark all nodes
// 		for p = top; p < n; p++ {
// 			w[s[p]] = -noarch.PtrdiffT((w[s[p]])) - 2
// 		}
// 	}
// 	{
// 		// unmark node k
// 		w[k] = -noarch.PtrdiffT((w[k])) - 2
// 	}
// 	// s [top..n-1] contains pattern of L(k,:)
// 	return noarch.PtrdiffT((top))
// }
//
// // cs_etree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_etree.c:3
// // compute the etree of A (using triu(A), or A'A without forming A'A
// func cs_etree(A []cs, ata noarch.PtrdiffT) []noarch.PtrdiffT {
// 	var i noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var m noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var inext noarch.PtrdiffT
// 	var Ap []noarch.PtrdiffT
// 	var Ai []noarch.PtrdiffT
// 	var w []noarch.PtrdiffT
// 	var parent []noarch.PtrdiffT
// 	var ancestor []noarch.PtrdiffT
// 	var prev []noarch.PtrdiffT
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) {
// 		// check inputs
// 		return nil
// 	}
// 	m = noarch.PtrdiffT(A[0].m)
// 	n = noarch.PtrdiffT(A[0].n)
// 	Ap = A[0].p
// 	Ai = A[0].i
// 	// allocate result
// 	parent = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
// 	// get workspace
// 	w = cs_malloc(n+noarch.PtrdiffT(func() int32 {
// 		if bool(noarch.PtrdiffT(ata)) {
// 			return int32(m)
// 		}
// 		return 0
// 	}()/8), uint(0)).([]noarch.PtrdiffT)
// 	if w == nil || parent == nil {
// 		return (cs_idone(parent, nil, w, 0))
// 	}
// 	ancestor = w
// 	prev = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
// 	if bool(noarch.PtrdiffT(ata)) {
// 		for i = 0; i < m; i++ {
// 			prev[i] = -1
// 		}
// 	}
// 	for k = 0; k < n; k++ {
// 		// node k has no parent yet
// 		parent[k] = -1
// 		// nor does k have an ancestor
// 		ancestor[k] = -1
// 		for p = Ap[k]; p < Ap[k+1]; p++ {
// 			i = noarch.PtrdiffT(func() int32 {
// 				if bool(noarch.PtrdiffT(ata)) {
// 					return int32(noarch.PtrdiffT((prev[Ai[p]])))
// 				}
// 				return int32(noarch.PtrdiffT((Ai[p])))
// 			}() / 8)
// 			{
// 				// traverse from i to k
// 				for ; i != -1 && i < k; i = inext {
// 					// inext = ancestor of i
// 					inext = ancestor[i]
// 					// path compression
// 					ancestor[i] = k
// 					if inext == -1 {
// 						// no anc., parent is k
// 						parent[i] = k
// 					}
// 				}
// 			}
// 			if bool(noarch.PtrdiffT(ata)) {
// 				prev[Ai[p]] = k
// 			}
// 		}
// 	}
// 	return (cs_idone(parent, nil, w, 1))
// }

// cs_fkeep - drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1
func cs_fkeep(A *cs, fkeep func(int, int, float64, interface{}) bool, other interface{}) int {
	var nz int
	if !(A != nil && A.nz == -1) || fkeep == nil {
		// check inputs
		return -1
	}
	n := A.n
	Ap := A.p
	Ai := A.i
	Ax := A.x
	for j := 0; j < n; j++ {
		// get current location of col j
		p := Ap[j]
		// record new location of col j
		Ap[j] = nz
		for ; p < Ap[j+1]; p++ {
			if fkeep(Ai[p], j, func() float64 {
				if Ax != nil {
					return Ax[p]
				}
				return 1
			}(), other) {
				if Ax != nil {
					// keep A(i,j)
					Ax[nz] = Ax[p]
				}
				Ai[func() int {
					defer func() {
						nz++
					}()
					return nz
				}()] = Ai[p]
			}
		}
	}
	// finalize A
	Ap[n] = nz
	// remove extra space from A
	cs_sprealloc(A, 0)
	return nz
}

// cs_gaxpy - y = A*x+y
func cs_gaxpy(A *cs, x []float64, y []float64) bool {
	if !(A != nil && A.nz == -1) || x == nil || y == nil {
		// check inputs
		return false
	}
	n := A.n
	Ap := A.p
	Ai := A.i
	Ax := A.x
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			y[Ai[p]] += Ax[p] * x[j]
		}
	}
	return true
}

// // cs_happly - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_happly.c:3
// // apply the ith Householder vector to x
// func cs_happly(V []cs, i noarch.PtrdiffT, beta float64, x []float64) noarch.PtrdiffT {
// 	var p noarch.PtrdiffT
// 	var Vp []noarch.PtrdiffT
// 	var Vi []noarch.PtrdiffT
// 	var Vx []float64
// 	var tau float64
// 	if !(V != nil && noarch.PtrdiffT(V[0].nz) == -1) || x == nil {
// 		// check inputs
// 		return 0
// 	}
// 	Vp = V[0].p
// 	Vi = V[0].i
// 	Vx = V[0].x
// 	{
// 		// tau = v'*x
// 		for p = Vp[i]; p < Vp[i+1]; p++ {
// 			tau += Vx[p] * x[Vi[p]]
// 		}
// 	}
// 	// tau = beta*(v'*x)
// 	tau *= beta
// 	{
// 		// x = x - v*tau
// 		for p = Vp[i]; p < Vp[i+1]; p++ {
// 			x[Vi[p]] -= Vx[p] * tau
// 		}
// 	}
// 	return 1
// }
//
// // cs_house - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_house.c:4
// // create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
// // * where (I-beta*v*v')*x = s*e1.  See Algo 5.1.1, Golub & Van Loan, 3rd ed.
// func cs_house(x []float64, beta []float64, n noarch.PtrdiffT) float64 {
// 	var s float64
// 	var sigma float64
// 	var i noarch.PtrdiffT
// 	if x == nil || beta == nil {
// 		// check inputs
// 		return float64((-1))
// 	}
// 	for i = 1; i < n; i++ {
// 		sigma += x[i] * x[i]
// 	}
// 	if sigma == 0 {
// 		// s = |x(0)|
// 		s = math.Abs(x[0])
// 		beta[0] = float64(func() int {
// 			if x[0] <= 0 {
// 				return 2
// 			}
// 			return 0
// 		}())
// 		x[0] = 1
// 	} else {
// 		// s = norm (x)
// 		s = math.Sqrt(x[0]*x[0] + sigma)
// 		x[0] = func() float64 {
// 			if x[0] <= 0 {
// 				return (x[0] - s)
// 			}
// 			return (-sigma / (x[0] + s))
// 		}()
// 		beta[0] = -1 / (s * x[0])
// 	}
// 	return (s)
// }
//
// // cs_ipvec - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_ipvec.c:3
// // x(p) = b, for dense vectors x and b; p=NULL denotes identity
// func cs_ipvec(p []noarch.PtrdiffT, b []float64, x []float64, n noarch.PtrdiffT) noarch.PtrdiffT {
// 	var k noarch.PtrdiffT
// 	if x == nil || b == nil {
// 		// check inputs
// 		return 0
// 	}
// 	for k = 0; k < n; k++ {
// 		x[func() int32 {
// 			if p != nil {
// 				return int32(noarch.PtrdiffT(p[k]))
// 			}
// 			return int32(noarch.PtrdiffT(k))
// 		}()] = b[k]
// 	}
// 	return 1
// }
//
// // cs_leaf - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_leaf.c:3
// // consider A(i,j), node j in ith row subtree and return lca(jprev,j)
// func cs_leaf(i noarch.PtrdiffT, j noarch.PtrdiffT, first []noarch.PtrdiffT, maxfirst []noarch.PtrdiffT, prevleaf []noarch.PtrdiffT, ancestor []noarch.PtrdiffT, jleaf []noarch.PtrdiffT) noarch.PtrdiffT {
// 	var q noarch.PtrdiffT
// 	var s noarch.PtrdiffT
// 	var sparent noarch.PtrdiffT
// 	var jprev noarch.PtrdiffT
// 	if first == nil || maxfirst == nil || prevleaf == nil || ancestor == nil || jleaf == nil {
// 		return -1
// 	}
// 	jleaf[0] = 0
// 	if i <= j || first[j] <= maxfirst[i] {
// 		// j not a leaf
// 		return -1
// 	}
// 	// update max first[j] seen so far
// 	maxfirst[i] = first[j]
// 	// jprev = previous leaf of ith subtree
// 	jprev = prevleaf[i]
// 	prevleaf[i] = j
// 	// j is first or subsequent leaf
// 	jleaf[0] = noarch.PtrdiffT(func() int {
// 		if jprev == -1 {
// 			return 1
// 		}
// 		return 2
// 	}())
// 	if jleaf[0] == 1 {
// 		// if 1st leaf, q = root of ith subtree
// 		return noarch.PtrdiffT((i))
// 	}
// 	for q = jprev; q != ancestor[q]; q = ancestor[q] {
// 	}
// 	for s = jprev; s != q; s = sparent {
// 		// path compression
// 		sparent = ancestor[s]
// 		ancestor[s] = q
// 	}
// 	// q = least common ancester (jprev,j)
// 	return noarch.PtrdiffT((q))
// }
//

// cs_load - load a triplet matrix from a file
func cs_load(f io.Reader) *cs {
	var T *cs
	if f == nil {
		// use double for integers to avoid csi conflicts
		// check inputs
		return nil
	}
	// allocate result
	T = cs_spalloc(0, 0, 1, true, true)
	for {
		var i, j int
		var x float64

		n, err := fmt.Fscanf(f, "%d %d %f\n", &i, &j, &x)
		if err == io.EOF {
			break
		}
		if err != nil || n != 3 {
			return nil
		}
		if !cs_entry(T, i, j, x) {
			return nil
		}
	}
	return T
}

//
// // cs_lsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_lsolve.c:3
// // solve Lx=b where x and b are dense.  x=b on input, solution on output.
// func cs_lsolve(L []cs, x []float64) noarch.PtrdiffT {
// 	var p noarch.PtrdiffT
// 	var j noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var Lp []noarch.PtrdiffT
// 	var Li []noarch.PtrdiffT
// 	var Lx []float64
// 	if !(L != nil && noarch.PtrdiffT(L[0].nz) == -1) || x == nil {
// 		// check inputs
// 		return 0
// 	}
// 	n = noarch.PtrdiffT(L[0].n)
// 	Lp = L[0].p
// 	Li = L[0].i
// 	Lx = L[0].x
// 	for j = 0; j < n; j++ {
// 		x[j] /= Lx[Lp[j]]
// 		for p = Lp[j] + 1; p < Lp[j+1]; p++ {
// 			x[Li[p]] -= Lx[p] * x[j]
// 		}
// 	}
// 	return 1
// }
//
// // cs_ltsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_ltsolve.c:3
// // solve L'x=b where x and b are dense.  x=b on input, solution on output.
// func cs_ltsolve(L []cs, x []float64) noarch.PtrdiffT {
// 	var p noarch.PtrdiffT
// 	var j noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var Lp []noarch.PtrdiffT
// 	var Li []noarch.PtrdiffT
// 	var Lx []float64
// 	if !(L != nil && noarch.PtrdiffT(L[0].nz) == -1) || x == nil {
// 		// check inputs
// 		return 0
// 	}
// 	n = noarch.PtrdiffT(L[0].n)
// 	Lp = L[0].p
// 	Li = L[0].i
// 	Lx = L[0].x
// 	for j = n - 1; j >= 0; j-- {
// 		for p = Lp[j] + 1; p < Lp[j+1]; p++ {
// 			x[j] -= Lx[p] * x[Li[p]]
// 		}
// 		x[j] /= Lx[Lp[j]]
// 	}
// 	return 1
// }
//
// // cs_lu - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_lu.c:3
// // [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guess
// func cs_lu(A []cs, S []css, tol float64) []csn {
// 	var L []cs
// 	var U []cs
// 	var N []csn
// 	var pivot float64
// 	var Lx []float64
// 	var Ux []float64
// 	var x []float64
// 	var a float64
// 	var t float64
// 	var Lp []noarch.PtrdiffT
// 	var Li []noarch.PtrdiffT
// 	var Up []noarch.PtrdiffT
// 	var Ui []noarch.PtrdiffT
// 	var pinv []noarch.PtrdiffT
// 	var xi []noarch.PtrdiffT
// 	var q []noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var ipiv noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var top noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var i noarch.PtrdiffT
// 	var col noarch.PtrdiffT
// 	var lnz noarch.PtrdiffT
// 	var unz noarch.PtrdiffT
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) || S == nil {
// 		// check inputs
// 		return nil
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	q = S[0].q
// 	lnz = noarch.PtrdiffT(S[0].lnz)
// 	unz = noarch.PtrdiffT(S[0].unz)
// 	// get double workspace
// 	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
// 	// get csi workspace
// 	xi = cs_malloc(noarch.PtrdiffT(2*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
// 	// allocate result
// 	N = cs_calloc(1, uint(32)).([]csn)
// 	if x == nil || xi == nil || N == nil {
// 		return (cs_ndone(N, nil, xi, x, 0))
// 	}
// 	L = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(lnz), 1, 0)
// 	// allocate result L
// 	N[0].L = L
// 	U = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(unz), 1, 0)
// 	// allocate result U
// 	N[0].U = U
// 	pinv = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
// 	// allocate result pinv
// 	N[0].pinv = pinv
// 	if L == nil || U == nil || pinv == nil {
// 		return (cs_ndone(N, nil, xi, x, 0))
// 	}
// 	Lp = L[0].p
// 	Up = U[0].p
// 	{
// 		// clear workspace
// 		for i = 0; i < n; i++ {
// 			x[i] = 0
// 		}
// 	}
// 	{
// 		// no rows pivotal yet
// 		for i = 0; i < n; i++ {
// 			pinv[i] = -1
// 		}
// 	}
// 	{
// 		// no cols of L yet
// 		for k = 0; k <= n; k++ {
// 			Lp[k] = 0
// 		}
// 	}
// 	unz = 0
// 	lnz = unz
// 	{
// 		// compute L(:,k) and U(:,k)
// 		for k = 0; k < n; k++ {
// 			// --- Triangular solve ---------------------------------------------
// 			// L(:,k) starts here
// 			Lp[k] = lnz
// 			// U(:,k) starts here
// 			Up[k] = unz
// 			if lnz+n > noarch.PtrdiffT(L[0].nzmax) && bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(L, noarch.PtrdiffT((2*int32(noarch.PtrdiffT(L[0].nzmax))+int32(n))/8)))) || unz+n > noarch.PtrdiffT(U[0].nzmax) && bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(U, noarch.PtrdiffT((2*int32(noarch.PtrdiffT(U[0].nzmax))+int32(n))/8)))) {
// 				return (cs_ndone(N, nil, xi, x, 0))
// 			}
// 			Li = L[0].i
// 			Lx = L[0].x
// 			Ui = U[0].i
// 			Ux = U[0].x
// 			col = noarch.PtrdiffT(func() int32 {
// 				if q != nil {
// 					return int32(noarch.PtrdiffT((q[k])))
// 				}
// 				return int32(noarch.PtrdiffT(k))
// 			}() / 8)
// 			// x = L\A(:,col)
// 			top = cs_spsolve(L, A, noarch.PtrdiffT(col), xi, x, pinv, 1)
// 			// --- Find pivot ---------------------------------------------------
// 			ipiv = -1
// 			a = float64(-1)
// 			for p = top; p < n; p++ {
// 				// x(i) is nonzero
// 				i = xi[p]
// 				if pinv[i] < 0 {
// 					if (func() float64 {
// 						t = math.Abs(x[i])
// 						return t
// 					}()) > a {
// 						// row i is not yet pivotal
// 						// largest pivot candidate so far
// 						a = t
// 						ipiv = i
// 					}
// 				} else {
// 					// x(i) is the entry U(pinv[i],k)
// 					Ui[unz] = pinv[i]
// 					Ux[func() noarch.PtrdiffT {
// 						defer func() {
// 							unz++
// 						}()
// 						return unz
// 					}()] = x[i]
// 				}
// 			}
// 			if ipiv == -1 || a <= 0 {
// 				return (cs_ndone(N, nil, xi, x, 0))
// 			}
// 			if pinv[col] < 0 && math.Abs(x[col]) >= a*tol {
// 				// tol=1 for  partial pivoting; tol<1 gives preference to diagonal
// 				ipiv = col
// 			}
// 			// --- Divide by pivot ----------------------------------------------
// 			// the chosen pivot
// 			pivot = x[ipiv]
// 			// last entry in U(:,k) is U(k,k)
// 			Ui[unz] = k
// 			Ux[func() noarch.PtrdiffT {
// 				defer func() {
// 					unz++
// 				}()
// 				return unz
// 			}()] = pivot
// 			// ipiv is the kth pivot row
// 			pinv[ipiv] = k
// 			// first entry in L(:,k) is L(k,k) = 1
// 			Li[lnz] = ipiv
// 			Lx[func() noarch.PtrdiffT {
// 				defer func() {
// 					lnz++
// 				}()
// 				return lnz
// 			}()] = 1
// 			{
// 				// L(k+1:n,k) = x / pivot
// 				for p = top; p < n; p++ {
// 					i = xi[p]
// 					if pinv[i] < 0 {
// 						// x(i) is an entry in L(:,k)
// 						// save unpermuted row in L
// 						Li[lnz] = i
// 						// scale pivot column
// 						Lx[func() noarch.PtrdiffT {
// 							defer func() {
// 								lnz++
// 							}()
// 							return lnz
// 						}()] = x[i] / pivot
// 					}
// 					// x [0..n-1] = 0 for next k
// 					x[i] = 0
// 				}
// 			}
// 		}
// 	}
// 	// --- Finalize L and U -------------------------------------------------
// 	Lp[n] = lnz
// 	Up[n] = unz
// 	// fix row indices of L for final pinv
// 	Li = L[0].i
// 	for p = 0; p < lnz; p++ {
// 		Li[p] = pinv[Li[p]]
// 	}
// 	// remove extra space from L and U
// 	cs_sprealloc(L, 0)
// 	cs_sprealloc(U, 0)
// 	// success
// 	return (cs_ndone(N, nil, xi, x, 1))
// }

// cs_lusol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_lusol.c:3
// x=A\b where A is unsymmetric; b overwritten with solution
func cs_lusol(order int, A *cs, b []float64, tol float64) bool {
	if !(A != nil && A.nz == -1) || b == nil {
		// check inputs
		return false
	}
	n := A.n
	// ordering and symbolic analysis
	S := cs_sqr(order, A, 0)
	// numeric LU factorization
	N := cs_lu(A, S, tol)
	// get workspace
	x := make([]float64, n)
	ok := (S != nil && N != nil && x != nil)
	if ok {
		// x = b(p)
		cs_ipvec(N.pinv, b, x, noarch.PtrdiffT(n))
		// x = L\x
		cs_lsolve(N.L, x)
		// x = U\x
		cs_usolve(N.U, x)
		// b(q) = x
		cs_ipvec(S.q, x, b, noarch.PtrdiffT(n))
	}
	cs_free(x)
	cs_sfree(S)
	cs_nfree(N)
	return ok
}

// // cs_malloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_malloc.c:10
// // wrapper for malloc
// func cs_malloc(n noarch.PtrdiffT, size uint) interface{} {
// 	return (make([]byte, uint32(func() int32 {
// 		if n > 1 {
// 			return int32(noarch.PtrdiffT((n)))
// 		}
// 		return 1
// 	}())*uint32(size)))
// }
//
// // cs_calloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_malloc.c:16
// // wrapper for calloc
// func cs_calloc(n noarch.PtrdiffT, size uint) interface{} {
// 	return (make([]byte, (size)*(uint(uint32((func() int32 {
// 		if n > 1 {
// 			return int32(noarch.PtrdiffT((n)))
// 		}
// 		return 1
// 	}()))))))
// }

// cs_free - wrapper for free
func cs_free(p interface{}) interface{} {
	if p != nil {
		_ = p
		// free p if it is not already NULL
	}
	// return NULL to simplify the use of cs_free
	return nil
}

// cs_realloc - wrapper for realloc
func cs_realloc(p interface{}, n int, ok *bool) interface{} {
	//
	// TODO (KI): redesign
	//

	switch v := p.(type) {
	case []int:
		if len(v) <= n {
			v = append(v, make([]int, n-len(v))...)
		}
		*ok = true
		return v

	case []float64:
		if len(v) <= n {
			v = append(v, make([]float64, n-len(v))...)
		}
		*ok = true
		return v
	}
	return nil

}

// cs_augment - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_maxtrans.c:3
// find an augmenting path starting at column k and extend the match if found
func cs_augment(k int,
	A *cs,
	jmatch []int,
	cheap []int,
	w []int,
	js []int,
	is []int,
	ps []int) {

	var found bool = false
	var p int
	var i int = -1
	Ap := A.p
	Ai := A.i
	var head int
	var j int
	// start with just node k in jstack
	js[0] = k
	for head >= 0 {
		// --- Start (or continue) depth-first-search at node j -------------
		// get j from top of jstack
		j = js[head]
		if w[j] != k {
			// 1st time j visited for kth path
			// mark j as visited for kth path
			w[j] = k
			for p = cheap[j]; p < Ap[j+1] && !found; p++ {
				// try a cheap assignment (i,j)
				i = Ai[p]
				found = (jmatch[i] == -1)
			}
			// start here next time j is traversed
			cheap[j] = p
			if found {
				// column j matched with row i
				is[head] = i
				// end of augmenting path
				break
			}
			// no cheap match: start dfs for j
			ps[head] = Ap[j]
		}

		// --- Depth-first-search of neighbors of j -------------------------
		for p = ps[head]; p < Ap[j+1]; p++ {
			// consider row i
			i = Ai[p]
			if w[jmatch[i]] == k {
				// skip jmatch [i] if marked
				continue
			}
			// pause dfs of node j
			ps[head] = p + 1
			// i will be matched with j if found
			is[head] = i
			// start dfs at column jmatch [i]
			js[func() int {
				head++
				return head
			}()] = jmatch[i]
			break
		}

		if p == Ap[j+1] {
			// node j is done; pop from stack
			head--
		}
	}
	if found {
		// augment the match if path found:
		for p = head; p >= 0; p-- {
			jmatch[is[p]] = js[p]
		}
	}
}

// cs_maxtrans - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_maxtrans.c:44
// find a maximum transveral
//[jmatch [0..m-1]; imatch [0..n-1]]
func cs_maxtrans(A *cs, seed int) []int {
	var i int
	var j int
	var k int
	// var n noarch.PtrdiffT
	// var m noarch.PtrdiffT
	var p int
	// var n2 int
	var m2 int
	// var Ap []noarch.PtrdiffT
	// var jimatch []int
	// var w []noarch.PtrdiffT
	var cheap []int
	var js []int
	var is []int
	var ps []int
	// var Ai []noarch.PtrdiffT
	var Cp []int
	var jmatch []int
	var imatch []int
	var q []int
	// var C []cs
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n := A.n
	m := A.m
	Ap := A.p
	Ai := A.i
	jimatch := make([]int, m+n) // cs_calloc(m+n, uint(0)).([]noarch.PtrdiffT)
	// allocate result
	w := jimatch
	if jimatch == nil {
		return nil
	}

	// count nonempty rows and columns
	n2 := 0
	for j, k = 0, 0; j < n; j++ {
		if Ap[j] < Ap[j+1] {
			n2++
		}
		for p = Ap[j]; p < Ap[j+1]; p++ {
			w[Ai[p]] = 1
			// count entries already on diagonal
			if j == Ai[p] {
				k++
			}
		}
	}

	if k == func() int {
		if m < n {
			return m
		}
		return n
	}() {
		// quick return if diagonal zero-free
		jmatch = jimatch
		imatch = jimatch[m:]
		for i = 0; i < k; i++ {
			jmatch[i] = i
		}
		for ; i < m; i++ {
			jmatch[i] = -1
		}
		for j = 0; j < k; j++ {
			imatch[j] = j
		}
		for ; j < n; j++ {
			imatch[j] = -1
		}
		return (cs_idone(jimatch, nil, nil, true))
	}
	for i = 0; i < m; i++ {
		m2 += w[i]
	}
	// transpose if needed
	C := func() *cs {
		if m2 < n2 {
			return cs_transpose(A, false)
		}
		return (A)
	}()
	if C == nil {
		return (cs_idone(jimatch, func() *cs {
			if m2 < n2 {
				return C
			}
			return nil
		}(), nil, false))
	}
	n = C.n
	m = C.m
	Cp = C.p
	jmatch = func() []int {
		if m2 < n2 {
			return jmatch[n:]
		}
		return jimatch
	}()
	imatch = func() []int {
		if m2 < n2 {
			return jimatch
		}
		return jimatch[m:]
	}()
	// get workspace
	w = make([]int, 5*n) // cs_malloc(noarch.PtrdiffT(5*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
	if w == nil {
		return (cs_idone(jimatch, func() *cs {
			if m2 < n2 {
				return C
			}
			return nil
		}(), w, false))
	}
	cheap = w[n:]
	js = w[2*n:]
	is = w[3*n:]
	ps = w[4*n:]

	// for cheap assignment
	for j = 0; j < n; j++ {
		cheap[j] = Cp[j]
	}

	// all columns unflagged
	for j = 0; j < n; j++ {
		w[j] = -1
	}

	// nothing matched yet
	for i = 0; i < m; i++ {
		jmatch[i] = -1
	}

	// q = random permutation
	q = cs_randperm(n, seed)

	// augment, starting at column q[k]
	for k = 0; k < n; k++ {
		cs_augment(func() int {
			if q != nil {
				return q[k]
			}
			return k
		}(), C, jmatch, cheap, w, js, is, ps)
	}

	cs_free(q)

	// find row match
	for j = 0; j < n; j++ {
		imatch[j] = -1
	}

	for i = 0; i < m; i++ {
		if jmatch[i] >= 0 {
			imatch[jmatch[i]] = i
		}
	}
	return (cs_idone(jimatch, func() *cs {
		if m2 < n2 {
			return C
		}
		return nil
	}(), w, true))
}

// cs_multiply - C = A*B
func cs_multiply(A *cs, B *cs) *cs {
	var p int
	var nz int
	var anz int
	var Cp []int
	var Ci []int
	var Bp []int
	var m int
	var n int
	var bnz int
	var w []int
	var Bi []int
	var x []float64
	var Bx []float64
	var Cx []float64
	var C *cs
	if !(A != nil && A.nz == -1) || !(B != nil && B.nz == -1) {
		// check inputs
		return nil
	}
	if A.n != B.m {
		return nil
	}
	m = A.m
	anz = A.p[A.n]
	n = B.n
	Bp = B.p
	Bi = B.i
	Bx = B.x
	bnz = Bp[n]
	// get workspace
	w = make([]int, m)
	values := (A.x != nil && Bx != nil)
	// get workspace
	if values {
		x = make([]float64, m)
	}
	// allocate result
	C = cs_spalloc(m, n, anz+bnz, values, false)
	if C == nil || w == nil || bool(values) && x == nil {
		return cs_done(C, w, x, false)
	}
	Cp = C.p
	for j := 0; j < n; j++ {
		if nz+m > C.nzmax && !cs_sprealloc(C, 2*C.nzmax+m) {
			// out of memory
			return cs_done(C, w, x, false)
		}
		// C->i and C->x may be reallocated
		Ci = C.i
		Cx = C.x
		// column j of C starts here
		Cp[j] = nz
		for p = Bp[j]; p < Bp[j+1]; p++ {
			nz = cs_scatter(A, Bi[p], func() float64 {
				if Bx != nil {
					return Bx[p]
				}
				return 1
			}(), w, x, j+1, C, nz)
		}
		if values {
			for p := Cp[j]; p < nz; p++ {
				Cx[p] = x[Ci[p]]
			}
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	// remove extra space from C
	cs_sprealloc(C, 0)
	// success; free workspace, return C
	return cs_done(C, w, x, true)
}

// cs_norm - 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum
func cs_norm(A *cs) float64 {
	var norm float64
	if !(A != nil && A.nz == -1) || A.x == nil {
		// check inputs
		return -1
	}
	n := A.n
	Ap := A.p
	Ax := A.x
	for j := 0; j < n; j++ {

		s := 0.0
		for p := Ap[j]; p < Ap[j+1]; p++ {
			s += math.Abs(Ax[p])
		}

		norm = func() float64 {
			if norm > s {
				return (norm)
			}
			return s
		}()
	}
	return norm
}

// cs_permute - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_permute.c:3
// C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1.
func cs_permute(A *cs, pinv []int, q []int, values bool) *cs {
	nz := 0
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	m := A.m
	n := A.n
	Ap := A.p
	Ai := A.i
	Ax := A.x
	// alloc result
	C := cs_spalloc(m, n, Ap[n], values && Ax != nil, false)
	if C == nil {
		// out of memory
		return cs_done(C, nil, nil, false)
	}
	Cp := C.p
	Ci := C.i
	Cx := C.x
	for k := 0; k < n; k++ {
		// column k of C is column q[k] of A
		Cp[k] = nz
		j := func() int {
			if q != nil {
				return q[k]
			}
			return k
		}()
		for t := Ap[j]; t < Ap[j+1]; t++ {
			if Cx != nil {
				// row i of A is row pinv[i] of C
				Cx[nz] = Ax[t]
			}
			Ci[func() int {
				defer func() {
					nz++
				}()
				return nz
			}()] = func() int {
				if pinv != nil {
					return pinv[Ai[t]]
				}
				return Ai[t]
			}()
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	return cs_done(C, nil, nil, true)
}

// cs_pinv - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_pinv.c:3
// pinv = p', or p = pinv'
func cs_pinv(p []int, n int) []int {
	var pinv []int
	if p == nil {
		// p = NULL denotes identity
		return nil
	}
	// allocate result
	pinv = make([]int, n) // cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	if pinv == nil {
		// out of memory
		return nil
	}

	// invert the permutation
	for k := 0; k < n; k++ {
		pinv[p[k]] = k
	}

	// return result
	return (pinv)
}

// // cs_post - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_post.c:3
// // post order a forest
// func cs_post(parent []noarch.PtrdiffT, n noarch.PtrdiffT) []noarch.PtrdiffT {
// 	var j noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var post []noarch.PtrdiffT
// 	var w []noarch.PtrdiffT
// 	var head []noarch.PtrdiffT
// 	var next []noarch.PtrdiffT
// 	var stack []noarch.PtrdiffT
// 	if parent == nil {
// 		// check inputs
// 		return nil
// 	}
// 	// allocate result
// 	post = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
// 	// get workspace
// 	w = cs_malloc(noarch.PtrdiffT(3*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
// 	if w == nil || post == nil {
// 		return (cs_idone(post, nil, w, 0))
// 	}
// 	head = w
// 	next = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
// 	stack = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof(w[0]))))[:]
// 	{
// 		// empty linked lists
// 		for j = 0; j < n; j++ {
// 			head[j] = -1
// 		}
// 	}
// 	{
// 		// traverse nodes in reverse order
// 		for j = n - 1; j >= 0; j-- {
// 			if parent[j] == -1 {
// 				// j is a root
// 				continue
// 			}
// 			// add j to list of its parent
// 			next[j] = head[parent[j]]
// 			head[parent[j]] = j
// 		}
// 	}
// 	for j = 0; j < n; j++ {
// 		if parent[j] != -1 {
// 			// skip j if it is not a root
// 			continue
// 		}
// 		k = cs_tdfs(noarch.PtrdiffT(j), noarch.PtrdiffT(k), head, next, post, stack)
// 	}
// 	// success; free w, return post
// 	return (cs_idone(post, nil, w, 1))
// }

// cs_print - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_print.c:3
// print a sparse matrix; use %g for integers to avoid differences with csi
func cs_print(A *cs, brief bool) bool {
	var p int
	var m int
	var n int
	var nzmax int
	var nz int
	var Ap []int
	var Ai []int
	var Ax []float64
	if A == nil {
		fmt.Printf("(null)\n")
		return false
	}
	m = A.m
	n = A.n
	Ap = A.p
	Ai = A.i
	Ax = A.x
	nzmax = A.nzmax
	nz = A.nz
	fmt.Printf("CSparse Version %d.%d.%d, %s.  %s\n", 3, 2, 0, "Sept 12, 2017", "Copyright (c) Timothy A. Davis, 2006-2016")
	if nz < 0 {
		fmt.Printf("%d-by-%d, nzmax: %d nnz: %d, 1-norm: %g\n", m, n, nzmax, Ap[n], cs_norm(A))
		for j := 0; j < n; j++ {
			fmt.Printf("    col %d : locations %d to %d\n", j, Ap[j], Ap[j+1]-1)
			for p = Ap[j]; p < Ap[j+1]; p++ {
				fmt.Printf("      %d : %v\n", Ai[p], func() float64 {
					if Ax != nil {
						return Ax[p]
					}
					return 1
				}())
				if brief && p > 20 {
					fmt.Printf("  ...\n")
					return true
				}
			}
		}
	} else {
		fmt.Printf("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz)
		for p = 0; p < nz; p++ {
			fmt.Printf("    %d %d : %v\n", Ai[p], Ap[p], func() float64 {
				if Ax != nil {
					return Ax[p]
				}
				return 1
			}())
			if brief && p > 20 {
				fmt.Printf("  ...\n")
				return true
			}
		}
	}
	return true
}

//
// // cs_pvec - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_pvec.c:3
// // x = b(p), for dense vectors x and b; p=NULL denotes identity
// func cs_pvec(p []noarch.PtrdiffT, b []float64, x []float64, n noarch.PtrdiffT) noarch.PtrdiffT {
// 	var k noarch.PtrdiffT
// 	if x == nil || b == nil {
// 		// check inputs
// 		return 0
// 	}
// 	for k = 0; k < n; k++ {
// 		x[k] = b[func() int32 {
// 			if p != nil {
// 				return int32(noarch.PtrdiffT(p[k]))
// 			}
// 			return int32(noarch.PtrdiffT(k))
// 		}()]
// 	}
// 	return 1
// }
//
// // cs_qr - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_qr.c:3
// // sparse QR factorization [V,beta,pinv,R] = qr (A)
// func cs_qr(A []cs, S []css) []csn {
// 	var Rx []float64
// 	var Vx []float64
// 	var Ax []float64
// 	var x []float64
// 	var Beta []float64
// 	var i noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var m noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var vnz noarch.PtrdiffT
// 	var p1 noarch.PtrdiffT
// 	var top noarch.PtrdiffT
// 	var m2 noarch.PtrdiffT
// 	var len noarch.PtrdiffT
// 	var col noarch.PtrdiffT
// 	var rnz noarch.PtrdiffT
// 	var s []noarch.PtrdiffT
// 	var leftmost []noarch.PtrdiffT
// 	var Ap []noarch.PtrdiffT
// 	var Ai []noarch.PtrdiffT
// 	var parent []noarch.PtrdiffT
// 	var Rp []noarch.PtrdiffT
// 	var Ri []noarch.PtrdiffT
// 	var Vp []noarch.PtrdiffT
// 	var Vi []noarch.PtrdiffT
// 	var w []noarch.PtrdiffT
// 	var pinv []noarch.PtrdiffT
// 	var q []noarch.PtrdiffT
// 	var R []cs
// 	var V []cs
// 	var N []csn
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) || S == nil {
// 		return nil
// 	}
// 	m = noarch.PtrdiffT(A[0].m)
// 	n = noarch.PtrdiffT(A[0].n)
// 	Ap = A[0].p
// 	Ai = A[0].i
// 	Ax = A[0].x
// 	q = S[0].q
// 	parent = S[0].parent
// 	pinv = S[0].pinv
// 	m2 = noarch.PtrdiffT(S[0].m2)
// 	vnz = noarch.PtrdiffT(S[0].lnz)
// 	rnz = noarch.PtrdiffT(S[0].unz)
// 	leftmost = S[0].leftmost
// 	// get csi workspace
// 	w = cs_malloc(m2+n, uint(0)).([]noarch.PtrdiffT)
// 	// get double workspace
// 	x = cs_malloc(noarch.PtrdiffT(m2), uint(8)).([]float64)
// 	// allocate result
// 	N = cs_calloc(1, uint(32)).([]csn)
// 	if w == nil || x == nil || N == nil {
// 		return (cs_ndone(N, nil, w, x, 0))
// 	}
// 	// s is size n
// 	s = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m2))*unsafe.Sizeof(w[0]))))[:]
// 	{
// 		// clear workspace x
// 		for k = 0; k < m2; k++ {
// 			x[k] = 0
// 		}
// 	}
// 	V = cs_spalloc(noarch.PtrdiffT(m2), noarch.PtrdiffT(n), noarch.PtrdiffT(vnz), 1, 0)
// 	// allocate result V
// 	N[0].L = V
// 	R = cs_spalloc(noarch.PtrdiffT(m2), noarch.PtrdiffT(n), noarch.PtrdiffT(rnz), 1, 0)
// 	// allocate result R
// 	N[0].U = R
// 	Beta = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
// 	// allocate result Beta
// 	N[0].B = Beta
// 	if R == nil || V == nil || Beta == nil {
// 		return (cs_ndone(N, nil, w, x, 0))
// 	}
// 	Rp = R[0].p
// 	Ri = R[0].i
// 	Rx = R[0].x
// 	Vp = V[0].p
// 	Vi = V[0].i
// 	Vx = V[0].x
// 	{
// 		// clear w, to mark nodes
// 		for i = 0; i < m2; i++ {
// 			w[i] = -1
// 		}
// 	}
// 	rnz = 0
// 	vnz = 0
// 	{
// 		// compute V and R
// 		for k = 0; k < n; k++ {
// 			// R(:,k) starts here
// 			Rp[k] = rnz
// 			p1 = vnz
// 			// V(:,k) starts here
// 			Vp[k] = p1
// 			// add V(k,k) to pattern of V
// 			w[k] = k
// 			Vi[func() noarch.PtrdiffT {
// 				defer func() {
// 					vnz++
// 				}()
// 				return vnz
// 			}()] = k
// 			top = n
// 			col = noarch.PtrdiffT(func() int32 {
// 				if q != nil {
// 					return int32(noarch.PtrdiffT(q[k]))
// 				}
// 				return int32(noarch.PtrdiffT(k))
// 			}() / 8)
// 			{
// 				// find R(:,k) pattern
// 				for p = Ap[col]; p < Ap[col+1]; p++ {
// 					// i = min(find(A(i,q)))
// 					i = leftmost[Ai[p]]
// 					{
// 						// traverse up to k
// 						for len = 0; w[i] != k; i = parent[i] {
// 							s[func() noarch.PtrdiffT {
// 								defer func() {
// 									len++
// 								}()
// 								return len
// 							}()] = i
// 							w[i] = k
// 						}
// 					}
// 					for len > 0 {
// 						// push path on stack
// 						s[func() noarch.PtrdiffT {
// 							top--
// 							return top
// 						}()] = s[func() noarch.PtrdiffT {
// 							len--
// 							return len
// 						}()]
// 					}
// 					// i = permuted row of A(:,col)
// 					i = pinv[Ai[p]]
// 					// x (i) = A(:,col)
// 					x[i] = Ax[p]
// 					if i > k && w[i] < k {
// 						// pattern of V(:,k) = x (k+1:m)
// 						// add i to pattern of V(:,k)
// 						Vi[func() noarch.PtrdiffT {
// 							defer func() {
// 								vnz++
// 							}()
// 							return vnz
// 						}()] = i
// 						w[i] = k
// 					}
// 				}
// 			}
// 			{
// 				// for each i in pattern of R(:,k)
// 				for p = top; p < n; p++ {
// 					// R(i,k) is nonzero
// 					i = s[p]
// 					// apply (V(i),Beta(i)) to x
// 					cs_happly(V, noarch.PtrdiffT(i), Beta[i], x)
// 					// R(i,k) = x(i)
// 					Ri[rnz] = i
// 					Rx[func() noarch.PtrdiffT {
// 						defer func() {
// 							rnz++
// 						}()
// 						return rnz
// 					}()] = x[i]
// 					x[i] = 0
// 					if parent[i] == k {
// 						vnz = cs_scatter(V, noarch.PtrdiffT(i), 0, w, nil, noarch.PtrdiffT(k), V, noarch.PtrdiffT(vnz))
// 					}
// 				}
// 			}
// 			{
// 				// gather V(:,k) = x
// 				for p = p1; p < vnz; p++ {
// 					Vx[p] = x[Vi[p]]
// 					x[Vi[p]] = 0
// 				}
// 			}
// 			// R(k,k) = norm (x)
// 			Ri[rnz] = k
// 			// [v,beta]=house(x)
// 			Rx[func() noarch.PtrdiffT {
// 				defer func() {
// 					rnz++
// 				}()
// 				return rnz
// 			}()] = cs_house((*(*[1000000000]float64)(unsafe.Pointer(uintptr(unsafe.Pointer(&Vx[0])) + (uintptr)(int(p1))*unsafe.Sizeof(Vx[0]))))[:], (*(*[1000000000]float64)(unsafe.Pointer(uintptr(unsafe.Pointer(&Beta[0])) + (uintptr)(int(k))*unsafe.Sizeof(Beta[0]))))[:], vnz-p1)
// 		}
// 	}
// 	// finalize R
// 	Rp[n] = rnz
// 	// finalize V
// 	Vp[n] = vnz
// 	// success
// 	return (cs_ndone(N, nil, w, x, 1))
// }
//
// // cs_qrsol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_qrsol.c:3
// // x=A\b where A can be rectangular; b overwritten with solution
// func cs_qrsol(order noarch.PtrdiffT, A []cs, b []float64) noarch.PtrdiffT {
// 	var x []float64
// 	var S []css
// 	var N []csn
// 	var AT []cs
// 	var k noarch.PtrdiffT
// 	var m noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var ok noarch.PtrdiffT
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) || b == nil {
// 		// check inputs
// 		return 0
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	m = noarch.PtrdiffT(A[0].m)
// 	if m >= n {
// 		// ordering and symbolic analysis
// 		S = cs_sqr(noarch.PtrdiffT(order), A, 1)
// 		// numeric QR factorization
// 		N = cs_qr(A, S)
// 		// get workspace
// 		x = cs_calloc(noarch.PtrdiffT(func() int32 {
// 			if S != nil {
// 				return int32(noarch.PtrdiffT(S[0].m2))
// 			}
// 			return 1
// 		}()/8), uint(8)).([]float64)
// 		ok = noarch.PtrdiffT(S != nil && N != nil && x != nil)
// 		if bool(ok) {
// 			// x(0:m-1) = b(p(0:m-1)
// 			cs_ipvec(S[0].pinv, b, x, m)
// 			{
// 				// apply Householder refl. to x
// 				for k = 0; k < n; k++ {
// 					cs_happly(N[0].L, noarch.PtrdiffT(k), N[0].B[k], x)
// 				}
// 			}
// 			// x = R\x
// 			cs_usolve(N[0].U, x)
// 			// b(q(0:n-1)) = x(0:n-1)
// 			cs_ipvec(S[0].q, x, b, noarch.PtrdiffT(n))
// 		}
// 	} else {
// 		// Ax=b is underdetermined
// 		AT = cs_transpose(A, 1)
// 		// ordering and symbolic analysis
// 		S = cs_sqr(noarch.PtrdiffT(order), AT, 1)
// 		// numeric QR factorization of A'
// 		N = cs_qr(AT, S)
// 		// get workspace
// 		x = cs_calloc(noarch.PtrdiffT(func() int32 {
// 			if S != nil {
// 				return int32(noarch.PtrdiffT(S[0].m2))
// 			}
// 			return 1
// 		}()/8), uint(8)).([]float64)
// 		ok = noarch.PtrdiffT(AT != nil && S != nil && N != nil && x != nil)
// 		if bool(ok) {
// 			// x(q(0:m-1)) = b(0:m-1)
// 			cs_pvec(S[0].q, b, x, m)
// 			// x = R'\x
// 			cs_utsolve(N[0].U, x)
// 			{
// 				// apply Householder refl. to x
// 				for k = m - 1; k >= 0; k-- {
// 					cs_happly(N[0].L, noarch.PtrdiffT(k), N[0].B[k], x)
// 				}
// 			}
// 			// b(0:n-1) = x(p(0:n-1))
// 			cs_pvec(S[0].pinv, x, b, noarch.PtrdiffT(n))
// 		}
// 	}
// 	cs_free(x)
// 	cs_sfree(S)
// 	cs_nfree(N)
// 	cs_spfree(AT) // TODO (KI) : remove
// 	return noarch.PtrdiffT((ok))
// }

// cs_randperm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_randperm.c:5
// return a random permutation vector, the identity perm, or p = n-1:-1:0.
// * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
// * p = random permutation.
func cs_randperm(n int, seed int) []int {
	if seed == 0 {
		// return p = NULL (identity)
		return nil
	}
	// allocate result
	p := make([]int, n)
	if p == nil {
		// out of memory
		return nil
	}
	for k := 0; k < n; k++ {
		p[k] = n - k - 1
	}
	if seed == -1 {
		// return reverse permutation
		return (p)
	}
	// get new random number seed
	rand.Seed(uint64(seed))
	for k := 0; k < n; k++ {
		// j = rand integer in range k to n-1
		j := k + (rand.Int() % (n - k))
		// swap p[k] and p[j]
		p[k], p[j] = p[j], p[k]
	}
	return p
}

// // cs_reach - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_reach.c:4
// // xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
// // * xi [n...2n-1] used as workspace
// func cs_reach(G []cs, B []cs, k noarch.PtrdiffT, xi []noarch.PtrdiffT, pinv []noarch.PtrdiffT) noarch.PtrdiffT {
// 	var p noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var top noarch.PtrdiffT
// 	var Bp []noarch.PtrdiffT
// 	var Bi []noarch.PtrdiffT
// 	var Gp []noarch.PtrdiffT
// 	if !(G != nil && noarch.PtrdiffT(G[0].nz) == -1) || !(B != nil && noarch.PtrdiffT(B[0].nz) == -1) || xi == nil {
// 		// check inputs
// 		return -1
// 	}
// 	n = noarch.PtrdiffT(G[0].n)
// 	Bp = B[0].p
// 	Bi = B[0].i
// 	Gp = G[0].p
// 	top = n
// 	for p = Bp[k]; p < Bp[k+1]; p++ {
// 		if !(Gp[Bi[p]] < 0) {
// 			// start a dfs at unmarked node i
// 			top = cs_dfs(noarch.PtrdiffT(Bi[p]), G, noarch.PtrdiffT(top), xi, (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&xi[0])) + (uintptr)(int(n))*unsafe.Sizeof(xi[0]))))[:], pinv)
// 		}
// 	}
// 	{
// 		// restore G
// 		for p = top; p < n; p++ {
// 			Gp[xi[p]] = -noarch.PtrdiffT((Gp[xi[p]])) - 2
// 		}
// 	}
// 	return noarch.PtrdiffT((top))
// }

// cs_scatter - x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse
func cs_scatter(A *cs, j int, beta float64, w []int, x []float64, mark int, C *cs, nz int) int {
	var i int
	var p int
	var Ap []int
	var Ai []int
	var Ci []int
	var Ax []float64
	if !(A != nil && A.nz == -1) || w == nil || !(C != nil && C.nz == -1) {
		// check inputs
		return -1
	}
	Ap = A.p
	Ai = A.i
	Ax = A.x
	Ci = C.i
	for p = Ap[j]; p < Ap[j+1]; p++ {
		// A(i,j) is nonzero
		i = Ai[p]
		if w[i] < mark {
			// i is new entry in column j
			w[i] = mark
			// add i to pattern of C(:,j)
			Ci[func() int {
				defer func() {
					nz++
				}()
				return nz
			}()] = i
			if x != nil {
				// x(i) = beta*A(i,j)
				x[i] = beta * Ax[p]
			}
		} else if x != nil {
			// i exists in C(:,j) already
			x[i] += beta * Ax[p]
		}
	}
	return nz
}

// cs_scc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_scc.c:3
// find the strongly connected components of a square matrix
// matrix A temporarily modified, then restored
func cs_scc(A *cs) *csd {
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n := A.n
	Ap := A.p
	// allocate result
	D := cs_dalloc(n, 0)
	// AT = A'
	AT := cs_transpose(A, false)
	// get workspace
	xi := make([]int, 2*n+1)
	if D == nil || AT == nil || xi == nil {
		return cs_ddone(D, AT, xi, false)
	}
	Blk := xi
	pstack := xi[n:]
	rcopy := pstack
	p := D.p
	r := D.r
	ATp := AT.p
	top := n

	// first dfs(A) to find finish times (xi)
	for i := 0; i < n; i++ {
		if !(Ap[i] < 0) {
			top = cs_dfs(i, A, top, xi, pstack, nil)
		}
	}

	// restore A; unmark all nodes
	for i := 0; i < n; i++ {
		Ap[i] = -Ap[i] - 2
	}

	top = n
	nb := n

	// dfs(A') to find strongly connnected comp
	for k := 0; k < n; k++ {
		// get i in reverse order of finish times
		i := xi[k]
		if ATp[i] < 0 {
			// skip node i if already ordered
			continue
		}
		// node i is the start of a component in p
		r[func() int {
			defer func() {
				nb--
			}()
			return nb
		}()] = top
		top = cs_dfs(i, AT, top, p, pstack, nil)
	}

	// first block starts at zero; shift r up
	r[nb] = 0
	for k := nb; k <= n; k++ {
		r[k-nb] = r[k]
	}
	nb = n - nb
	// nb = # of strongly connected components
	D.nb = nb

	// sort each block in natural order
	for b := 0; b < nb; b++ {
		for k := r[b]; k < r[b+1]; k++ {
			Blk[p[k]] = b
		}
	}

	for b := 0; b <= nb; b++ {
		rcopy[b] = r[b]
	}
	for i := 0; i < n; i++ {
		p[func() int {
			tempVar := &rcopy[Blk[i]]
			defer func() {
				*tempVar++
			}()
			return *tempVar
		}()] = i
	}
	return (cs_ddone(D, AT, xi, true))
}

//
// // cs_schol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_schol.c:3
// // ordering and symbolic analysis for a Cholesky factorization
// func cs_schol(order noarch.PtrdiffT, A []cs) []css {
// 	var n noarch.PtrdiffT
// 	var c []noarch.PtrdiffT
// 	var post []noarch.PtrdiffT
// 	var P []noarch.PtrdiffT
// 	var C []cs
// 	var S []css
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) {
// 		// check inputs
// 		return nil
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	// allocate result S
// 	S = cs_calloc(1, uint(0)).([]css)
// 	if S == nil {
// 		// out of memory
// 		return nil
// 	}
// 	// P = amd(A+A'), or natural
// 	P = cs_amd(noarch.PtrdiffT(order), A)
// 	// find inverse permutation
// 	S[0].pinv = cs_pinv(P, noarch.PtrdiffT(n))
// 	cs_free(P)
// 	if bool(order) && S[0].pinv == nil {
// 		return (cs_sfree(S))
// 	}
// 	// C = spones(triu(A(P,P)))
// 	C = cs_symperm(A, S[0].pinv, 0)
// 	// find etree of C
// 	S[0].parent = cs_etree(C, 0)
// 	// postorder the etree
// 	post = cs_post(S[0].parent, noarch.PtrdiffT(n))
// 	// find column counts of chol(C)
// 	c = cs_counts(C, S[0].parent, post, 0)
// 	cs_free(post)
// 	cs_spfree(C) // TODO (KI) : remove
// 	// allocate result S->cp
// 	S[0].cp = cs_malloc(n+1, uint(0)).([]noarch.PtrdiffT)
// 	S[0].lnz = cs_cumsum(S[0].cp, c, noarch.PtrdiffT(n))
// 	// find column pointers for L
// 	S[0].unz = S[0].lnz
// 	cs_free(c)
// 	return (func() []css {
// 		if S[0].lnz >= 0 {
// 			return S
// 		}
// 		return cs_sfree(S)
// 	}())
// }
//
// // cs_spsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_spsolve.c:3
// // solve Gx=b(:,k), where G is either upper (lo=0) or lower (lo=1) triangular
// func cs_spsolve(G []cs, B []cs, k noarch.PtrdiffT, xi []noarch.PtrdiffT, x []float64, pinv []noarch.PtrdiffT, lo noarch.PtrdiffT) noarch.PtrdiffT {
// 	var j noarch.PtrdiffT
// 	var J noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var q noarch.PtrdiffT
// 	var px noarch.PtrdiffT
// 	var top noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var Gp []noarch.PtrdiffT
// 	var Gi []noarch.PtrdiffT
// 	var Bp []noarch.PtrdiffT
// 	var Bi []noarch.PtrdiffT
// 	var Gx []float64
// 	var Bx []float64
// 	if !(G != nil && noarch.PtrdiffT(G[0].nz) == -1) || !(B != nil && noarch.PtrdiffT(B[0].nz) == -1) || xi == nil || x == nil {
// 		return -1
// 	}
// 	Gp = G[0].p
// 	Gi = G[0].i
// 	Gx = G[0].x
// 	n = noarch.PtrdiffT(G[0].n)
// 	Bp = B[0].p
// 	Bi = B[0].i
// 	Bx = B[0].x
// 	// xi[top..n-1]=Reach(B(:,k))
// 	top = cs_reach(G, B, noarch.PtrdiffT(k), xi, pinv)
// 	{
// 		// clear x
// 		for p = top; p < n; p++ {
// 			x[xi[p]] = 0
// 		}
// 	}
// 	{
// 		// scatter B
// 		for p = Bp[k]; p < Bp[k+1]; p++ {
// 			x[Bi[p]] = Bx[p]
// 		}
// 	}
// 	for px = top; px < n; px++ {
// 		// x(j) is nonzero
// 		j = xi[px]
// 		// j maps to col J of G
// 		J = noarch.PtrdiffT(func() int32 {
// 			if pinv != nil {
// 				return int32(noarch.PtrdiffT((pinv[j])))
// 			}
// 			return int32(noarch.PtrdiffT(j))
// 		}() / 8)
// 		if J < 0 {
// 			// column J is empty
// 			continue
// 		}
// 		// x(j) /= G(j,j)
// 		x[j] /= Gx[func() int32 {
// 			if bool(noarch.PtrdiffT(lo)) {
// 				return int32(noarch.PtrdiffT((Gp[J])))
// 			}
// 			return (int32(Gp[J+1] - 1))
// 		}()]
// 		// lo: L(j,j) 1st entry
// 		p = noarch.PtrdiffT(func() int32 {
// 			if bool(noarch.PtrdiffT(lo)) {
// 				return (int32(Gp[J] + 1))
// 			}
// 			return int32(noarch.PtrdiffT((Gp[J])))
// 		}() / 8)
// 		// up: U(j,j) last entry
// 		q = noarch.PtrdiffT(func() int32 {
// 			if bool(noarch.PtrdiffT(lo)) {
// 				return int32(noarch.PtrdiffT((Gp[J+1])))
// 			}
// 			return (int32(Gp[J+1] - 1))
// 		}() / 8)
// 		for ; p < q; p++ {
// 			// x(i) -= G(i,j) * x(j)
// 			x[Gi[p]] -= Gx[p] * x[j]
// 		}
// 	}
// 	// return top of stack
// 	return noarch.PtrdiffT((top))
// }
//
// // cs_vcount - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_sqr.c:3
// // compute nnz(V) = S->lnz, S->pinv, S->leftmost, S->m2 from A and S->parent
// func cs_vcount(A []cs, S []css) noarch.PtrdiffT {
// 	var i noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var pa noarch.PtrdiffT
// 	var n noarch.PtrdiffT = noarch.PtrdiffT(A[0].n)
// 	var m noarch.PtrdiffT = noarch.PtrdiffT(A[0].m)
// 	var Ap []noarch.PtrdiffT = A[0].p
// 	var Ai []noarch.PtrdiffT = A[0].i
// 	var next []noarch.PtrdiffT
// 	var head []noarch.PtrdiffT
// 	var tail []noarch.PtrdiffT
// 	var nque []noarch.PtrdiffT
// 	var pinv []noarch.PtrdiffT
// 	var leftmost []noarch.PtrdiffT
// 	var w []noarch.PtrdiffT
// 	var parent []noarch.PtrdiffT = S[0].parent
// 	pinv = cs_malloc(m+n, uint(0)).([]noarch.PtrdiffT)
// 	// allocate pinv,
// 	S[0].pinv = pinv
// 	leftmost = cs_malloc(m, uint(0)).([]noarch.PtrdiffT)
// 	// and leftmost
// 	S[0].leftmost = leftmost
// 	// get workspace
// 	w = cs_malloc(m+noarch.PtrdiffT(3*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
// 	if pinv == nil || w == nil || leftmost == nil {
// 		// pinv and leftmost freed later
// 		cs_free(w)
// 		// out of memory
// 		return 0
// 	}
// 	next = w
// 	head = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:]
// 	tail = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&(*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0])) + (uintptr)(int(n))*unsafe.Sizeof((*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0]))))[:]
// 	nque = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&(*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof((*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0]))))[:]
// 	{
// 		// queue k is empty
// 		for k = 0; k < n; k++ {
// 			head[k] = -1
// 		}
// 	}
// 	for k = 0; k < n; k++ {
// 		tail[k] = -1
// 	}
// 	for k = 0; k < n; k++ {
// 		nque[k] = 0
// 	}
// 	for i = 0; i < m; i++ {
// 		leftmost[i] = -1
// 	}
// 	for k = n - 1; k >= 0; k-- {
// 		for p = Ap[k]; p < Ap[k+1]; p++ {
// 			// leftmost[i] = min(find(A(i,:)))
// 			leftmost[Ai[p]] = k
// 		}
// 	}
// 	{
// 		// scan rows in reverse order
// 		for i = m - 1; i >= 0; i-- {
// 			// row i is not yet ordered
// 			pinv[i] = -1
// 			k = leftmost[i]
// 			if k == -1 {
// 				// row i is empty
// 				continue
// 			}
// 			if func() noarch.PtrdiffT {
// 				tempVar := &nque[k]
// 				defer func() {
// 					*tempVar++
// 				}()
// 				return *tempVar
// 			}() == 0 {
// 				// first row in queue k
// 				tail[k] = i
// 			}
// 			// put i at head of queue k
// 			next[i] = head[k]
// 			head[k] = i
// 		}
// 	}
// 	S[0].lnz = 0
// 	S[0].m2 = m
// 	{
// 		// find row permutation and nnz(V)
// 		for k = 0; k < n; k++ {
// 			// remove row i from queue k
// 			i = head[k]
// 			// count V(k,k) as nonzero
// 			S[0].lnz++
// 			if i < 0 {
// 				// add a fictitious row
// 				i = func() noarch.PtrdiffT {
// 					tempVar := &S[0].m2
// 					defer func() {
// 						*tempVar++
// 					}()
// 					return *tempVar
// 				}()
// 			}
// 			// associate row i with V(:,k)
// 			pinv[i] = k
// 			if func() noarch.PtrdiffT {
// 				tempVar := &nque[k]
// 				*tempVar--
// 				return *tempVar
// 			}() <= 0 {
// 				// skip if V(k+1:m,k) is empty
// 				continue
// 			}
// 			// nque [k] is nnz (V(k+1:m,k))
// 			S[0].lnz += float64(nque[k])
// 			if (func() noarch.PtrdiffT {
// 				pa = parent[k]
// 				return pa
// 			}()) != -1 {
// 				if nque[pa] == 0 {
// 					// move all rows to parent of k
// 					tail[pa] = tail[k]
// 				}
// 				next[tail[k]] = head[pa]
// 				head[pa] = next[i]
// 				nque[pa] += nque[k]
// 			}
// 		}
// 	}
// 	for i = 0; i < m; i++ {
// 		if pinv[i] < 0 {
// 			pinv[i] = func() noarch.PtrdiffT {
// 				defer func() {
// 					k++
// 				}()
// 				return k
// 			}()
// 		}
// 	}
// 	cs_free(w)
// 	return 1
// }
//
// // cs_sqr - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_sqr.c:60
// // symbolic ordering and analysis for QR or LU
// func cs_sqr(order noarch.PtrdiffT, A []cs, qr noarch.PtrdiffT) []css {
// 	var n noarch.PtrdiffT
// 	var k noarch.PtrdiffT
// 	var ok noarch.PtrdiffT = 1
// 	var post []noarch.PtrdiffT
// 	var S []css
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) {
// 		// check inputs
// 		return nil
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	// allocate result S
// 	S = cs_calloc(1, uint(0)).([]css)
// 	if S == nil {
// 		// out of memory
// 		return nil
// 	}
// 	// fill-reducing ordering
// 	S[0].q = cs_amd(noarch.PtrdiffT(order), A)
// 	if bool(order) && S[0].q == nil {
// 		return (cs_sfree(S))
// 	}
// 	if bool(noarch.PtrdiffT(qr)) {
// 		var C []cs = func() []cs {
// 			if bool(noarch.PtrdiffT(order)) {
// 				return cs_permute(A, nil, S[0].q, 0)
// 			}
// 			return (A)
// 		}()
// 		// QR symbolic analysis
// 		// etree of C'*C, where C=A(:,q)
// 		S[0].parent = cs_etree(C, 1)
// 		post = cs_post(S[0].parent, noarch.PtrdiffT(n))
// 		// col counts chol(C'*C)
// 		S[0].cp = cs_counts(C, S[0].parent, post, 1)
// 		cs_free(post)
// 		ok = noarch.PtrdiffT(C != nil && S[0].parent != nil && S[0].cp != nil && bool(cs_vcount(C, S)))
// 		if bool(ok) {
// 			S[0].unz = 0
// 			k = 0
// 			for k = 0; k < n; k++ {
// 				S[0].unz += float64(S[0].cp[k])
// 			}
// 		}
// 		if bool(noarch.PtrdiffT(order)) {
// 			cs_spfree(C) // TODO (KI) : remove
// 		}
// 	} else {
// 		// for LU factorization only,
// 		S[0].unz = float64(4*int32(A[0].p[n]) + int32(n))
// 		// guess nnz(L) and nnz(U)
// 		S[0].lnz = S[0].unz
// 	}
// 	// return result S
// 	return (func() []css {
// 		if bool(ok) {
// 			return S
// 		}
// 		return cs_sfree(S)
// 	}())
// }
//
// // cs_symperm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_symperm.c:3
// // C = A(p,p) where A and C are symmetric the upper part stored; pinv not p
// func cs_symperm(A []cs, pinv []noarch.PtrdiffT, values noarch.PtrdiffT) []cs {
// 	var i noarch.PtrdiffT
// 	var j noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var q noarch.PtrdiffT
// 	var i2 noarch.PtrdiffT
// 	var j2 noarch.PtrdiffT
// 	var n noarch.PtrdiffT
// 	var Ap []noarch.PtrdiffT
// 	var Ai []noarch.PtrdiffT
// 	var Cp []noarch.PtrdiffT
// 	var Ci []noarch.PtrdiffT
// 	var w []noarch.PtrdiffT
// 	var Cx []float64
// 	var Ax []float64
// 	var C []cs
// 	if !(A != nil && noarch.PtrdiffT(A[0].nz) == -1) {
// 		// check inputs
// 		return nil
// 	}
// 	n = noarch.PtrdiffT(A[0].n)
// 	Ap = A[0].p
// 	Ai = A[0].i
// 	Ax = A[0].x
// 	// alloc result
// 	C = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(Ap[n]), noarch.PtrdiffT(bool(values) && Ax != nil), 0)
// 	// get workspace
// 	w = cs_calloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
// 	if C == nil || w == nil {
// 		// out of memory
// 		return (cs_done(C, w, nil, 0))
// 	}
// 	Cp = C[0].p
// 	Ci = C[0].i
// 	Cx = C[0].x
// 	{
// 		// count entries in each column of C
// 		for j = 0; j < n; j++ {
// 			// column j of A is column j2 of C
// 			j2 = noarch.PtrdiffT(func() int32 {
// 				if pinv != nil {
// 					return int32(noarch.PtrdiffT(pinv[j]))
// 				}
// 				return int32(noarch.PtrdiffT(j))
// 			}() / 8)
// 			for p = Ap[j]; p < Ap[j+1]; p++ {
// 				i = Ai[p]
// 				if i > j {
// 					// skip lower triangular part of A
// 					continue
// 				}
// 				// row i of A is row i2 of C
// 				i2 = noarch.PtrdiffT(func() int32 {
// 					if pinv != nil {
// 						return int32(noarch.PtrdiffT(pinv[i]))
// 					}
// 					return int32(noarch.PtrdiffT(i))
// 				}() / 8)
// 				// column count of C
// 				w[func() int32 {
// 					if i2 > j2 {
// 						return int32(noarch.PtrdiffT((i2)))
// 					}
// 					return int32(noarch.PtrdiffT((j2)))
// 				}()]++
// 			}
// 		}
// 	}
// 	// compute column pointers of C
// 	cs_cumsum(Cp, w, noarch.PtrdiffT(n))
// 	for j = 0; j < n; j++ {
// 		// column j of A is column j2 of C
// 		j2 = noarch.PtrdiffT(func() int32 {
// 			if pinv != nil {
// 				return int32(noarch.PtrdiffT(pinv[j]))
// 			}
// 			return int32(noarch.PtrdiffT(j))
// 		}() / 8)
// 		for p = Ap[j]; p < Ap[j+1]; p++ {
// 			i = Ai[p]
// 			if i > j {
// 				// skip lower triangular part of A
// 				continue
// 			}
// 			// row i of A is row i2 of C
// 			i2 = noarch.PtrdiffT(func() int32 {
// 				if pinv != nil {
// 					return int32(noarch.PtrdiffT(pinv[i]))
// 				}
// 				return int32(noarch.PtrdiffT(i))
// 			}() / 8)
// 			Ci[(func() noarch.PtrdiffT {
// 				q = func() noarch.PtrdiffT {
// 					tempVar := &w[func() int32 {
// 						if i2 > j2 {
// 							return int32(noarch.PtrdiffT((i2)))
// 						}
// 						return int32(noarch.PtrdiffT((j2)))
// 					}()]
// 					defer func() {
// 						*tempVar++
// 					}()
// 					return *tempVar
// 				}()
// 				return q
// 			}())] = noarch.PtrdiffT(func() int32 {
// 				if i2 < j2 {
// 					return int32(noarch.PtrdiffT((i2)))
// 				}
// 				return int32(noarch.PtrdiffT((j2)))
// 			}() / 8)
// 			if Cx != nil {
// 				Cx[q] = Ax[p]
// 			}
// 		}
// 	}
// 	// success; free workspace, return C
// 	return (cs_done(C, w, nil, 1))
// }
//
// // cs_tdfs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_tdfs.c:3
// // depth-first search and postorder of a tree rooted at node j
// func cs_tdfs(j, k int, head []int, next []int, post []int, stack []int) noarch.PtrdiffT {
// 	var i noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var top noarch.PtrdiffT
// 	if head == nil || next == nil || post == nil || stack == nil {
// 		// check inputs
// 		return -1
// 	}
// 	// place j on the stack
// 	stack[0] = j
// 	for top >= 0 {
// 		// while (stack is not empty)
// 		// p = top of stack
// 		p = stack[top]
// 		// i = youngest child of p
// 		i = head[p]
// 		if i == -1 {
// 			// p has no unordered children left
// 			top--
// 			// node p is the kth postordered node
// 			post[func() noarch.PtrdiffT {
// 				defer func() {
// 					k++
// 				}()
// 				return k
// 			}()] = p
// 		} else {
// 			// remove i from children of p
// 			head[p] = next[i]
// 			// start dfs on child node i
// 			stack[func() noarch.PtrdiffT {
// 				top++
// 				return top
// 			}()] = i
// 		}
// 	}
// 	return noarch.PtrdiffT((k))
// }

// cs_transpose - C = A'
func cs_transpose(A *cs, values bool) *cs {
	var p int
	var q int
	var j int
	var Cp []int
	var Ci []int
	var n int
	var m int
	var Ap []int
	var Ai []int
	var w []int
	var Cx []float64
	var Ax []float64
	var C *cs
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	m = A.m
	n = A.n
	Ap = A.p
	Ai = A.i
	Ax = A.x
	// allocate result
	C = cs_spalloc(n, m, Ap[n], values && Ax != nil, false)
	// get workspace
	w = make([]int, m)
	if C == nil || w == nil {
		// out of memory
		return cs_done(C, w, nil, false)
	}
	Cp = C.p
	Ci = C.i
	Cx = C.x

	// row counts
	for p = 0; p < Ap[n]; p++ {
		w[Ai[p]]++
	}

	// row pointers
	cs_cumsum(Cp, w, m)
	for j = 0; j < n; j++ {
		for p = Ap[j]; p < Ap[j+1]; p++ {
			// place A(i,j) as entry C(j,i)
			Ci[(func() int {
				q = func() int {
					tempVar := &w[Ai[p]]
					defer func() {
						*tempVar++
					}()
					return *tempVar
				}()
				return q
			}())] = j
			if Cx != nil {
				Cx[q] = Ax[p]
			}
		}
	}
	// success; free w and return C
	return cs_done(C, w, nil, true)
}

// // cs_updown - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_updown.c:3
// // sparse Cholesky update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1)
// func cs_updown(L []cs, sigma noarch.PtrdiffT, C []cs, parent []noarch.PtrdiffT) noarch.PtrdiffT {
// 	var n noarch.PtrdiffT
// 	var p noarch.PtrdiffT
// 	var f noarch.PtrdiffT
// 	var j noarch.PtrdiffT
// 	var Lp []noarch.PtrdiffT
// 	var Li []noarch.PtrdiffT
// 	var Cp []noarch.PtrdiffT
// 	var Ci []noarch.PtrdiffT
// 	var Lx []float64
// 	var Cx []float64
// 	var alpha float64
// 	var beta float64 = 1
// 	var delta float64
// 	var gamma float64
// 	var w1 float64
// 	var w2 float64
// 	var w []float64
// 	var beta2 float64 = 1
// 	if !(L != nil && noarch.PtrdiffT(L[0].nz) == -1) || !(C != nil && noarch.PtrdiffT(C[0].nz) == -1) || parent == nil {
// 		// check inputs
// 		return 0
// 	}
// 	Lp = L[0].p
// 	Li = L[0].i
// 	Lx = L[0].x
// 	n = noarch.PtrdiffT(L[0].n)
// 	Cp = C[0].p
// 	Ci = C[0].i
// 	Cx = C[0].x
// 	if (func() noarch.PtrdiffT {
// 		p = Cp[0]
// 		return p
// 	}()) >= Cp[1] {
// 		// return if C empty
// 		return 1
// 	}
// 	// get workspace
// 	w = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
// 	if w == nil {
// 		// out of memory
// 		return 0
// 	}
// 	f = Ci[p]
// 	for ; p < Cp[1]; p++ {
// 		// f = min (find (C))
// 		f = noarch.PtrdiffT(func() int32 {
// 			if f < Ci[p] {
// 				return int32(f)
// 			}
// 			return int32(noarch.PtrdiffT((Ci[p])))
// 		}() / 8)
// 	}
// 	{
// 		// clear workspace w
// 		for j = f; j != -1; j = parent[j] {
// 			w[j] = 0
// 		}
// 	}
// 	{
// 		// w = C
// 		for p = Cp[0]; p < Cp[1]; p++ {
// 			w[Ci[p]] = Cx[p]
// 		}
// 	}
// 	{
// 		// walk path f up to root
// 		for j = f; j != -1; j = parent[j] {
// 			p = Lp[j]
// 			// alpha = w(j) / L(j,j)
// 			alpha = w[j] / Lx[p]
// 			beta2 = beta*beta + float64(sigma)*alpha*alpha
// 			if beta2 <= 0 {
// 				// not positive definite
// 				break
// 			}
// 			beta2 = math.Sqrt(beta2)
// 			delta = func() float64 {
// 				if sigma > 0 {
// 					return (beta / beta2)
// 				}
// 				return (beta2 / beta)
// 			}()
// 			gamma = float64(sigma) * alpha / (beta2 * beta)
// 			Lx[p] = delta*Lx[p] + func() float64 {
// 				if sigma > 0 {
// 					return (gamma * w[j])
// 				}
// 				return 0
// 			}()
// 			beta = beta2
// 			for p += 1; p < Lp[j+1]; p++ {
// 				w1 = w[Li[p]]
// 				w2 = w1 - alpha*Lx[p]
// 				w[Li[p]] = w2
// 				Lx[p] = delta*Lx[p] + gamma*func() float64 {
// 					if sigma > 0 {
// 						return w1
// 					}
// 					return w2
// 				}()
// 			}
// 		}
// 	}
// 	cs_free(w)
// 	return noarch.PtrdiffT((beta2 > 0))
// }
//
// // cs_usolve - solve Ux=b where x and b are dense.  x=b on input, solution on output.
// func cs_usolve(U *cs, x []float64) bool {
// 	if !(U != nil && U.nz == -1) || x == nil {
// 		// check inputs
// 		return false
// 	}
// 	var (
// 		n  = U.n
// 		Up = U.p
// 		Ui = U.i
// 		Ux = U.x
// 	)
// 	for j := n - 1; j >= 0; j-- {
// 		x[j] /= Ux[Up[j+1]-1]
// 		for p := Up[j]; p < Up[j+1]-1; p++ {
// 			x[Ui[p]] -= Ux[p] * x[j]
// 		}
// 	}
// 	//
// 	// TODO (KI) : probably var `x` is return var
// 	//
// 	return true
// }

// cs_spalloc - allocate a sparse matrix (triplet form or compressed-column form)
func cs_spalloc(m, n, nzmax int, values, triplet bool) *cs {
	A := new(cs)
	if A == nil {
		// allocate the cs struct
		// out of memory
		return nil
	}
	// define dimensions and nzmax
	A.m = m
	A.n = n
	nzmax = func() int {
		if nzmax > 1 {
			return nzmax
		}
		return 1
	}()
	A.nzmax = nzmax
	// allocate triplet or comp.col
	A.nz = func() int {
		if triplet {
			return 0
		}
		return -1
	}()
	if triplet {
		A.p = make([]int, nzmax)
	} else {
		A.p = make([]int, n+1)
	}
	A.i = make([]int, nzmax) // cs_malloc(nzmax, uint(0)).([]noarch.PtrdiffT)
	A.x = make([]float64, nzmax)
	// func() interface{} {
	// 	if values {
	// 		return cs_malloc(nzmax, uint(8))
	// 	}
	// 	return nil
	// }().([]float64)
	return (func() *cs {
		if A.p == nil || A.i == nil || values && A.x == nil {
			return nil //cs_spfree(A) // TODO (KI) : remove
		}
		return A
	}())
}

// cs_sprealloc - change the max # of entries sparse matrix
func cs_sprealloc(A *cs, nzmax int) bool {
	var oki bool
	var okj bool = true
	var okx bool = true
	if A == nil {
		return false
	}
	if nzmax <= 0 {
		nzmax = func() int {
			if A != nil && A.nz == -1 {
				return A.p[A.n]
			}
			return A.nz
		}()
	}
	nzmax = func() int {
		if nzmax > 1 {
			return nzmax
		}
		return 1
	}()

	A.i = cs_realloc(A.i, nzmax, &oki).([]int)
	if A != nil && A.nz >= 0 {
		A.p = cs_realloc(A.p, nzmax, &okj).([]int)
	}
	if A.x != nil {
		A.x = cs_realloc(A.x, nzmax, &okx).([]float64)
	}
	ok := oki && okj && okx
	if ok {
		A.nzmax = nzmax
	}
	return ok
}

// cs_spfree - free a sparse matrix
func cs_spfree(A *cs) *cs {
	// free the cs struct and return NULL
	return nil
}

// // cs_nfree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:43
// // free a numeric factorization
// func cs_nfree(N []csn) []csn {
// 	if N == nil {
// 		// do nothing if N already NULL
// 		return nil
// 	}
// 	cs_spfree(N[0].L) // TODO (KI) : remove
// 	cs_spfree(N[0].U) // TODO (KI) : remove
// 	cs_free(N[0].pinv)
// 	cs_free(N[0].B)
// 	// free the csn struct and return NULL
// 	return (cs_free(N).([]csn))
// }
//
// // cs_sfree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:54
// // free a symbolic factorization
// func cs_sfree(S []css) []css {
// 	if S == nil {
// 		// do nothing if S already NULL
// 		return nil
// 	}
// 	cs_free(S[0].pinv)
// 	cs_free(S[0].q)
// 	cs_free(S[0].parent)
// 	cs_free(S[0].cp)
// 	cs_free(S[0].leftmost)
// 	// free the css struct and return NULL
// 	return (cs_free(S).([]css))
// }

// cs_dalloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:66
// allocate a cs_dmperm or cs_scc result
func cs_dalloc(m, n int) *csd {
	D := new(csd)
	if D == nil {
		return nil
	}
	D.p = make([]int, m)
	D.r = make([]int, m+6)
	D.q = make([]int, n)
	D.s = make([]int, n+6)
	return (func() *csd {
		if D.p == nil || D.r == nil || D.q == nil || D.s == nil {
			return cs_dfree(D)
		}
		return D
	}())
}

// cs_dfree - free a cs_dmperm or cs_scc result
func cs_dfree(D *csd) *csd {
	// free the csd struct and return NULL
	return nil
}

// cs_done - free workspace and return a sparse matrix result
func cs_done(C *cs, w []int, x []float64, ok bool) *cs {
	//
	// TODO(KI): remove w,x
	//

	// return result if OK, else free it
	if ok {
		return C
	}
	return nil
}

// cs_idone - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:98
// free workspace and return csi array result
func cs_idone(p []int, C *cs, w interface{}, ok bool) []int {
	//
	// TODO (KI) : remove C, w
	//

	// return result, or free it
	if ok {
		return p
	}
	return nil
}

// // cs_ndone - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:106
// // free workspace and return a numeric factorization (Cholesky, LU, or QR)
// func cs_ndone(N *csn, C *cs, w interface{}, x interface{}, ok bool) *csn {
// 	// return result if OK, else free it
// 	if ok {
// 		return N
// 	}
// 	return nil
// }

// cs_ddone - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:115
// free workspace and return a csd result
func cs_ddone(D *csd, C *cs, w interface{}, ok bool) *csd {
	// return result if OK, else free it
	if ok {
		return D
	}
	return nil
}

// // cs_utsolve - solve U'x=b where x and b are dense.  x=b on input, solution on output.
// func cs_utsolve(U *cs, x []float64) bool {
// 	if !(U != nil && U.nz == -1) || x == nil {
// 		// check inputs
// 		return false
// 	}
// 	var (
// 		n  = U.n
// 		Up = U.p
// 		Ui = U.i
// 		Ux = U.x
// 	)
// 	for j := 0; j < n; j++ {
// 		for p := Up[j]; p < Up[j+1]-1; p++ {
// 			x[j] -= Ux[p] * x[Ui[p]]
// 		}
// 		x[j] /= Ux[Up[j+1]-1]
// 	}
// 	//
// 	// TODO (KI) : probably variable `x` is return var
// 	//
// 	return true
// }
