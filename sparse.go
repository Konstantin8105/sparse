package main

import (
	"fmt"
	"github.com/Konstantin8105/c4go/noarch"
	"math"
	"math/rand"
	"unsafe"
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
	p  *int   // size m, row permutation
	q  *int   // size n, column permutation
	r  *int   // size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q)
	s  *int   // size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q)
	nb int    // # of blocks in fine dmperm decomposition
	rr [5]int // coarse row decomposition
	cc [5]int // coarse column decomposition
}

// TODO(KI) : remove comments like "transpiled function from ..."

// cs_add - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_add.c:3
// C = alpha*A + beta*B
func cs_add(A *cs, B *cs, alpha float64, beta float64) *cs {
	var p int
	var j int
	var nz int
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
		w = new(int)
	)
	values := (A.x != nil && Bx != nil)
	// get workspace
	if values {
		x = make([]float64, m)
	}
	// allocate result
	C := cs_spalloc(m, n, anz+bnz, values, 0)
	if C == nil || w == nil || values && x == nil {
		return cs_done(C, w, x, false)
	}
	var (
		Cp = C.p
		Ci = C.i
		Cx = C.x
	)
	for j = 0; j < n; j++ {
		// column j of C starts here
		Cp[j] = nz
		// alpha*A(:,j)
		nz = cs_scatter(A, j, alpha, w, x, j+1, C, nz)
		// beta*B(:,j)
		nz = cs_scatter(B, j, beta, w, x, j+1, C, nz)
		if values {
			for p = Cp[j]; p < nz; p++ {
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

// cs_wclear - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_amd.c:3
// clear w
func cs_wclear(mark noarch.PtrdiffT, lemax noarch.PtrdiffT, w []noarch.PtrdiffT, n noarch.PtrdiffT) noarch.PtrdiffT {
	var k noarch.PtrdiffT
	if mark < noarch.PtrdiffT(2/8) || mark+lemax < noarch.PtrdiffT(0/8) {
		for k = 0; k < n; k++ {
			if w[k] != noarch.PtrdiffT(0/8) {
				w[k] = 1
			}
		}
		mark = 2
	}
	// at this point, w [0..n-1] < mark holds
	return noarch.PtrdiffT((mark))
}

// cs_diag - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_amd.c:15
// keep off-diagonal entries; drop diagonal entries
func cs_diag(i noarch.PtrdiffT, j noarch.PtrdiffT, aij float64, other interface{}) noarch.PtrdiffT {
	return noarch.PtrdiffT((i != j))
}

// cs_amd - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_amd.c:18
// p = amd(A+A') if symmetric is true, or amd(A'A) otherwise
// order 0:natural, 1:Chol, 2:LU, 3:QR
func cs_amd(order noarch.PtrdiffT, A []cs) []noarch.PtrdiffT {
	var C []cs
	var A2 []cs
	var AT []cs
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var last []noarch.PtrdiffT
	var W []noarch.PtrdiffT
	var len []noarch.PtrdiffT
	var nv []noarch.PtrdiffT
	var next []noarch.PtrdiffT
	var P []noarch.PtrdiffT
	var head []noarch.PtrdiffT
	var elen []noarch.PtrdiffT
	var degree []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var hhead []noarch.PtrdiffT
	var ATp []noarch.PtrdiffT
	var ATi []noarch.PtrdiffT
	var d noarch.PtrdiffT
	var dk noarch.PtrdiffT
	var dext noarch.PtrdiffT
	var lemax noarch.PtrdiffT
	var e noarch.PtrdiffT
	var elenk noarch.PtrdiffT
	var eln noarch.PtrdiffT
	var i noarch.PtrdiffT
	var j noarch.PtrdiffT
	var k noarch.PtrdiffT
	var k1 noarch.PtrdiffT
	var k2 noarch.PtrdiffT
	var k3 noarch.PtrdiffT
	var jlast noarch.PtrdiffT
	var ln noarch.PtrdiffT
	var dense noarch.PtrdiffT
	var nzmax noarch.PtrdiffT
	var mindeg noarch.PtrdiffT
	var nvi noarch.PtrdiffT
	var nvj noarch.PtrdiffT
	var nvk noarch.PtrdiffT
	var mark noarch.PtrdiffT
	var wnvi noarch.PtrdiffT
	var ok noarch.PtrdiffT
	var cnz noarch.PtrdiffT
	var nel noarch.PtrdiffT
	var p noarch.PtrdiffT
	var p1 noarch.PtrdiffT
	var p2 noarch.PtrdiffT
	var p3 noarch.PtrdiffT
	var p4 noarch.PtrdiffT
	var pj noarch.PtrdiffT
	var pk noarch.PtrdiffT
	var pk1 noarch.PtrdiffT
	var pk2 noarch.PtrdiffT
	var pn noarch.PtrdiffT
	var q noarch.PtrdiffT
	var n noarch.PtrdiffT
	var m noarch.PtrdiffT
	var t noarch.PtrdiffT
	var h noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || order <= 0 || order > noarch.PtrdiffT(3/8) {
		// --- Construct matrix C -----------------------------------------------
		// check
		return nil
	}
	// compute A'
	AT = cs_transpose(A, 0)
	if AT == nil {
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	// find dense threshold
	dense = noarch.PtrdiffT(func() float64 {
		if float64(16) > 10*math.Sqrt(float64(noarch.PtrdiffT(n))) {
			return float64((16))
		}
		return (10 * math.Sqrt(float64(noarch.PtrdiffT(n))))
	}())
	dense = noarch.PtrdiffT(func() int32 {
		if n-noarch.PtrdiffT(2/8) < dense {
			return (int32(n - noarch.PtrdiffT(2/8)))
		}
		return int32(noarch.PtrdiffT((dense)))
	}() / 8)
	if order == noarch.PtrdiffT(1/8) && n == m {
		// C = A+A'
		C = cs_add(A, AT, 0, 0)
	} else if order == noarch.PtrdiffT(2/8) {
		// drop dense columns from AT
		ATp = AT[0].p
		ATi = AT[0].i
		{
			p2 = 0
			j = 0
			for j = 0; j < m; j++ {
				// column j of AT starts here
				p = ATp[j]
				// new column j starts here
				ATp[j] = p2
				if ATp[j+noarch.PtrdiffT(1/8)]-p > dense {
					// skip dense col j
					continue
				}
				for ; p < ATp[j+noarch.PtrdiffT(1/8)]; p++ {
					ATi[func() noarch.PtrdiffT {
						defer func() {
							p2++
						}()
						return p2
					}()] = ATi[p]
				}
			}
		}
		// finalize AT
		ATp[m] = p2
		// A2 = AT'
		A2 = cs_transpose(AT, 0)
		// C=A'*A with no dense rows
		C = func() []cs {
			if A2 != nil {
				return cs_multiply(AT, A2)
			}
			return nil
		}()
		cs_spfree(A2)
	} else {
		// C=A'*A
		C = cs_multiply(AT, A)
	}
	cs_spfree(AT)
	if C == nil {
		return nil
	}
	// drop diagonal entries
	cs_fkeep(C, cs_diag, nil)
	Cp = C[0].p
	cnz = Cp[n]
	// allocate result
	P = cs_malloc(n+noarch.PtrdiffT(1/8), uint(0)).([]noarch.PtrdiffT)
	// get workspace
	W = cs_malloc(noarch.PtrdiffT(8*int32(n+noarch.PtrdiffT(1/8))/8), uint(0)).([]noarch.PtrdiffT)
	// add elbow room to C
	t = cnz + cnz/noarch.PtrdiffT(5/8) + noarch.PtrdiffT(2*int32(n)/8)
	if P == nil || W == nil || bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(C, noarch.PtrdiffT(t)))) {
		return (cs_idone(P, C, W, 0))
	}
	len = W
	nv = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(n+noarch.PtrdiffT(1/8)))*unsafe.Sizeof(W[0]))))[:]
	next = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(2*int32(n+noarch.PtrdiffT(1/8))))*unsafe.Sizeof(W[0]))))[:]
	head = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(3*int32(n+noarch.PtrdiffT(1/8))))*unsafe.Sizeof(W[0]))))[:]
	elen = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(4*int32(n+noarch.PtrdiffT(1/8))))*unsafe.Sizeof(W[0]))))[:]
	degree = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(5*int32(n+noarch.PtrdiffT(1/8))))*unsafe.Sizeof(W[0]))))[:]
	w = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(6*int32(n+noarch.PtrdiffT(1/8))))*unsafe.Sizeof(W[0]))))[:]
	hhead = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&W[0])) + (uintptr)(int(7*int32(n+noarch.PtrdiffT(1/8))))*unsafe.Sizeof(W[0]))))[:]
	// use P as workspace for last
	last = P
	{
		// --- Initialize quotient graph ----------------------------------------
		for k = 0; k < n; k++ {
			len[k] = Cp[k+noarch.PtrdiffT(1/8)] - Cp[k]
		}
	}
	len[n] = 0
	nzmax = noarch.PtrdiffT(C[0].nzmax)
	Ci = C[0].i
	for i = 0; i <= n; i++ {
		// degree list i is empty
		head[i] = noarch.PtrdiffT(-1)
		last[i] = noarch.PtrdiffT(-1)
		next[i] = noarch.PtrdiffT(-1)
		// hash list i is empty
		hhead[i] = noarch.PtrdiffT(-1)
		// node i is just one node
		nv[i] = 1
		// node i is alive
		w[i] = 1
		// Ek of node i is empty
		elen[i] = 0
		// degree of node i
		degree[i] = len[i]
	}
	// clear w
	mark = cs_wclear(0, 0, w, noarch.PtrdiffT(n))
	// n is a dead element
	elen[n] = noarch.PtrdiffT(-2)
	// n is a root of assembly tree
	Cp[n] = noarch.PtrdiffT(-1)
	// n is a dead element
	w[n] = 0
	{
		// --- Initialize degree lists ------------------------------------------
		for i = 0; i < n; i++ {
			d = degree[i]
			if d == noarch.PtrdiffT(0/8) {
				// node i is empty
				// element i is dead
				elen[i] = noarch.PtrdiffT(-2)
				nel++
				// i is a root of assembly tree
				Cp[i] = noarch.PtrdiffT(-1)
				w[i] = 0
			} else if d > dense {
				// node i is dense
				// absorb i into element n
				nv[i] = 0
				// node i is dead
				elen[i] = noarch.PtrdiffT(-1)
				nel++
				Cp[i] = -noarch.PtrdiffT((n)) - noarch.PtrdiffT(2/8)
				nv[n]++
			} else {
				if head[d] != noarch.PtrdiffT(int32(-1)/8) {
					last[head[d]] = i
				}
				// put node i in degree list d
				next[i] = head[d]
				head[d] = i
			}
		}
	}
	for nel < n {
		{
			// while (selecting pivots) do
			// --- Select node of minimum approximate degree --------------------
			for k = noarch.PtrdiffT(-1); mindeg < n && (func() noarch.PtrdiffT {
				k = head[mindeg]
				return k
			}()) == noarch.PtrdiffT(int32(-1)/8); mindeg++ {
			}
		}
		if next[k] != noarch.PtrdiffT(int32(-1)/8) {
			last[next[k]] = noarch.PtrdiffT(-1)
		}
		// remove k from degree list
		head[mindeg] = next[k]
		// elenk = |Ek|
		elenk = elen[k]
		// # of nodes k represents
		nvk = nv[k]
		// nv[k] nodes of A eliminated
		nel += nvk
		if elenk > noarch.PtrdiffT(0/8) && cnz+mindeg >= nzmax {
			{
				// --- Garbage collection -------------------------------------------
				for j = 0; j < n; j++ {
					if (func() noarch.PtrdiffT {
						p = Cp[j]
						return p
					}()) >= 0 {
						// j is a live node or element
						// save first entry of object
						Cp[j] = Ci[p]
						// first entry is now CS_FLIP(j)
						Ci[p] = -noarch.PtrdiffT((j)) - noarch.PtrdiffT(2/8)
					}
				}
			}
			{
				// scan all of memory
				q = 0
				p = 0
				for p = 0; p < cnz; {
					if (func() noarch.PtrdiffT {
						j = -noarch.PtrdiffT((Ci[func() noarch.PtrdiffT {
							defer func() {
								p++
							}()
							return p
						}()])) - noarch.PtrdiffT(2/8)
						return j
					}()) >= 0 {
						// found object j
						// restore first entry of object
						Ci[q] = Cp[j]
						// new pointer to object j
						Cp[j] = func() noarch.PtrdiffT {
							defer func() {
								q++
							}()
							return q
						}()
						for k3 = 0; k3 < len[j]-noarch.PtrdiffT(1/8); k3++ {
							Ci[func() noarch.PtrdiffT {
								defer func() {
									q++
								}()
								return q
							}()] = Ci[func() noarch.PtrdiffT {
								defer func() {
									p++
								}()
								return p
							}()]
						}
					}
				}
			}
			// Ci [cnz...nzmax-1] now free
			cnz = q
		}
		// --- Construct new element ----------------------------------------
		dk = 0
		// flag k as in Lk
		nv[k] = -noarch.PtrdiffT(nvk)
		p = Cp[k]
		// do in place if elen[k] == 0
		pk1 = noarch.PtrdiffT(func() int32 {
			if elenk == noarch.PtrdiffT(0/8) {
				return int32(noarch.PtrdiffT(p))
			}
			return int32(noarch.PtrdiffT(cnz))
		}() / 8)
		pk2 = pk1
		for k1 = 1; k1 <= elenk+noarch.PtrdiffT(1/8); k1++ {
			if k1 > elenk {
				// search the nodes in k
				e = k
				// list of nodes starts at Ci[pj]
				pj = p
				// length of list of nodes in k
				ln = len[k] - elenk
			} else {
				// search the nodes in e
				e = Ci[func() noarch.PtrdiffT {
					defer func() {
						p++
					}()
					return p
				}()]
				pj = Cp[e]
				// length of list of nodes in e
				ln = len[e]
			}
			for k2 = 1; k2 <= ln; k2++ {
				i = Ci[func() noarch.PtrdiffT {
					defer func() {
						pj++
					}()
					return pj
				}()]
				if (func() noarch.PtrdiffT {
					nvi = nv[i]
					return nvi
				}()) <= 0 {
					// node i dead, or seen
					continue
				}
				// degree[Lk] += size of node i
				dk += nvi
				// negate nv[i] to denote i in Lk
				nv[i] = -noarch.PtrdiffT(nvi)
				// place i in Lk
				Ci[func() noarch.PtrdiffT {
					defer func() {
						pk2++
					}()
					return pk2
				}()] = i
				if next[i] != noarch.PtrdiffT(int32(-1)/8) {
					last[next[i]] = last[i]
				}
				if last[i] != noarch.PtrdiffT(int32(-1)/8) {
					// remove i from degree list
					next[last[i]] = next[i]
				} else {
					head[degree[i]] = next[i]
				}
			}
			if e != k {
				// absorb e into k
				Cp[e] = -noarch.PtrdiffT((k)) - noarch.PtrdiffT(2/8)
				// e is now a dead element
				w[e] = 0
			}
		}
		if elenk != noarch.PtrdiffT(0/8) {
			// Ci [cnz...nzmax] is free
			cnz = pk2
		}
		// external degree of k - |Lk\i|
		degree[k] = dk
		// element k is in Ci[pk1..pk2-1]
		Cp[k] = pk1
		len[k] = pk2 - pk1
		// k is now an element
		elen[k] = noarch.PtrdiffT(-2)
		// --- Find set differences -----------------------------------------
		// clear w if necessary
		mark = cs_wclear(noarch.PtrdiffT(mark), noarch.PtrdiffT(lemax), w, noarch.PtrdiffT(n))
		{
			// scan 1: find |Le\Lk|
			for pk = pk1; pk < pk2; pk++ {
				i = Ci[pk]
				if (func() noarch.PtrdiffT {
					eln = elen[i]
					return eln
				}()) <= 0 {
					// skip if elen[i] empty
					continue
				}
				// nv [i] was negated
				nvi = -noarch.PtrdiffT(nv[i])
				wnvi = mark - nvi
				{
					// scan Ei
					for p = Cp[i]; p <= Cp[i]+eln-noarch.PtrdiffT(1/8); p++ {
						e = Ci[p]
						if w[e] >= mark {
							// decrement |Le\Lk|
							w[e] -= nvi
						} else if w[e] != noarch.PtrdiffT(0/8) {
							// ensure e is a live element
							// 1st time e seen in scan 1
							w[e] = degree[e] + wnvi
						}
					}
				}
			}
		}
		{
			// --- Degree update ------------------------------------------------
			// scan2: degree update
			for pk = pk1; pk < pk2; pk++ {
				// consider node i in Lk
				i = Ci[pk]
				p1 = Cp[i]
				p2 = p1 + elen[i] - noarch.PtrdiffT(1/8)
				pn = p1
				{
					// scan Ei
					h = 0
					d = 0
					p = p1
					for p = p1; p <= p2; p++ {
						e = Ci[p]
						if w[e] != noarch.PtrdiffT(0/8) {
							// e is an unabsorbed element
							// dext = |Le\Lk|
							dext = w[e] - mark
							if dext > noarch.PtrdiffT(0/8) {
								// sum up the set differences
								d += dext
								// keep e in Ei
								Ci[func() noarch.PtrdiffT {
									defer func() {
										pn++
									}()
									return pn
								}()] = e
								// compute the hash of node i
								h += e
							} else {
								// aggressive absorb. e->k
								Cp[e] = -noarch.PtrdiffT((k)) - noarch.PtrdiffT(2/8)
								// e is a dead element
								w[e] = 0
							}
						}
					}
				}
				// elen[i] = |Ei|
				elen[i] = pn - p1 + noarch.PtrdiffT(1/8)
				p3 = pn
				p4 = p1 + len[i]
				{
					// prune edges in Ai
					for p = p2 + noarch.PtrdiffT(1/8); p < p4; p++ {
						j = Ci[p]
						if (func() noarch.PtrdiffT {
							nvj = nv[j]
							return nvj
						}()) <= 0 {
							// node j dead or in Lk
							continue
						}
						// degree(i) += |j|
						d += nvj
						// place j in node list of i
						Ci[func() noarch.PtrdiffT {
							defer func() {
								pn++
							}()
							return pn
						}()] = j
						// compute hash for node i
						h += j
					}
				}
				if d == noarch.PtrdiffT(0/8) {
					// check for mass elimination
					// absorb i into k
					Cp[i] = -noarch.PtrdiffT((k)) - noarch.PtrdiffT(2/8)
					nvi = -noarch.PtrdiffT(nv[i])
					// |Lk| -= |i|
					dk -= nvi
					// |k| += nv[i]
					nvk += nvi
					nel += nvi
					nv[i] = 0
					// node i is dead
					elen[i] = noarch.PtrdiffT(-1)
				} else {
					// update degree(i)
					degree[i] = noarch.PtrdiffT(func() int32 {
						if degree[i] < d {
							return int32(noarch.PtrdiffT((degree[i])))
						}
						return int32(noarch.PtrdiffT((d)))
					}() / 8)
					// move first node to end
					Ci[pn] = Ci[p3]
					// move 1st el. to end of Ei
					Ci[p3] = Ci[p1]
					// add k as 1st element in of Ei
					Ci[p1] = k
					// new len of adj. list of node i
					len[i] = pn - p1 + noarch.PtrdiffT(1/8)
					// finalize hash of i
					h = noarch.PtrdiffT(func() int32 {
						if h < noarch.PtrdiffT(0/8) {
							return int32((-noarch.PtrdiffT(h)))
						}
						return int32(noarch.PtrdiffT(h))
					}() % int32(n) / 8)
					// place i in hash bucket
					next[i] = hhead[h]
					hhead[h] = i
					// save hash of i in last[i]
					last[i] = h
				}
			}
		}
		// scan2 is done
		// finalize |Lk|
		degree[k] = dk
		lemax = noarch.PtrdiffT(func() int32 {
			if lemax > dk {
				return int32(noarch.PtrdiffT((lemax)))
			}
			return int32(noarch.PtrdiffT((dk)))
		}() / 8)
		// clear w
		mark = cs_wclear(mark+lemax, noarch.PtrdiffT(lemax), w, noarch.PtrdiffT(n))
		{
			// --- Supernode detection ------------------------------------------
			for pk = pk1; pk < pk2; pk++ {
				i = Ci[pk]
				if nv[i] >= 0 {
					// skip if i is dead
					continue
				}
				// scan hash bucket of node i
				h = last[i]
				i = hhead[h]
				// hash bucket will be empty
				hhead[h] = noarch.PtrdiffT(-1)
				for i != noarch.PtrdiffT(int32(-1)/8) && next[i] != noarch.PtrdiffT(int32(-1)/8) {
					ln = len[i]
					eln = elen[i]
					for p = Cp[i] + noarch.PtrdiffT(1/8); p <= Cp[i]+ln-noarch.PtrdiffT(1/8); p++ {
						w[Ci[p]] = mark
					}
					jlast = i
					{
						// compare i with all j
						for j = next[i]; j != noarch.PtrdiffT(int32(-1)/8); {
							ok = noarch.PtrdiffT(len[j] == ln && elen[j] == eln)
							for p = Cp[j] + noarch.PtrdiffT(1/8); bool(ok) && p <= Cp[j]+ln-noarch.PtrdiffT(1/8); p++ {
								if w[Ci[p]] != mark {
									// compare i and j
									ok = 0
								}
							}
							if bool(noarch.PtrdiffT(ok)) {
								// i and j are identical
								// absorb j into i
								Cp[j] = -noarch.PtrdiffT((i)) - noarch.PtrdiffT(2/8)
								nv[i] += nv[j]
								nv[j] = 0
								// node j is dead
								elen[j] = noarch.PtrdiffT(-1)
								// delete j from hash bucket
								j = next[j]
								next[jlast] = j
							} else {
								// j and i are different
								jlast = j
								j = next[j]
							}
						}
					}
					i = next[i]
					mark++
				}
			}
		}
		{
			// --- Finalize new element------------------------------------------
			// finalize Lk
			p = pk1
			pk = pk1
			for pk = pk1; pk < pk2; pk++ {
				i = Ci[pk]
				if (func() noarch.PtrdiffT {
					nvi = -noarch.PtrdiffT(nv[i])
					return nvi
				}()) <= 0 {
					// skip if i is dead
					continue
				}
				// restore nv[i]
				nv[i] = nvi
				// compute external degree(i)
				d = degree[i] + dk - nvi
				d = noarch.PtrdiffT(func() int32 {
					if d < n-nel-nvi {
						return int32(noarch.PtrdiffT((d)))
					}
					return (int32(n - nel - nvi))
				}() / 8)
				if head[d] != noarch.PtrdiffT(int32(-1)/8) {
					last[head[d]] = i
				}
				// put i back in degree list
				next[i] = head[d]
				last[i] = noarch.PtrdiffT(-1)
				head[d] = i
				// find new minimum degree
				mindeg = noarch.PtrdiffT(func() int32 {
					if mindeg < d {
						return int32(noarch.PtrdiffT((mindeg)))
					}
					return int32(noarch.PtrdiffT((d)))
				}() / 8)
				degree[i] = d
				// place i in Lk
				Ci[func() noarch.PtrdiffT {
					defer func() {
						p++
					}()
					return p
				}()] = i
			}
		}
		// # nodes absorbed into k
		nv[k] = nvk
		if (func() noarch.PtrdiffT {
			len[k] = p - pk1
			return len[k]
		}()) == noarch.PtrdiffT(0/8) {
			// length of adj list of element k
			// k is a root of the tree
			Cp[k] = noarch.PtrdiffT(-1)
			// k is now a dead element
			w[k] = 0
		}
		if elenk != noarch.PtrdiffT(0/8) {
			// free unused space in Lk
			cnz = p
		}
	}
	{
		// --- Postordering -----------------------------------------------------
		// fix assembly tree
		for i = 0; i < n; i++ {
			Cp[i] = -noarch.PtrdiffT((Cp[i])) - noarch.PtrdiffT(2/8)
		}
	}
	for j = 0; j <= n; j++ {
		head[j] = noarch.PtrdiffT(-1)
	}
	{
		// place unordered nodes in lists
		for j = n; j >= 0; j-- {
			if nv[j] > noarch.PtrdiffT(0/8) {
				// skip if j is an element
				continue
			}
			// place j in list of its parent
			next[j] = head[Cp[j]]
			head[Cp[j]] = j
		}
	}
	{
		// place elements in lists
		for e = n; e >= 0; e-- {
			if nv[e] <= 0 {
				// skip unless e is an element
				continue
			}
			if Cp[e] != noarch.PtrdiffT(int32(-1)/8) {
				// place e in list of its parent
				next[e] = head[Cp[e]]
				head[Cp[e]] = e
			}
		}
	}
	{
		// postorder the assembly tree
		k = 0
		i = 0
		for i = 0; i <= n; i++ {
			if Cp[i] == noarch.PtrdiffT(int32(-1)/8) {
				k = cs_tdfs(noarch.PtrdiffT(i), noarch.PtrdiffT(k), head, next, P, w)
			}
		}
	}
	return (cs_idone(P, C, W, 1))
}

// cs_chol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_chol.c:3
// L = chol (A, [pinv parent cp]), pinv is optional
func cs_chol(A []cs, S []css) []csn {
	var d float64
	var lki float64
	var Lx []float64
	var x []float64
	var Cx []float64
	var top noarch.PtrdiffT
	var i noarch.PtrdiffT
	var p noarch.PtrdiffT
	var k noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Li []noarch.PtrdiffT
	var Lp []noarch.PtrdiffT
	var cp []noarch.PtrdiffT
	var pinv []noarch.PtrdiffT
	var s []noarch.PtrdiffT
	var c []noarch.PtrdiffT
	var parent []noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var L []cs
	var C []cs
	var E []cs
	var N []csn
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || S == nil || S[0].cp == nil || S[0].parent == nil {
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	// allocate result
	N = cs_calloc(1, uint(32)).([]csn)
	// get csi workspace
	c = cs_malloc(noarch.PtrdiffT(2*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
	// get double workspace
	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
	cp = S[0].cp
	pinv = S[0].pinv
	parent = S[0].parent
	C = func() []cs {
		if pinv != nil {
			return cs_symperm(A, pinv, 1)
		}
		return (A)
	}()
	// E is alias for A, or a copy E=A(p,p)
	E = func() []cs {
		if pinv != nil {
			return C
		}
		return nil
	}()
	if N == nil || c == nil || x == nil || C == nil {
		return (cs_ndone(N, E, c, x, 0))
	}
	s = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&c[0])) + (uintptr)(int(n))*unsafe.Sizeof(c[0]))))[:]
	Cp = C[0].p
	Ci = C[0].i
	Cx = C[0].x
	L = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(cp[n]), 1, 0)
	// allocate result
	N[0].L = L
	if L == nil {
		return (cs_ndone(N, E, c, x, 0))
	}
	Lp = L[0].p
	Li = L[0].i
	Lx = L[0].x
	for k = 0; k < n; k++ {
		c[k] = cp[k]
		Lp[k] = c[k]
	}
	{
		// compute L(k,:) for L*L' = C
		for k = 0; k < n; k++ {
			// --- Nonzero pattern of L(k,:) ------------------------------------
			// find pattern of L(k,:)
			top = cs_ereach(C, noarch.PtrdiffT(k), parent, s, c)
			// x (0:k) is now zero
			x[k] = 0
			{
				// x = full(triu(C(:,k)))
				for p = Cp[k]; p < Cp[k+noarch.PtrdiffT(1/8)]; p++ {
					if Ci[p] <= k {
						x[Ci[p]] = Cx[p]
					}
				}
			}
			// d = C(k,k)
			d = x[k]
			// clear x for k+1st iteration
			x[k] = 0
			for ; top < n; top++ {
				// --- Triangular solve ---------------------------------------------
				// solve L(0:k-1,0:k-1) * x = C(:,k)
				// s [top..n-1] is pattern of L(k,:)
				i = s[top]
				// L(k,i) = x (i) / L(i,i)
				lki = x[i] / Lx[Lp[i]]
				// clear x for k+1st iteration
				x[i] = 0
				for p = Lp[i] + noarch.PtrdiffT(1/8); p < c[i]; p++ {
					x[Li[p]] -= Lx[p] * lki
				}
				// d = d - L(k,i)*L(k,i)
				d -= lki * lki
				p = func() noarch.PtrdiffT {
					tempVar := &c[i]
					defer func() {
						*tempVar++
					}()
					return *tempVar
				}()
				// store L(k,i) in column i
				Li[p] = k
				Lx[p] = lki
			}
			if d <= 0 {
				// --- Compute L(k,k) -----------------------------------------------
				// not pos def
				return (cs_ndone(N, E, c, x, 0))
			}
			p = func() noarch.PtrdiffT {
				tempVar := &c[k]
				defer func() {
					*tempVar++
				}()
				return *tempVar
			}()
			// store L(k,k) = sqrt (d) in column k
			Li[p] = k
			Lx[p] = math.Sqrt(d)
		}
	}
	// finalize L
	Lp[n] = cp[n]
	// success: free E,s,x; return N
	return (cs_ndone(N, E, c, x, 1))
}

// cs_cholsol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_cholsol.c:3
// x=A\b where A is symmetric positive definite; b overwritten with solution
func cs_cholsol(order noarch.PtrdiffT, A []cs, b []float64) noarch.PtrdiffT {
	var x []float64
	var S []css
	var N []csn
	var n noarch.PtrdiffT
	var ok noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || b == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(A[0].n)
	// ordering and symbolic analysis
	S = cs_schol(noarch.PtrdiffT(order), A)
	// numeric Cholesky factorization
	N = cs_chol(A, S)
	// get workspace
	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
	ok = noarch.PtrdiffT(S != nil && N != nil && x != nil)
	if bool(noarch.PtrdiffT(ok)) {
		// x = P*b
		cs_ipvec(S[0].pinv, b, x, noarch.PtrdiffT(n))
		// x = L\x
		cs_lsolve(N[0].L, x)
		// x = L'\x
		cs_ltsolve(N[0].L, x)
		// b = P'*x
		cs_pvec(S[0].pinv, x, b, noarch.PtrdiffT(n))
	}
	cs_free(x)
	cs_sfree(S)
	cs_nfree(N)
	return noarch.PtrdiffT((ok))
}

// cs_compress - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_compress.c:3
// C = compressed-column form of a triplet matrix T
func cs_compress(T []cs) []cs {
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var p noarch.PtrdiffT
	var k noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var Ti []noarch.PtrdiffT
	var Tj []noarch.PtrdiffT
	var Cx []float64
	var Tx []float64
	var C []cs
	if !(T != nil && noarch.PtrdiffT(T[0].nz) >= 0) {
		// check inputs
		return nil
	}
	m = noarch.PtrdiffT(T[0].m)
	n = noarch.PtrdiffT(T[0].n)
	Ti = T[0].i
	Tj = T[0].p
	Tx = T[0].x
	nz = noarch.PtrdiffT(T[0].nz)
	// allocate result
	C = cs_spalloc(noarch.PtrdiffT(m), noarch.PtrdiffT(n), noarch.PtrdiffT(nz), noarch.PtrdiffT(Tx != nil), 0)
	// get workspace
	w = cs_calloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	if C == nil || w == nil {
		// out of memory
		return (cs_done(C, w, nil, 0))
	}
	Cp = C[0].p
	Ci = C[0].i
	Cx = C[0].x
	{
		// column counts
		for k = 0; k < nz; k++ {
			w[Tj[k]]++
		}
	}
	// column pointers
	cs_cumsum(Cp, w, noarch.PtrdiffT(n))
	for k = 0; k < nz; k++ {
		// A(i,j) is the pth entry in C
		Ci[(func() noarch.PtrdiffT {
			p = func() noarch.PtrdiffT {
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
	return (cs_done(C, w, nil, 1))
}

// init_ata - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_counts.c:5
// column counts of LL'=A or LL'=A'A, given parent & post ordering
func init_ata(AT []cs, post []noarch.PtrdiffT, w []noarch.PtrdiffT, head [][]noarch.PtrdiffT, next [][]noarch.PtrdiffT) {
	var i noarch.PtrdiffT
	var k noarch.PtrdiffT
	var p noarch.PtrdiffT
	var m noarch.PtrdiffT = noarch.PtrdiffT(AT[0].n)
	var n noarch.PtrdiffT = noarch.PtrdiffT(AT[0].m)
	var ATp []noarch.PtrdiffT = AT[0].p
	var ATi []noarch.PtrdiffT = AT[0].i
	head[0] = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(4*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	next[0] = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(5*int32(n)))*unsafe.Sizeof(w[0]))))[:][1:]
	{
		// invert post
		for k = 0; k < n; k++ {
			w[post[k]] = k
		}
	}
	for i = 0; i < m; i++ {
		{
			k = n
			p = ATp[i]
			for p = ATp[i]; p < ATp[i+noarch.PtrdiffT(1/8)]; p++ {
				k = noarch.PtrdiffT(func() int32 {
					if k < w[ATi[p]] {
						return int32(noarch.PtrdiffT((k)))
					}
					return int32(noarch.PtrdiffT((w[ATi[p]])))
				}() / 8)
			}
		}
		// place row i in linked list k
		(next[0])[i] = (head[0])[k]
		(head[0])[k] = i
	}
}

// cs_counts - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_counts.c:17
func cs_counts(A []cs, parent []noarch.PtrdiffT, post []noarch.PtrdiffT, ata noarch.PtrdiffT) []noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var j noarch.PtrdiffT
	var k noarch.PtrdiffT
	var n noarch.PtrdiffT
	var m noarch.PtrdiffT
	var J noarch.PtrdiffT
	var s noarch.PtrdiffT
	var p noarch.PtrdiffT
	var q noarch.PtrdiffT
	var jleaf noarch.PtrdiffT
	var ATp []noarch.PtrdiffT
	var ATi []noarch.PtrdiffT
	var maxfirst []noarch.PtrdiffT
	var prevleaf []noarch.PtrdiffT
	var ancestor []noarch.PtrdiffT
	var head []noarch.PtrdiffT
	var next []noarch.PtrdiffT
	var colcount []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var first []noarch.PtrdiffT
	var delta []noarch.PtrdiffT
	var AT []cs
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || parent == nil || post == nil {
		// check inputs
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	s = noarch.PtrdiffT((4*int32(n) + func() int32 {
		if bool(noarch.PtrdiffT(ata)) {
			return (int32(n + m + noarch.PtrdiffT(1/8)))
		}
		return 0
	}()) / 8)
	colcount = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	// allocate result
	delta = colcount
	// get workspace
	w = cs_malloc(noarch.PtrdiffT(s), uint(0)).([]noarch.PtrdiffT)
	// AT = A'
	AT = cs_transpose(A, 0)
	if AT == nil || colcount == nil || w == nil {
		return (cs_idone(colcount, AT, w, 0))
	}
	ancestor = w
	maxfirst = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
	prevleaf = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	first = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(3*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	{
		// clear workspace w [0..s-1]
		for k = 0; k < s; k++ {
			w[k] = noarch.PtrdiffT(-1)
		}
	}
	{
		// find first [j]
		for k = 0; k < n; k++ {
			j = post[k]
			// delta[j]=1 if j is a leaf
			delta[j] = noarch.PtrdiffT(func() int {
				if first[j] == noarch.PtrdiffT(int32(-1)/8) {
					return 1
				}
				return 0
			}())
			for ; j != noarch.PtrdiffT(int32(-1)/8) && first[j] == noarch.PtrdiffT(int32(-1)/8); j = parent[j] {
				first[j] = k
			}
		}
	}
	ATp = AT[0].p
	ATi = AT[0].i
	if bool(noarch.PtrdiffT(ata)) {
		init_ata(AT, post, w, (*[100000000][]noarch.PtrdiffT)(unsafe.Pointer(&head))[:], (*[100000000][]noarch.PtrdiffT)(unsafe.Pointer(&next))[:])
	}
	{
		// each node in its own set
		for i = 0; i < n; i++ {
			ancestor[i] = i
		}
	}
	for k = 0; k < n; k++ {
		// j is the kth node in postordered etree
		j = post[k]
		if parent[j] != noarch.PtrdiffT(int32(-1)/8) {
			// j is not a root
			delta[parent[j]]--
		}
		{
			// J=j for LL'=A case
			for J = noarch.PtrdiffT(func() int32 {
				if bool(noarch.PtrdiffT(ata)) {
					return int32(noarch.PtrdiffT(head[k]))
				}
				return int32(noarch.PtrdiffT(j))
			}() / 8); J != noarch.PtrdiffT(int32(-1)/8); J = noarch.PtrdiffT(func() int32 {
				if bool(noarch.PtrdiffT(ata)) {
					return int32(noarch.PtrdiffT(next[J]))
				}
				return int32(-1)
			}() / 8) {
				for p = ATp[J]; p < ATp[J+noarch.PtrdiffT(1/8)]; p++ {
					i = ATi[p]
					q = cs_leaf(noarch.PtrdiffT(i), noarch.PtrdiffT(j), first, maxfirst, prevleaf, ancestor, (*[100000000]noarch.PtrdiffT)(unsafe.Pointer(&jleaf))[:])
					if jleaf >= 1 {
						// A(i,j) is in skeleton
						delta[j]++
					}
					if jleaf == noarch.PtrdiffT(2/8) {
						// account for overlap in q
						delta[q]--
					}
				}
			}
		}
		if parent[j] != noarch.PtrdiffT(int32(-1)/8) {
			ancestor[j] = parent[j]
		}
	}
	{
		// sum up delta's of each child
		for j = 0; j < n; j++ {
			if parent[j] != noarch.PtrdiffT(int32(-1)/8) {
				colcount[parent[j]] += colcount[j]
			}
		}
	}
	// success: free workspace
	return (cs_idone(colcount, AT, w, 1))
}

// cs_cumsum - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_cumsum.c:3
// p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
func cs_cumsum(p []noarch.PtrdiffT, c []noarch.PtrdiffT, n noarch.PtrdiffT) float64 {
	var i noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var nz2 float64
	if p == nil || c == nil {
		// check inputs
		return float64((-1))
	}
	for i = 0; i < n; i++ {
		p[i] = nz
		nz += c[i]
		// also in double to avoid csi overflow
		nz2 += float64(c[i])
		// also copy p[0..n-1] back into c[0..n-1]
		c[i] = p[i]
	}
	p[n] = nz
	// return sum (c [0..n-1])
	return (nz2)
}

// cs_dfs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dfs.c:3
// depth-first-search of the graph of a matrix, starting at node j
func cs_dfs(j noarch.PtrdiffT, G []cs, top noarch.PtrdiffT, xi []noarch.PtrdiffT, pstack []noarch.PtrdiffT, pinv []noarch.PtrdiffT) noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var p noarch.PtrdiffT
	var p2 noarch.PtrdiffT
	var done noarch.PtrdiffT
	var jnew noarch.PtrdiffT
	var head noarch.PtrdiffT
	var Gp []noarch.PtrdiffT
	var Gi []noarch.PtrdiffT
	if !(G != nil && noarch.PtrdiffT(G[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || xi == nil || pstack == nil {
		// check inputs
		return noarch.PtrdiffT((-1))
	}
	Gp = G[0].p
	Gi = G[0].i
	// initialize the recursion stack
	xi[0] = j
	for head >= 0 {
		// get j from the top of the recursion stack
		j = xi[head]
		jnew = noarch.PtrdiffT(func() int32 {
			if pinv != nil {
				return int32(noarch.PtrdiffT((pinv[j])))
			}
			return int32(noarch.PtrdiffT(j))
		}() / 8)
		if !(Gp[j] < noarch.PtrdiffT(0/8)) {
			{
				// mark node j as visited
				Gp[j] = -noarch.PtrdiffT((Gp[j])) - noarch.PtrdiffT(2/8)
			}
			pstack[head] = noarch.PtrdiffT(func() int32 {
				if jnew < noarch.PtrdiffT(0/8) {
					return 0
				}
				return (func() int32 {
					if Gp[jnew] < noarch.PtrdiffT(0/8) {
						return (int32(-noarch.PtrdiffT((Gp[jnew])) - noarch.PtrdiffT(2/8)))
					}
					return int32(noarch.PtrdiffT((Gp[jnew])))
				}())
			}() / 8)
		}
		// node j done if no unvisited neighbors
		done = 1
		p2 = noarch.PtrdiffT(func() int32 {
			if jnew < noarch.PtrdiffT(0/8) {
				return 0
			}
			return (func() int32 {
				if Gp[jnew+noarch.PtrdiffT(1/8)] < noarch.PtrdiffT(0/8) {
					return (int32(-noarch.PtrdiffT((Gp[jnew+noarch.PtrdiffT(1/8)])) - noarch.PtrdiffT(2/8)))
				}
				return int32(noarch.PtrdiffT((Gp[jnew+noarch.PtrdiffT(1/8)])))
			}())
		}() / 8)
		{
			// examine all neighbors of j
			for p = pstack[head]; p < p2; p++ {
				// consider neighbor node i
				i = Gi[p]
				if Gp[i] < noarch.PtrdiffT(0/8) {
					// skip visited node i
					continue
				}
				// pause depth-first search of node j
				pstack[head] = p
				// start dfs at node i
				xi[func() noarch.PtrdiffT {
					head++
					return head
				}()] = i
				// node j is not done
				done = 0
				// break, to start dfs (i)
				break
			}
		}
		if bool(noarch.PtrdiffT(done)) {
			// depth-first search at node j is done
			// remove j from the recursion stack
			head--
			// and place in the output stack
			xi[func() noarch.PtrdiffT {
				top--
				return top
			}()] = j
		}
	}
	return noarch.PtrdiffT((top))
}

// cs_bfs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:3
// breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3)
func cs_bfs(A []cs, n noarch.PtrdiffT, wi []noarch.PtrdiffT, wj []noarch.PtrdiffT, queue []noarch.PtrdiffT, imatch []noarch.PtrdiffT, jmatch []noarch.PtrdiffT, mark noarch.PtrdiffT) noarch.PtrdiffT {
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var head noarch.PtrdiffT
	var tail noarch.PtrdiffT
	var j noarch.PtrdiffT
	var i noarch.PtrdiffT
	var p noarch.PtrdiffT
	var j2 noarch.PtrdiffT
	var C []cs
	{
		// place all unmatched nodes in queue
		for j = 0; j < n; j++ {
			if imatch[j] >= 0 {
				// skip j if matched
				continue
			}
			// j in set C0 (R0 if transpose)
			wj[j] = 0
			// place unmatched col j in queue
			queue[func() noarch.PtrdiffT {
				defer func() {
					tail++
				}()
				return tail
			}()] = j
		}
	}
	if tail == noarch.PtrdiffT(0/8) {
		// quick return if no unmatched nodes
		return noarch.PtrdiffT((1))
	}
	C = func() []cs {
		if mark == noarch.PtrdiffT(1/8) {
			return (A)
		}
		return cs_transpose(A, 0)
	}()
	if C == nil {
		// bfs of C=A' to find R3,C3 from R0
		return noarch.PtrdiffT((0))
	}
	Ap = C[0].p
	Ai = C[0].i
	for head < tail {
		// while queue is not empty
		// get the head of the queue
		j = queue[func() noarch.PtrdiffT {
			defer func() {
				head++
			}()
			return head
		}()]
		for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			i = Ai[p]
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
			queue[func() noarch.PtrdiffT {
				defer func() {
					tail++
				}()
				return tail
			}()] = j2
		}
	}
	if mark != noarch.PtrdiffT(1/8) {
		// free A' if it was created
		cs_spfree(C)
	}
	return noarch.PtrdiffT((1))
}

// cs_matched - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:37
// collect matched rows and columns into p and q
func cs_matched(n noarch.PtrdiffT, wj []noarch.PtrdiffT, imatch []noarch.PtrdiffT, p []noarch.PtrdiffT, q []noarch.PtrdiffT, cc []noarch.PtrdiffT, rr []noarch.PtrdiffT, set noarch.PtrdiffT, mark noarch.PtrdiffT) {
	var kc noarch.PtrdiffT = cc[set]
	var j noarch.PtrdiffT
	var kr noarch.PtrdiffT = rr[set-noarch.PtrdiffT(1/8)]
	for j = 0; j < n; j++ {
		if wj[j] != mark {
			// skip if j is not in C set
			continue
		}
		p[func() noarch.PtrdiffT {
			defer func() {
				kr++
			}()
			return kr
		}()] = imatch[j]
		q[func() noarch.PtrdiffT {
			defer func() {
				kc++
			}()
			return kc
		}()] = j
	}
	cc[set+noarch.PtrdiffT(1/8)] = kc
	rr[set] = kr
}

// cs_unmatched - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:53
// collect unmatched rows into the permutation vector p
func cs_unmatched(m noarch.PtrdiffT, wi []noarch.PtrdiffT, p []noarch.PtrdiffT, rr []noarch.PtrdiffT, set noarch.PtrdiffT) {
	var i noarch.PtrdiffT
	var kr noarch.PtrdiffT = rr[set]
	for i = 0; i < m; i++ {
		if wi[i] == noarch.PtrdiffT(0/8) {
			p[func() noarch.PtrdiffT {
				defer func() {
					kr++
				}()
				return kr
			}()] = i
		}
	}
	rr[set+noarch.PtrdiffT(1/8)] = kr
}

// cs_rprune - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:61
// return 1 if row i is in R2
func cs_rprune(i noarch.PtrdiffT, j noarch.PtrdiffT, aij float64, other interface{}) noarch.PtrdiffT {
	var rr []noarch.PtrdiffT = other.([]noarch.PtrdiffT)
	return noarch.PtrdiffT((i >= rr[1] && i < rr[2]))
}

// cs_dmperm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dmperm.c:68
// Given A, compute coarse and then fine dmperm
func cs_dmperm(A []cs, seed noarch.PtrdiffT) []csd {
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var i noarch.PtrdiffT
	var j noarch.PtrdiffT
	var k noarch.PtrdiffT
	var cnz noarch.PtrdiffT
	var nc noarch.PtrdiffT
	var jmatch []noarch.PtrdiffT
	var imatch []noarch.PtrdiffT
	var wi []noarch.PtrdiffT
	var wj []noarch.PtrdiffT
	var pinv []noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var ps []noarch.PtrdiffT
	var rs []noarch.PtrdiffT
	var nb1 noarch.PtrdiffT
	var nb2 noarch.PtrdiffT
	var p []noarch.PtrdiffT
	var q []noarch.PtrdiffT
	var cc []noarch.PtrdiffT
	var rr []noarch.PtrdiffT
	var r []noarch.PtrdiffT
	var s []noarch.PtrdiffT
	var ok noarch.PtrdiffT
	var C []cs
	var D []csd
	var scc []csd
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// --- Maximum matching -------------------------------------------------
		// check inputs
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	// allocate result
	D = cs_dalloc(noarch.PtrdiffT(m), noarch.PtrdiffT(n))
	if D == nil {
		return nil
	}
	p = D[0].p
	q = D[0].q
	r = D[0].r
	s = D[0].s
	cc = D[0].cc[:]
	rr = D[0].rr[:]
	// max transversal
	jmatch = cs_maxtrans(A, noarch.PtrdiffT(seed))
	// imatch = inverse of jmatch
	imatch = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&jmatch[0])) + (uintptr)(int(m))*unsafe.Sizeof(jmatch[0]))))[:]
	if jmatch == nil {
		return (cs_ddone(D, nil, jmatch, 0))
	}
	// --- Coarse decomposition ---------------------------------------------
	// use r and s as workspace
	wi = r
	wj = s
	{
		// unmark all cols for bfs
		for j = 0; j < n; j++ {
			wj[j] = noarch.PtrdiffT(-1)
		}
	}
	{
		// unmark all rows for bfs
		for i = 0; i < m; i++ {
			wi[i] = noarch.PtrdiffT(-1)
		}
	}
	// find C1, R1 from C0
	cs_bfs(A, noarch.PtrdiffT(n), wi, wj, q, imatch, jmatch, 1)
	// find R3, C3 from R0
	ok = cs_bfs(A, noarch.PtrdiffT(m), wj, wi, p, jmatch, imatch, 3)
	if bool(noarch.NotNoarch.PtrdiffT(noarch.PtrdiffT(ok))) {
		return (cs_ddone(D, nil, jmatch, 0))
	}
	// unmatched set C0
	cs_unmatched(noarch.PtrdiffT(n), wj, q, cc, 0)
	// set R1 and C1
	cs_matched(noarch.PtrdiffT(n), wj, imatch, p, q, cc, rr, 1, 1)
	// set R2 and C2
	cs_matched(noarch.PtrdiffT(n), wj, imatch, p, q, cc, rr, 2, noarch.PtrdiffT(-1))
	// set R3 and C3
	cs_matched(noarch.PtrdiffT(n), wj, imatch, p, q, cc, rr, 3, 3)
	// unmatched set R0
	cs_unmatched(noarch.PtrdiffT(m), wi, p, rr, 3)
	cs_free(jmatch)
	// --- Fine decomposition -----------------------------------------------
	// pinv=p'
	pinv = cs_pinv(p, noarch.PtrdiffT(m))
	if pinv == nil {
		return (cs_ddone(D, nil, nil, 0))
	}
	// C=A(p,q) (it will hold A(R2,C2))
	C = cs_permute(A, pinv, q, 0)
	cs_free(pinv)
	if C == nil {
		return (cs_ddone(D, nil, nil, 0))
	}
	Cp = C[0].p
	// delete cols C0, C1, and C3 from C
	nc = cc[3] - cc[2]
	if cc[2] > noarch.PtrdiffT(0/8) {
		for j = cc[2]; j <= cc[3]; j++ {
			Cp[j-cc[2]] = Cp[j]
		}
	}
	C[0].n = nc
	if rr[2]-rr[1] < m {
		// delete rows R0, R1, and R3 from C
		cs_fkeep(C, cs_rprune, rr)
		cnz = Cp[nc]
		Ci = C[0].i
		if rr[1] > noarch.PtrdiffT(0/8) {
			for k = 0; k < cnz; k++ {
				Ci[k] -= rr[1]
			}
		}
	}
	C[0].m = nc
	// find strongly connected components of C
	scc = cs_scc(C)
	if scc == nil {
		return (cs_ddone(D, C, nil, 0))
	}
	// --- Combine coarse and fine decompositions ---------------------------
	// C(ps,ps) is the permuted matrix
	ps = scc[0].p
	// kth block is rs[k]..rs[k+1]-1
	rs = scc[0].r
	// # of blocks of A(R2,C2)
	nb1 = noarch.PtrdiffT(scc[0].nb)
	for k = 0; k < nc; k++ {
		wj[k] = q[ps[k]+cc[2]]
	}
	for k = 0; k < nc; k++ {
		q[k+cc[2]] = wj[k]
	}
	for k = 0; k < nc; k++ {
		wi[k] = p[ps[k]+rr[1]]
	}
	for k = 0; k < nc; k++ {
		p[k+rr[1]] = wi[k]
	}
	// create the fine block partitions
	nb2 = 0
	s[0] = 0
	r[0] = s[0]
	if cc[2] > noarch.PtrdiffT(0/8) {
		// leading coarse block A (R1, [C0 C1])
		nb2++
	}
	{
		// coarse block A (R2,C2)
		for k = 0; k < nb1; k++ {
			// A (R2,C2) splits into nb1 fine blocks
			r[nb2] = rs[k] + rr[1]
			s[nb2] = rs[k] + cc[2]
			nb2++
		}
	}
	if rr[2] < m {
		// trailing coarse block A ([R3 R0], C3)
		r[nb2] = rr[2]
		s[nb2] = cc[3]
		nb2++
	}
	r[nb2] = m
	s[nb2] = n
	D[0].nb = nb2
	cs_dfree(scc)
	return (cs_ddone(D, C, nil, 1))
}

// cs_tol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_droptol.c:2
func cs_tol(i noarch.PtrdiffT, j noarch.PtrdiffT, aij float64, tol interface{}) noarch.PtrdiffT {
	return noarch.PtrdiffT((math.Abs(aij) > (tol.([]float64))[(0)]))
}

// cs_droptol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_droptol.c:6
func cs_droptol(A []cs, tol float64) noarch.PtrdiffT {
	// keep all large entries
	return (cs_fkeep(A, cs_tol, (*[100000000]float64)(unsafe.Pointer(&tol))[:]))
}

// cs_nonzero - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dropzeros.c:2
func cs_nonzero(i noarch.PtrdiffT, j noarch.PtrdiffT, aij float64, other interface{}) noarch.PtrdiffT {
	return noarch.PtrdiffT((aij != 0))
}

// cs_dropzeros - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dropzeros.c:6
func cs_dropzeros(A []cs) noarch.PtrdiffT {
	// keep all nonzero entries
	return (cs_fkeep(A, cs_nonzero, nil))
}

// cs_dupl - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_dupl.c:3
// remove duplicate entries from A
func cs_dupl(A []cs) noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var j noarch.PtrdiffT
	var p noarch.PtrdiffT
	var q noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var n noarch.PtrdiffT
	var m noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var Ax []float64
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	// get workspace
	w = cs_malloc(noarch.PtrdiffT(m), uint(0)).([]noarch.PtrdiffT)
	if w == nil {
		// out of memory
		return noarch.PtrdiffT((0))
	}
	{
		// row i not yet seen
		for i = 0; i < m; i++ {
			w[i] = noarch.PtrdiffT(-1)
		}
	}
	for j = 0; j < n; j++ {
		// column j will start at q
		q = nz
		for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			// A(i,j) is nonzero
			i = Ai[p]
			if w[i] >= q {
				// A(i,j) is a duplicate
				Ax[w[i]] += Ax[p]
			} else {
				// record where row i occurs
				w[i] = nz
				// keep A(i,j)
				Ai[nz] = i
				Ax[func() noarch.PtrdiffT {
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

// cs_entry - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_entry.c:3
// add an entry to a triplet matrix; return 1 if ok, 0 otherwise
func cs_entry(T []cs, i noarch.PtrdiffT, j noarch.PtrdiffT, x float64) noarch.PtrdiffT {
	if !(T != nil && noarch.PtrdiffT(T[0].nz) >= 0) || i < noarch.PtrdiffT(0/8) || j < noarch.PtrdiffT(0/8) {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	if noarch.PtrdiffT(T[0].nz) >= noarch.PtrdiffT(T[0].nzmax) && bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(T, noarch.PtrdiffT(2*int32(T[0].nzmax)/8)))) {
		return noarch.PtrdiffT((0))
	}
	if T[0].x != nil {
		T[0].x[noarch.PtrdiffT(T[0].nz)] = x
	}
	T[0].i[noarch.PtrdiffT(T[0].nz)] = i
	T[0].p[func() noarch.PtrdiffT {
		tempVar := &T[0].nz
		defer func() {
			*tempVar++
		}()
		return *tempVar
	}()] = j
	T[0].m = noarch.PtrdiffT(func() int32 {
		if T[0].m > i+noarch.PtrdiffT(1/8) {
			return int32(noarch.PtrdiffT((T[0].m)))
		}
		return (int32(i + noarch.PtrdiffT(1/8)))
	}() / 8)
	T[0].n = noarch.PtrdiffT(func() int32 {
		if T[0].n > j+noarch.PtrdiffT(1/8) {
			return int32(noarch.PtrdiffT((T[0].n)))
		}
		return (int32(j + noarch.PtrdiffT(1/8)))
	}() / 8)
	return noarch.PtrdiffT((1))
}

// cs_ereach - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_ereach.c:3
// find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k))
func cs_ereach(A []cs, k noarch.PtrdiffT, parent []noarch.PtrdiffT, s []noarch.PtrdiffT, w []noarch.PtrdiffT) noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var p noarch.PtrdiffT
	var n noarch.PtrdiffT
	var len noarch.PtrdiffT
	var top noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || parent == nil || s == nil || w == nil {
		// check inputs
		return noarch.PtrdiffT((-1))
	}
	n = noarch.PtrdiffT(A[0].n)
	top = n
	Ap = A[0].p
	Ai = A[0].i
	{
		// mark node k as visited
		w[k] = -noarch.PtrdiffT((w[k])) - noarch.PtrdiffT(2/8)
	}
	for p = Ap[k]; p < Ap[k+noarch.PtrdiffT(1/8)]; p++ {
		// A(i,k) is nonzero
		i = Ai[p]
		if i > k {
			// only use upper triangular part of A
			continue
		}
		{
			// traverse up etree
			for len = 0; !(w[i] < noarch.PtrdiffT(0/8)); i = parent[i] {
				// L(k,i) is nonzero
				s[func() noarch.PtrdiffT {
					defer func() {
						len++
					}()
					return len
				}()] = i
				{
					// mark i as visited
					w[i] = -noarch.PtrdiffT((w[i])) - noarch.PtrdiffT(2/8)
				}
			}
		}
		for len > noarch.PtrdiffT(0/8) {
			// push path onto stack
			s[func() noarch.PtrdiffT {
				top--
				return top
			}()] = s[func() noarch.PtrdiffT {
				len--
				return len
			}()]
		}
	}
	{
		// unmark all nodes
		for p = top; p < n; p++ {
			w[s[p]] = -noarch.PtrdiffT((w[s[p]])) - noarch.PtrdiffT(2/8)
		}
	}
	{
		// unmark node k
		w[k] = -noarch.PtrdiffT((w[k])) - noarch.PtrdiffT(2/8)
	}
	// s [top..n-1] contains pattern of L(k,:)
	return noarch.PtrdiffT((top))
}

// cs_etree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_etree.c:3
// compute the etree of A (using triu(A), or A'A without forming A'A
func cs_etree(A []cs, ata noarch.PtrdiffT) []noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var k noarch.PtrdiffT
	var p noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var inext noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var parent []noarch.PtrdiffT
	var ancestor []noarch.PtrdiffT
	var prev []noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	// allocate result
	parent = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	// get workspace
	w = cs_malloc(n+noarch.PtrdiffT(func() int32 {
		if bool(noarch.PtrdiffT(ata)) {
			return int32(noarch.PtrdiffT(m))
		}
		return 0
	}()/8), uint(0)).([]noarch.PtrdiffT)
	if w == nil || parent == nil {
		return (cs_idone(parent, nil, w, 0))
	}
	ancestor = w
	prev = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
	if bool(noarch.PtrdiffT(ata)) {
		for i = 0; i < m; i++ {
			prev[i] = noarch.PtrdiffT(-1)
		}
	}
	for k = 0; k < n; k++ {
		// node k has no parent yet
		parent[k] = noarch.PtrdiffT(-1)
		// nor does k have an ancestor
		ancestor[k] = noarch.PtrdiffT(-1)
		for p = Ap[k]; p < Ap[k+noarch.PtrdiffT(1/8)]; p++ {
			i = noarch.PtrdiffT(func() int32 {
				if bool(noarch.PtrdiffT(ata)) {
					return int32(noarch.PtrdiffT((prev[Ai[p]])))
				}
				return int32(noarch.PtrdiffT((Ai[p])))
			}() / 8)
			{
				// traverse from i to k
				for ; i != noarch.PtrdiffT(int32(-1)/8) && i < k; i = inext {
					// inext = ancestor of i
					inext = ancestor[i]
					// path compression
					ancestor[i] = k
					if inext == noarch.PtrdiffT(int32(-1)/8) {
						// no anc., parent is k
						parent[i] = k
					}
				}
			}
			if bool(noarch.PtrdiffT(ata)) {
				prev[Ai[p]] = k
			}
		}
	}
	return (cs_idone(parent, nil, w, 1))
}

// cs_fkeep - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_fkeep.c:3
// drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1
func cs_fkeep(A []cs, fkeep func(noarch.PtrdiffT, noarch.PtrdiffT, float64, interface{}) noarch.PtrdiffT, other interface{}) noarch.PtrdiffT {
	var j noarch.PtrdiffT
	var p noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Ax []float64
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || fkeep == nil {
		// check inputs
		return noarch.PtrdiffT((-1))
	}
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	for j = 0; j < n; j++ {
		// get current location of col j
		p = Ap[j]
		// record new location of col j
		Ap[j] = nz
		for ; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			if bool(fkeep(noarch.PtrdiffT(Ai[p]), noarch.PtrdiffT(j), func() float64 {
				if Ax != nil {
					return Ax[p]
				}
				return 1
			}(), other)) {
				if Ax != nil {
					// keep A(i,j)
					Ax[nz] = Ax[p]
				}
				Ai[func() noarch.PtrdiffT {
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
	return noarch.PtrdiffT((nz))
}

// cs_gaxpy - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_gaxpy.c:3
// y = A*x+y
func cs_gaxpy(A []cs, x []float64, y []float64) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Ax []float64
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || x == nil || y == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	for j = 0; j < n; j++ {
		for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			y[Ai[p]] += Ax[p] * x[j]
		}
	}
	return noarch.PtrdiffT((1))
}

// cs_happly - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_happly.c:3
// apply the ith Householder vector to x
func cs_happly(V []cs, i noarch.PtrdiffT, beta float64, x []float64) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var Vp []noarch.PtrdiffT
	var Vi []noarch.PtrdiffT
	var Vx []float64
	var tau float64
	if !(V != nil && noarch.PtrdiffT(V[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || x == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	Vp = V[0].p
	Vi = V[0].i
	Vx = V[0].x
	{
		// tau = v'*x
		for p = Vp[i]; p < Vp[i+noarch.PtrdiffT(1/8)]; p++ {
			tau += Vx[p] * x[Vi[p]]
		}
	}
	// tau = beta*(v'*x)
	tau *= beta
	{
		// x = x - v*tau
		for p = Vp[i]; p < Vp[i+noarch.PtrdiffT(1/8)]; p++ {
			x[Vi[p]] -= Vx[p] * tau
		}
	}
	return noarch.PtrdiffT((1))
}

// cs_house - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_house.c:4
// create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
// * where (I-beta*v*v')*x = s*e1.  See Algo 5.1.1, Golub & Van Loan, 3rd ed.
func cs_house(x []float64, beta []float64, n noarch.PtrdiffT) float64 {
	var s float64
	var sigma float64
	var i noarch.PtrdiffT
	if x == nil || beta == nil {
		// check inputs
		return float64((-1))
	}
	for i = 1; i < n; i++ {
		sigma += x[i] * x[i]
	}
	if sigma == 0 {
		// s = |x(0)|
		s = math.Abs(x[0])
		beta[0] = float64(func() int {
			if x[0] <= 0 {
				return 2
			}
			return 0
		}())
		x[0] = 1
	} else {
		// s = norm (x)
		s = math.Sqrt(x[0]*x[0] + sigma)
		x[0] = func() float64 {
			if x[0] <= 0 {
				return (x[0] - s)
			}
			return (-sigma / (x[0] + s))
		}()
		beta[0] = -1 / (s * x[0])
	}
	return (s)
}

// cs_ipvec - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_ipvec.c:3
// x(p) = b, for dense vectors x and b; p=NULL denotes identity
func cs_ipvec(p []noarch.PtrdiffT, b []float64, x []float64, n noarch.PtrdiffT) noarch.PtrdiffT {
	var k noarch.PtrdiffT
	if x == nil || b == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	for k = 0; k < n; k++ {
		x[func() int32 {
			if p != nil {
				return int32(noarch.PtrdiffT(p[k]))
			}
			return int32(noarch.PtrdiffT(k))
		}()] = b[k]
	}
	return noarch.PtrdiffT((1))
}

// cs_leaf - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_leaf.c:3
// consider A(i,j), node j in ith row subtree and return lca(jprev,j)
func cs_leaf(i noarch.PtrdiffT, j noarch.PtrdiffT, first []noarch.PtrdiffT, maxfirst []noarch.PtrdiffT, prevleaf []noarch.PtrdiffT, ancestor []noarch.PtrdiffT, jleaf []noarch.PtrdiffT) noarch.PtrdiffT {
	var q noarch.PtrdiffT
	var s noarch.PtrdiffT
	var sparent noarch.PtrdiffT
	var jprev noarch.PtrdiffT
	if first == nil || maxfirst == nil || prevleaf == nil || ancestor == nil || jleaf == nil {
		return noarch.PtrdiffT((-1))
	}
	jleaf[0] = 0
	if i <= j || first[j] <= maxfirst[i] {
		// j not a leaf
		return noarch.PtrdiffT((-1))
	}
	// update max first[j] seen so far
	maxfirst[i] = first[j]
	// jprev = previous leaf of ith subtree
	jprev = prevleaf[i]
	prevleaf[i] = j
	// j is first or subsequent leaf
	jleaf[0] = noarch.PtrdiffT(func() int {
		if jprev == noarch.PtrdiffT(int32(-1)/8) {
			return 1
		}
		return 2
	}())
	if jleaf[0] == noarch.PtrdiffT(1/8) {
		// if 1st leaf, q = root of ith subtree
		return noarch.PtrdiffT((i))
	}
	for q = jprev; q != ancestor[q]; q = ancestor[q] {
	}
	for s = jprev; s != q; s = sparent {
		// path compression
		sparent = ancestor[s]
		ancestor[s] = q
	}
	// q = least common ancester (jprev,j)
	return noarch.PtrdiffT((q))
}

// cs_load - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_load.c:3
// load a triplet matrix from a file
func cs_load(f *noarch.File) []cs {
	var i float64
	var j float64
	var x float64
	var T []cs
	if f == nil {
		// use double for integers to avoid csi conflicts
		// check inputs
		return nil
	}
	// allocate result
	T = cs_spalloc(0, 0, 1, 1, 1)
	for noarch.Fscanf(f, []byte("%lg %lg %lg\n\x00"), (*[100000000]float64)(unsafe.Pointer(&i))[:], (*[100000000]float64)(unsafe.Pointer(&j))[:], (*[100000000]float64)(unsafe.Pointer(&x))[:]) == 3 {
		if bool(noarch.NotNoarch.PtrdiffT(cs_entry(T, noarch.PtrdiffT(i), noarch.PtrdiffT(j), x))) {
			return (cs_spfree(T))
		}
	}
	return (T)
}

// cs_lsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_lsolve.c:3
// solve Lx=b where x and b are dense.  x=b on input, solution on output.
func cs_lsolve(L []cs, x []float64) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Lp []noarch.PtrdiffT
	var Li []noarch.PtrdiffT
	var Lx []float64
	if !(L != nil && noarch.PtrdiffT(L[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || x == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(L[0].n)
	Lp = L[0].p
	Li = L[0].i
	Lx = L[0].x
	for j = 0; j < n; j++ {
		x[j] /= Lx[Lp[j]]
		for p = Lp[j] + noarch.PtrdiffT(1/8); p < Lp[j+noarch.PtrdiffT(1/8)]; p++ {
			x[Li[p]] -= Lx[p] * x[j]
		}
	}
	return noarch.PtrdiffT((1))
}

// cs_ltsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_ltsolve.c:3
// solve L'x=b where x and b are dense.  x=b on input, solution on output.
func cs_ltsolve(L []cs, x []float64) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Lp []noarch.PtrdiffT
	var Li []noarch.PtrdiffT
	var Lx []float64
	if !(L != nil && noarch.PtrdiffT(L[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || x == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(L[0].n)
	Lp = L[0].p
	Li = L[0].i
	Lx = L[0].x
	for j = n - noarch.PtrdiffT(1/8); j >= 0; j-- {
		for p = Lp[j] + noarch.PtrdiffT(1/8); p < Lp[j+noarch.PtrdiffT(1/8)]; p++ {
			x[j] -= Lx[p] * x[Li[p]]
		}
		x[j] /= Lx[Lp[j]]
	}
	return noarch.PtrdiffT((1))
}

// cs_lu - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_lu.c:3
// [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guess
func cs_lu(A []cs, S []css, tol float64) []csn {
	var L []cs
	var U []cs
	var N []csn
	var pivot float64
	var Lx []float64
	var Ux []float64
	var x []float64
	var a float64
	var t float64
	var Lp []noarch.PtrdiffT
	var Li []noarch.PtrdiffT
	var Up []noarch.PtrdiffT
	var Ui []noarch.PtrdiffT
	var pinv []noarch.PtrdiffT
	var xi []noarch.PtrdiffT
	var q []noarch.PtrdiffT
	var n noarch.PtrdiffT
	var ipiv noarch.PtrdiffT
	var k noarch.PtrdiffT
	var top noarch.PtrdiffT
	var p noarch.PtrdiffT
	var i noarch.PtrdiffT
	var col noarch.PtrdiffT
	var lnz noarch.PtrdiffT
	var unz noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || S == nil {
		// check inputs
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	q = S[0].q
	lnz = noarch.PtrdiffT(S[0].lnz)
	unz = noarch.PtrdiffT(S[0].unz)
	// get double workspace
	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
	// get csi workspace
	xi = cs_malloc(noarch.PtrdiffT(2*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
	// allocate result
	N = cs_calloc(1, uint(32)).([]csn)
	if x == nil || xi == nil || N == nil {
		return (cs_ndone(N, nil, xi, x, 0))
	}
	L = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(lnz), 1, 0)
	// allocate result L
	N[0].L = L
	U = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(unz), 1, 0)
	// allocate result U
	N[0].U = U
	pinv = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	// allocate result pinv
	N[0].pinv = pinv
	if L == nil || U == nil || pinv == nil {
		return (cs_ndone(N, nil, xi, x, 0))
	}
	Lp = L[0].p
	Up = U[0].p
	{
		// clear workspace
		for i = 0; i < n; i++ {
			x[i] = 0
		}
	}
	{
		// no rows pivotal yet
		for i = 0; i < n; i++ {
			pinv[i] = noarch.PtrdiffT(-1)
		}
	}
	{
		// no cols of L yet
		for k = 0; k <= n; k++ {
			Lp[k] = 0
		}
	}
	unz = 0
	lnz = unz
	{
		// compute L(:,k) and U(:,k)
		for k = 0; k < n; k++ {
			// --- Triangular solve ---------------------------------------------
			// L(:,k) starts here
			Lp[k] = lnz
			// U(:,k) starts here
			Up[k] = unz
			if lnz+n > noarch.PtrdiffT(L[0].nzmax) && bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(L, noarch.PtrdiffT((2*int32(noarch.PtrdiffT(L[0].nzmax))+int32(n))/8)))) || unz+n > noarch.PtrdiffT(U[0].nzmax) && bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(U, noarch.PtrdiffT((2*int32(noarch.PtrdiffT(U[0].nzmax))+int32(n))/8)))) {
				return (cs_ndone(N, nil, xi, x, 0))
			}
			Li = L[0].i
			Lx = L[0].x
			Ui = U[0].i
			Ux = U[0].x
			col = noarch.PtrdiffT(func() int32 {
				if q != nil {
					return int32(noarch.PtrdiffT((q[k])))
				}
				return int32(noarch.PtrdiffT(k))
			}() / 8)
			// x = L\A(:,col)
			top = cs_spsolve(L, A, noarch.PtrdiffT(col), xi, x, pinv, 1)
			// --- Find pivot ---------------------------------------------------
			ipiv = noarch.PtrdiffT(-1)
			a = float64(-1)
			for p = top; p < n; p++ {
				// x(i) is nonzero
				i = xi[p]
				if pinv[i] < noarch.PtrdiffT(0/8) {
					if (func() float64 {
						t = math.Abs(x[i])
						return t
					}()) > a {
						// row i is not yet pivotal
						// largest pivot candidate so far
						a = t
						ipiv = i
					}
				} else {
					// x(i) is the entry U(pinv[i],k)
					Ui[unz] = pinv[i]
					Ux[func() noarch.PtrdiffT {
						defer func() {
							unz++
						}()
						return unz
					}()] = x[i]
				}
			}
			if ipiv == noarch.PtrdiffT(int32(-1)/8) || a <= 0 {
				return (cs_ndone(N, nil, xi, x, 0))
			}
			if pinv[col] < noarch.PtrdiffT(0/8) && math.Abs(x[col]) >= a*tol {
				// tol=1 for  partial pivoting; tol<1 gives preference to diagonal
				ipiv = col
			}
			// --- Divide by pivot ----------------------------------------------
			// the chosen pivot
			pivot = x[ipiv]
			// last entry in U(:,k) is U(k,k)
			Ui[unz] = k
			Ux[func() noarch.PtrdiffT {
				defer func() {
					unz++
				}()
				return unz
			}()] = pivot
			// ipiv is the kth pivot row
			pinv[ipiv] = k
			// first entry in L(:,k) is L(k,k) = 1
			Li[lnz] = ipiv
			Lx[func() noarch.PtrdiffT {
				defer func() {
					lnz++
				}()
				return lnz
			}()] = 1
			{
				// L(k+1:n,k) = x / pivot
				for p = top; p < n; p++ {
					i = xi[p]
					if pinv[i] < noarch.PtrdiffT(0/8) {
						// x(i) is an entry in L(:,k)
						// save unpermuted row in L
						Li[lnz] = i
						// scale pivot column
						Lx[func() noarch.PtrdiffT {
							defer func() {
								lnz++
							}()
							return lnz
						}()] = x[i] / pivot
					}
					// x [0..n-1] = 0 for next k
					x[i] = 0
				}
			}
		}
	}
	// --- Finalize L and U -------------------------------------------------
	Lp[n] = lnz
	Up[n] = unz
	// fix row indices of L for final pinv
	Li = L[0].i
	for p = 0; p < lnz; p++ {
		Li[p] = pinv[Li[p]]
	}
	// remove extra space from L and U
	cs_sprealloc(L, 0)
	cs_sprealloc(U, 0)
	// success
	return (cs_ndone(N, nil, xi, x, 1))
}

// cs_lusol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_lusol.c:3
// x=A\b where A is unsymmetric; b overwritten with solution
func cs_lusol(order noarch.PtrdiffT, A []cs, b []float64, tol float64) noarch.PtrdiffT {
	var x []float64
	var S []css
	var N []csn
	var n noarch.PtrdiffT
	var ok noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || b == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(A[0].n)
	// ordering and symbolic analysis
	S = cs_sqr(noarch.PtrdiffT(order), A, 0)
	// numeric LU factorization
	N = cs_lu(A, S, tol)
	// get workspace
	x = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
	ok = noarch.PtrdiffT(S != nil && N != nil && x != nil)
	if bool(noarch.PtrdiffT(ok)) {
		// x = b(p)
		cs_ipvec(N[0].pinv, b, x, noarch.PtrdiffT(n))
		// x = L\x
		cs_lsolve(N[0].L, x)
		// x = U\x
		cs_usolve(N[0].U, x)
		// b(q) = x
		cs_ipvec(S[0].q, x, b, noarch.PtrdiffT(n))
	}
	cs_free(x)
	cs_sfree(S)
	cs_nfree(N)
	return noarch.PtrdiffT((ok))
}

// cs_malloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_malloc.c:10
// wrapper for malloc
func cs_malloc(n noarch.PtrdiffT, size uint) interface{} {
	return (make([]byte, uint32(func() int32 {
		if n > noarch.PtrdiffT(int32(1)/8) {
			return int32(noarch.PtrdiffT((n)))
		}
		return int32((1))
	}())*uint32(size)))
}

// cs_calloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_malloc.c:16
// wrapper for calloc
func cs_calloc(n noarch.PtrdiffT, size uint) interface{} {
	return (make([]byte, (size)*(uint(uint32((func() int32 {
		if n > noarch.PtrdiffT(int32(1)/8) {
			return int32(noarch.PtrdiffT((n)))
		}
		return int32((1))
	}()))))))
}

// cs_free - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_malloc.c:22
// wrapper for free
func cs_free(p interface{}) interface{} {
	if p != nil {
		_ = p
		// free p if it is not already NULL
	}
	// return NULL to simplify the use of cs_free
	return nil
}

// cs_realloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_malloc.c:29
// wrapper for realloc
func cs_realloc(p interface{}, n noarch.PtrdiffT, size uint, ok []noarch.PtrdiffT) interface{} {
	var pnew interface{}
	// realloc the block
	pnew = make([]byte, uint32(func() int32 {
		if n > noarch.PtrdiffT(int32(1)/8) {
			return int32(noarch.PtrdiffT((n)))
		}
		return int32((1))
	}())*uint32(size)*1/1)
	// realloc fails if pnew is NULL
	ok[0] = noarch.PtrdiffT(pnew != nil)
	// return original p if failure
	return (func() interface{} {
		if bool(noarch.PtrdiffT((ok[0]))) {
			return pnew
		}
		return p
	}())
}

// cs_augment - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_maxtrans.c:3
// find an augmenting path starting at column k and extend the match if found
func cs_augment(k noarch.PtrdiffT, A []cs, jmatch []noarch.PtrdiffT, cheap []noarch.PtrdiffT, w []noarch.PtrdiffT, js []noarch.PtrdiffT, is []noarch.PtrdiffT, ps []noarch.PtrdiffT) {
	var found noarch.PtrdiffT
	var p noarch.PtrdiffT
	var i noarch.PtrdiffT = noarch.PtrdiffT(-1)
	var Ap []noarch.PtrdiffT = A[0].p
	var Ai []noarch.PtrdiffT = A[0].i
	var head noarch.PtrdiffT
	var j noarch.PtrdiffT
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
			for p = cheap[j]; p < Ap[j+noarch.PtrdiffT(1/8)] && bool(noarch.NotNoarch.PtrdiffT(noarch.PtrdiffT(found))); p++ {
				// try a cheap assignment (i,j)
				i = Ai[p]
				found = noarch.PtrdiffT(jmatch[i] == noarch.PtrdiffT(int32(-1)/8))
			}
			// start here next time j is traversed
			cheap[j] = p
			if bool(noarch.PtrdiffT(found)) {
				// column j matched with row i
				is[head] = i
				// end of augmenting path
				break
			}
			// no cheap match: start dfs for j
			ps[head] = Ap[j]
		}
		{
			// --- Depth-first-search of neighbors of j -------------------------
			for p = ps[head]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
				// consider row i
				i = Ai[p]
				if w[jmatch[i]] == k {
					// skip jmatch [i] if marked
					continue
				}
				// pause dfs of node j
				ps[head] = p + noarch.PtrdiffT(1/8)
				// i will be matched with j if found
				is[head] = i
				// start dfs at column jmatch [i]
				js[func() noarch.PtrdiffT {
					head++
					return head
				}()] = jmatch[i]
				break
			}
		}
		if p == Ap[j+noarch.PtrdiffT(1/8)] {
			// node j is done; pop from stack
			head--
		}
	}
	if bool(noarch.PtrdiffT(found)) {
		// augment the match if path found:
		for p = head; p >= 0; p-- {
			jmatch[is[p]] = js[p]
		}
	}
}

// cs_maxtrans - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_maxtrans.c:44
// find a maximum transveral
//[jmatch [0..m-1]; imatch [0..n-1]]
func cs_maxtrans(A []cs, seed noarch.PtrdiffT) []noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var j noarch.PtrdiffT
	var k noarch.PtrdiffT
	var n noarch.PtrdiffT
	var m noarch.PtrdiffT
	var p noarch.PtrdiffT
	var n2 noarch.PtrdiffT
	var m2 noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var jimatch []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var cheap []noarch.PtrdiffT
	var js []noarch.PtrdiffT
	var is []noarch.PtrdiffT
	var ps []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var jmatch []noarch.PtrdiffT
	var imatch []noarch.PtrdiffT
	var q []noarch.PtrdiffT
	var C []cs
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	m = noarch.PtrdiffT(A[0].m)
	Ap = A[0].p
	Ai = A[0].i
	jimatch = cs_calloc(m+n, uint(0)).([]noarch.PtrdiffT)
	// allocate result
	w = jimatch
	if jimatch == nil {
		return nil
	}
	{
		// count nonempty rows and columns
		k = 0
		j = 0
		for j = 0; j < n; j++ {
			n2 += noarch.PtrdiffT(int32(map[bool]int{false: 0, true: 1}[Ap[j] < Ap[j+noarch.PtrdiffT(1/8)]]) / 8)
			for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
				w[Ai[p]] = 1
				// count entries already on diagonal
				k += noarch.PtrdiffT(int32(map[bool]int{false: 0, true: 1}[j == Ai[p]]) / 8)
			}
		}
	}
	if k == noarch.PtrdiffT(func() int32 {
		if m < n {
			return int32(noarch.PtrdiffT((m)))
		}
		return int32(noarch.PtrdiffT((n)))
	}()/8) {
		// quick return if diagonal zero-free
		jmatch = jimatch
		imatch = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&jimatch[0])) + (uintptr)(int(m))*unsafe.Sizeof(jimatch[0]))))[:]
		for i = 0; i < k; i++ {
			jmatch[i] = i
		}
		for ; i < m; i++ {
			jmatch[i] = noarch.PtrdiffT(-1)
		}
		for j = 0; j < k; j++ {
			imatch[j] = j
		}
		for ; j < n; j++ {
			imatch[j] = noarch.PtrdiffT(-1)
		}
		return (cs_idone(jimatch, nil, nil, 1))
	}
	for i = 0; i < m; i++ {
		m2 += w[i]
	}
	// transpose if needed
	C = func() []cs {
		if m2 < n2 {
			return cs_transpose(A, 0)
		}
		return (A)
	}()
	if C == nil {
		return (cs_idone(jimatch, func() []cs {
			if m2 < n2 {
				return C
			}
			return nil
		}(), nil, 0))
	}
	n = noarch.PtrdiffT(C[0].n)
	m = noarch.PtrdiffT(C[0].m)
	Cp = C[0].p
	jmatch = func() []noarch.PtrdiffT {
		if m2 < n2 {
			return (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&jimatch[0])) + (uintptr)(int(n))*unsafe.Sizeof(jimatch[0]))))[:]
		}
		return jimatch
	}()
	imatch = func() []noarch.PtrdiffT {
		if m2 < n2 {
			return jimatch
		}
		return (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&jimatch[0])) + (uintptr)(int(m))*unsafe.Sizeof(jimatch[0]))))[:]
	}()
	// get workspace
	w = cs_malloc(noarch.PtrdiffT(5*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
	if w == nil {
		return (cs_idone(jimatch, func() []cs {
			if m2 < n2 {
				return C
			}
			return nil
		}(), w, 0))
	}
	cheap = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
	js = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	is = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(3*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	ps = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(4*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	{
		// for cheap assignment
		for j = 0; j < n; j++ {
			cheap[j] = Cp[j]
		}
	}
	{
		// all columns unflagged
		for j = 0; j < n; j++ {
			w[j] = noarch.PtrdiffT(-1)
		}
	}
	{
		// nothing matched yet
		for i = 0; i < m; i++ {
			jmatch[i] = noarch.PtrdiffT(-1)
		}
	}
	// q = random permutation
	q = cs_randperm(noarch.PtrdiffT(n), noarch.PtrdiffT(seed))
	{
		// augment, starting at column q[k]
		for k = 0; k < n; k++ {
			cs_augment(noarch.PtrdiffT(func() int32 {
				if q != nil {
					return int32(noarch.PtrdiffT(q[k]))
				}
				return int32(noarch.PtrdiffT(k))
			}()/8), C, jmatch, cheap, w, js, is, ps)
		}
	}
	cs_free(q)
	{
		// find row match
		for j = 0; j < n; j++ {
			imatch[j] = noarch.PtrdiffT(-1)
		}
	}
	for i = 0; i < m; i++ {
		if jmatch[i] >= 0 {
			imatch[jmatch[i]] = i
		}
	}
	return (cs_idone(jimatch, func() []cs {
		if m2 < n2 {
			return C
		}
		return nil
	}(), w, 1))
}

// cs_multiply - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_multiply.c:3
// C = A*B
func cs_multiply(A []cs, B []cs) []cs {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var anz noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var Bp []noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var bnz noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var values noarch.PtrdiffT
	var Bi []noarch.PtrdiffT
	var x []float64
	var Bx []float64
	var Cx []float64
	var C []cs
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || !(B != nil && noarch.PtrdiffT(B[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	if noarch.PtrdiffT(A[0].n) != noarch.PtrdiffT(B[0].m) {
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	anz = A[0].p[noarch.PtrdiffT(A[0].n)]
	n = noarch.PtrdiffT(B[0].n)
	Bp = B[0].p
	Bi = B[0].i
	Bx = B[0].x
	bnz = Bp[n]
	// get workspace
	w = cs_calloc(noarch.PtrdiffT(m), uint(0)).([]noarch.PtrdiffT)
	values = noarch.PtrdiffT(A[0].x != nil && Bx != nil)
	// get workspace
	x = func() interface{} {
		if bool(noarch.PtrdiffT(values)) {
			return cs_malloc(noarch.PtrdiffT(m), uint(8))
		}
		return nil
	}().([]float64)
	// allocate result
	C = cs_spalloc(noarch.PtrdiffT(m), noarch.PtrdiffT(n), anz+bnz, noarch.PtrdiffT(values), 0)
	if C == nil || w == nil || bool(values) && x == nil {
		return (cs_done(C, w, x, 0))
	}
	Cp = C[0].p
	for j = 0; j < n; j++ {
		if nz+m > noarch.PtrdiffT(C[0].nzmax) && bool(noarch.NotNoarch.PtrdiffT(cs_sprealloc(C, noarch.PtrdiffT((2*int32(C[0].nzmax)+int32(m))/8)))) {
			// out of memory
			return (cs_done(C, w, x, 0))
		}
		// C->i and C->x may be reallocated
		Ci = C[0].i
		Cx = C[0].x
		// column j of C starts here
		Cp[j] = nz
		for p = Bp[j]; p < Bp[j+noarch.PtrdiffT(1/8)]; p++ {
			nz = cs_scatter(A, noarch.PtrdiffT(Bi[p]), func() float64 {
				if Bx != nil {
					return Bx[p]
				}
				return 1
			}(), w, x, j+noarch.PtrdiffT(1/8), C, noarch.PtrdiffT(nz))
		}
		if bool(noarch.PtrdiffT(values)) {
			for p = Cp[j]; p < nz; p++ {
				Cx[p] = x[Ci[p]]
			}
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	// remove extra space from C
	cs_sprealloc(C, 0)
	// success; free workspace, return C
	return (cs_done(C, w, x, 1))
}

// cs_norm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_norm.c:3
// 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum
func cs_norm(A []cs) float64 {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ax []float64
	var norm float64
	var s float64
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || A[0].x == nil {
		// check inputs
		return float64((-1))
	}
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ax = A[0].x
	for j = 0; j < n; j++ {
		{
			s = 0
			p = Ap[j]
			for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
				s += math.Abs(Ax[p])
			}
		}
		norm = func() float64 {
			if norm > s {
				return (norm)
			}
			return (s)
		}()
	}
	return (norm)
}

// cs_permute - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_permute.c:3
// C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1.
func cs_permute(A []cs, pinv []noarch.PtrdiffT, q []noarch.PtrdiffT, values noarch.PtrdiffT) []cs {
	var t noarch.PtrdiffT
	var j noarch.PtrdiffT
	var k noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var Cx []float64
	var Ax []float64
	var C []cs
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	// alloc result
	C = cs_spalloc(noarch.PtrdiffT(m), noarch.PtrdiffT(n), noarch.PtrdiffT(Ap[n]), noarch.PtrdiffT(bool(values) && Ax != nil), 0)
	if C == nil {
		// out of memory
		return (cs_done(C, nil, nil, 0))
	}
	Cp = C[0].p
	Ci = C[0].i
	Cx = C[0].x
	for k = 0; k < n; k++ {
		// column k of C is column q[k] of A
		Cp[k] = nz
		j = noarch.PtrdiffT(func() int32 {
			if q != nil {
				return int32(noarch.PtrdiffT((q[k])))
			}
			return int32(noarch.PtrdiffT(k))
		}() / 8)
		for t = Ap[j]; t < Ap[j+noarch.PtrdiffT(1/8)]; t++ {
			if Cx != nil {
				// row i of A is row pinv[i] of C
				Cx[nz] = Ax[t]
			}
			Ci[func() noarch.PtrdiffT {
				defer func() {
					nz++
				}()
				return nz
			}()] = noarch.PtrdiffT(func() int32 {
				if pinv != nil {
					return int32(noarch.PtrdiffT((pinv[Ai[t]])))
				}
				return int32(noarch.PtrdiffT(Ai[t]))
			}() / 8)
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	return (cs_done(C, nil, nil, 1))
}

// cs_pinv - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_pinv.c:3
// pinv = p', or p = pinv'
func cs_pinv(p []noarch.PtrdiffT, n noarch.PtrdiffT) []noarch.PtrdiffT {
	var k noarch.PtrdiffT
	var pinv []noarch.PtrdiffT
	if p == nil {
		// p = NULL denotes identity
		return nil
	}
	// allocate result
	pinv = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	if pinv == nil {
		// out of memory
		return nil
	}
	{
		// invert the permutation
		for k = 0; k < n; k++ {
			pinv[p[k]] = k
		}
	}
	// return result
	return (pinv)
}

// cs_post - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_post.c:3
// post order a forest
func cs_post(parent []noarch.PtrdiffT, n noarch.PtrdiffT) []noarch.PtrdiffT {
	var j noarch.PtrdiffT
	var k noarch.PtrdiffT
	var post []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var head []noarch.PtrdiffT
	var next []noarch.PtrdiffT
	var stack []noarch.PtrdiffT
	if parent == nil {
		// check inputs
		return nil
	}
	// allocate result
	post = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	// get workspace
	w = cs_malloc(noarch.PtrdiffT(3*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
	if w == nil || post == nil {
		return (cs_idone(post, nil, w, 0))
	}
	head = w
	next = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(n))*unsafe.Sizeof(w[0]))))[:]
	stack = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof(w[0]))))[:]
	{
		// empty linked lists
		for j = 0; j < n; j++ {
			head[j] = noarch.PtrdiffT(-1)
		}
	}
	{
		// traverse nodes in reverse order
		for j = n - noarch.PtrdiffT(1/8); j >= 0; j-- {
			if parent[j] == noarch.PtrdiffT(int32(-1)/8) {
				// j is a root
				continue
			}
			// add j to list of its parent
			next[j] = head[parent[j]]
			head[parent[j]] = j
		}
	}
	for j = 0; j < n; j++ {
		if parent[j] != noarch.PtrdiffT(int32(-1)/8) {
			// skip j if it is not a root
			continue
		}
		k = cs_tdfs(noarch.PtrdiffT(j), noarch.PtrdiffT(k), head, next, post, stack)
	}
	// success; free w, return post
	return (cs_idone(post, nil, w, 1))
}

// cs_print - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_print.c:3
// print a sparse matrix; use %g for integers to avoid differences with csi
func cs_print(A []cs, brief noarch.PtrdiffT) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var nzmax noarch.PtrdiffT
	var nz noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Ax []float64
	if A == nil {
		fmt.Printf("(null)\n")
		return noarch.PtrdiffT((0))
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	nzmax = noarch.PtrdiffT(A[0].nzmax)
	nz = noarch.PtrdiffT(A[0].nz)
	noarch.Printf([]byte("CSparse Version %d.%d.%d, %s.  %s\n\x00"), 3, 2, 0, []byte("Sept 12, 2017\x00"), []byte("Copyright (c) Timothy A. Davis, 2006-2016\x00"))
	if nz < noarch.PtrdiffT(0/8) {
		noarch.Printf([]byte("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n\x00"), float64(noarch.PtrdiffT(m)), float64(noarch.PtrdiffT(n)), float64(noarch.PtrdiffT(nzmax)), float64(noarch.PtrdiffT((Ap[n]))), cs_norm(A))
		for j = 0; j < n; j++ {
			noarch.Printf([]byte("    col %g : locations %g to %g\n\x00"), float64(noarch.PtrdiffT(j)), float64(noarch.PtrdiffT((Ap[j]))), float64((int32(Ap[j+noarch.PtrdiffT(1/8)] - noarch.PtrdiffT(1/8)))))
			for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
				noarch.Printf([]byte("      %g : %g\n\x00"), float64(noarch.PtrdiffT((Ai[p]))), func() float64 {
					if Ax != nil {
						return Ax[p]
					}
					return 1
				}())
				if bool(brief) && p > noarch.PtrdiffT(20/8) {
					fmt.Printf("  ...\n")
					return noarch.PtrdiffT((1))
				}
			}
		}
	} else {
		noarch.Printf([]byte("triplet: %g-by-%g, nzmax: %g nnz: %g\n\x00"), float64(noarch.PtrdiffT(m)), float64(noarch.PtrdiffT(n)), float64(noarch.PtrdiffT(nzmax)), float64(noarch.PtrdiffT(nz)))
		for p = 0; p < nz; p++ {
			noarch.Printf([]byte("    %g %g : %g\n\x00"), float64(noarch.PtrdiffT((Ai[p]))), float64(noarch.PtrdiffT((Ap[p]))), func() float64 {
				if Ax != nil {
					return Ax[p]
				}
				return 1
			}())
			if bool(brief) && p > noarch.PtrdiffT(20/8) {
				fmt.Printf("  ...\n")
				return noarch.PtrdiffT((1))
			}
		}
	}
	return noarch.PtrdiffT((1))
}

// cs_pvec - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_pvec.c:3
// x = b(p), for dense vectors x and b; p=NULL denotes identity
func cs_pvec(p []noarch.PtrdiffT, b []float64, x []float64, n noarch.PtrdiffT) noarch.PtrdiffT {
	var k noarch.PtrdiffT
	if x == nil || b == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	for k = 0; k < n; k++ {
		x[k] = b[func() int32 {
			if p != nil {
				return int32(noarch.PtrdiffT(p[k]))
			}
			return int32(noarch.PtrdiffT(k))
		}()]
	}
	return noarch.PtrdiffT((1))
}

// cs_qr - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_qr.c:3
// sparse QR factorization [V,beta,pinv,R] = qr (A)
func cs_qr(A []cs, S []css) []csn {
	var Rx []float64
	var Vx []float64
	var Ax []float64
	var x []float64
	var Beta []float64
	var i noarch.PtrdiffT
	var k noarch.PtrdiffT
	var p noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var vnz noarch.PtrdiffT
	var p1 noarch.PtrdiffT
	var top noarch.PtrdiffT
	var m2 noarch.PtrdiffT
	var len noarch.PtrdiffT
	var col noarch.PtrdiffT
	var rnz noarch.PtrdiffT
	var s []noarch.PtrdiffT
	var leftmost []noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var parent []noarch.PtrdiffT
	var Rp []noarch.PtrdiffT
	var Ri []noarch.PtrdiffT
	var Vp []noarch.PtrdiffT
	var Vi []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var pinv []noarch.PtrdiffT
	var q []noarch.PtrdiffT
	var R []cs
	var V []cs
	var N []csn
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || S == nil {
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	q = S[0].q
	parent = S[0].parent
	pinv = S[0].pinv
	m2 = noarch.PtrdiffT(S[0].m2)
	vnz = noarch.PtrdiffT(S[0].lnz)
	rnz = noarch.PtrdiffT(S[0].unz)
	leftmost = S[0].leftmost
	// get csi workspace
	w = cs_malloc(m2+n, uint(0)).([]noarch.PtrdiffT)
	// get double workspace
	x = cs_malloc(noarch.PtrdiffT(m2), uint(8)).([]float64)
	// allocate result
	N = cs_calloc(1, uint(32)).([]csn)
	if w == nil || x == nil || N == nil {
		return (cs_ndone(N, nil, w, x, 0))
	}
	// s is size n
	s = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m2))*unsafe.Sizeof(w[0]))))[:]
	{
		// clear workspace x
		for k = 0; k < m2; k++ {
			x[k] = 0
		}
	}
	V = cs_spalloc(noarch.PtrdiffT(m2), noarch.PtrdiffT(n), noarch.PtrdiffT(vnz), 1, 0)
	// allocate result V
	N[0].L = V
	R = cs_spalloc(noarch.PtrdiffT(m2), noarch.PtrdiffT(n), noarch.PtrdiffT(rnz), 1, 0)
	// allocate result R
	N[0].U = R
	Beta = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
	// allocate result Beta
	N[0].B = Beta
	if R == nil || V == nil || Beta == nil {
		return (cs_ndone(N, nil, w, x, 0))
	}
	Rp = R[0].p
	Ri = R[0].i
	Rx = R[0].x
	Vp = V[0].p
	Vi = V[0].i
	Vx = V[0].x
	{
		// clear w, to mark nodes
		for i = 0; i < m2; i++ {
			w[i] = noarch.PtrdiffT(-1)
		}
	}
	rnz = 0
	vnz = 0
	{
		// compute V and R
		for k = 0; k < n; k++ {
			// R(:,k) starts here
			Rp[k] = rnz
			p1 = vnz
			// V(:,k) starts here
			Vp[k] = p1
			// add V(k,k) to pattern of V
			w[k] = k
			Vi[func() noarch.PtrdiffT {
				defer func() {
					vnz++
				}()
				return vnz
			}()] = k
			top = n
			col = noarch.PtrdiffT(func() int32 {
				if q != nil {
					return int32(noarch.PtrdiffT(q[k]))
				}
				return int32(noarch.PtrdiffT(k))
			}() / 8)
			{
				// find R(:,k) pattern
				for p = Ap[col]; p < Ap[col+noarch.PtrdiffT(1/8)]; p++ {
					// i = min(find(A(i,q)))
					i = leftmost[Ai[p]]
					{
						// traverse up to k
						for len = 0; w[i] != k; i = parent[i] {
							s[func() noarch.PtrdiffT {
								defer func() {
									len++
								}()
								return len
							}()] = i
							w[i] = k
						}
					}
					for len > noarch.PtrdiffT(0/8) {
						// push path on stack
						s[func() noarch.PtrdiffT {
							top--
							return top
						}()] = s[func() noarch.PtrdiffT {
							len--
							return len
						}()]
					}
					// i = permuted row of A(:,col)
					i = pinv[Ai[p]]
					// x (i) = A(:,col)
					x[i] = Ax[p]
					if i > k && w[i] < k {
						// pattern of V(:,k) = x (k+1:m)
						// add i to pattern of V(:,k)
						Vi[func() noarch.PtrdiffT {
							defer func() {
								vnz++
							}()
							return vnz
						}()] = i
						w[i] = k
					}
				}
			}
			{
				// for each i in pattern of R(:,k)
				for p = top; p < n; p++ {
					// R(i,k) is nonzero
					i = s[p]
					// apply (V(i),Beta(i)) to x
					cs_happly(V, noarch.PtrdiffT(i), Beta[i], x)
					// R(i,k) = x(i)
					Ri[rnz] = i
					Rx[func() noarch.PtrdiffT {
						defer func() {
							rnz++
						}()
						return rnz
					}()] = x[i]
					x[i] = 0
					if parent[i] == k {
						vnz = cs_scatter(V, noarch.PtrdiffT(i), 0, w, nil, noarch.PtrdiffT(k), V, noarch.PtrdiffT(vnz))
					}
				}
			}
			{
				// gather V(:,k) = x
				for p = p1; p < vnz; p++ {
					Vx[p] = x[Vi[p]]
					x[Vi[p]] = 0
				}
			}
			// R(k,k) = norm (x)
			Ri[rnz] = k
			// [v,beta]=house(x)
			Rx[func() noarch.PtrdiffT {
				defer func() {
					rnz++
				}()
				return rnz
			}()] = cs_house((*(*[1000000000]float64)(unsafe.Pointer(uintptr(unsafe.Pointer(&Vx[0])) + (uintptr)(int(p1))*unsafe.Sizeof(Vx[0]))))[:], (*(*[1000000000]float64)(unsafe.Pointer(uintptr(unsafe.Pointer(&Beta[0])) + (uintptr)(int(k))*unsafe.Sizeof(Beta[0]))))[:], vnz-p1)
		}
	}
	// finalize R
	Rp[n] = rnz
	// finalize V
	Vp[n] = vnz
	// success
	return (cs_ndone(N, nil, w, x, 1))
}

// cs_qrsol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_qrsol.c:3
// x=A\b where A can be rectangular; b overwritten with solution
func cs_qrsol(order noarch.PtrdiffT, A []cs, b []float64) noarch.PtrdiffT {
	var x []float64
	var S []css
	var N []csn
	var AT []cs
	var k noarch.PtrdiffT
	var m noarch.PtrdiffT
	var n noarch.PtrdiffT
	var ok noarch.PtrdiffT
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || b == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(A[0].n)
	m = noarch.PtrdiffT(A[0].m)
	if m >= n {
		// ordering and symbolic analysis
		S = cs_sqr(noarch.PtrdiffT(order), A, 1)
		// numeric QR factorization
		N = cs_qr(A, S)
		// get workspace
		x = cs_calloc(noarch.PtrdiffT(func() int32 {
			if S != nil {
				return int32(noarch.PtrdiffT(S[0].m2))
			}
			return 1
		}()/8), uint(8)).([]float64)
		ok = noarch.PtrdiffT(S != nil && N != nil && x != nil)
		if bool(noarch.PtrdiffT(ok)) {
			// x(0:m-1) = b(p(0:m-1)
			cs_ipvec(S[0].pinv, b, x, noarch.PtrdiffT(m))
			{
				// apply Householder refl. to x
				for k = 0; k < n; k++ {
					cs_happly(N[0].L, noarch.PtrdiffT(k), N[0].B[k], x)
				}
			}
			// x = R\x
			cs_usolve(N[0].U, x)
			// b(q(0:n-1)) = x(0:n-1)
			cs_ipvec(S[0].q, x, b, noarch.PtrdiffT(n))
		}
	} else {
		// Ax=b is underdetermined
		AT = cs_transpose(A, 1)
		// ordering and symbolic analysis
		S = cs_sqr(noarch.PtrdiffT(order), AT, 1)
		// numeric QR factorization of A'
		N = cs_qr(AT, S)
		// get workspace
		x = cs_calloc(noarch.PtrdiffT(func() int32 {
			if S != nil {
				return int32(noarch.PtrdiffT(S[0].m2))
			}
			return 1
		}()/8), uint(8)).([]float64)
		ok = noarch.PtrdiffT(AT != nil && S != nil && N != nil && x != nil)
		if bool(noarch.PtrdiffT(ok)) {
			// x(q(0:m-1)) = b(0:m-1)
			cs_pvec(S[0].q, b, x, noarch.PtrdiffT(m))
			// x = R'\x
			cs_utsolve(N[0].U, x)
			{
				// apply Householder refl. to x
				for k = m - noarch.PtrdiffT(1/8); k >= 0; k-- {
					cs_happly(N[0].L, noarch.PtrdiffT(k), N[0].B[k], x)
				}
			}
			// b(0:n-1) = x(p(0:n-1))
			cs_pvec(S[0].pinv, x, b, noarch.PtrdiffT(n))
		}
	}
	cs_free(x)
	cs_sfree(S)
	cs_nfree(N)
	cs_spfree(AT)
	return noarch.PtrdiffT((ok))
}

// cs_randperm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_randperm.c:5
// return a random permutation vector, the identity perm, or p = n-1:-1:0.
// * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
// * p = random permutation.
func cs_randperm(n noarch.PtrdiffT, seed noarch.PtrdiffT) []noarch.PtrdiffT {
	var p []noarch.PtrdiffT
	var k noarch.PtrdiffT
	var j noarch.PtrdiffT
	var t noarch.PtrdiffT
	if seed == noarch.PtrdiffT(0/8) {
		// return p = NULL (identity)
		return nil
	}
	// allocate result
	p = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	if p == nil {
		// out of memory
		return nil
	}
	for k = 0; k < n; k++ {
		p[k] = n - k - noarch.PtrdiffT(1/8)
	}
	if seed == noarch.PtrdiffT(int32(-1)/8) {
		// return reverse permutation
		return (p)
	}
	// get new random number seed
	rand.Seed(int64(uint32(noarch.PtrdiffT(seed))))
	for k = 0; k < n; k++ {
		// j = rand integer in range k to n-1
		j = k + noarch.PtrdiffT(int32(rand.Int())%int32(n-k)/8)
		// swap p[k] and p[j]
		t = p[j]
		p[j] = p[k]
		p[k] = t
	}
	return (p)
}

// cs_reach - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_reach.c:4
// xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
// * xi [n...2n-1] used as workspace
func cs_reach(G []cs, B []cs, k noarch.PtrdiffT, xi []noarch.PtrdiffT, pinv []noarch.PtrdiffT) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var n noarch.PtrdiffT
	var top noarch.PtrdiffT
	var Bp []noarch.PtrdiffT
	var Bi []noarch.PtrdiffT
	var Gp []noarch.PtrdiffT
	if !(G != nil && noarch.PtrdiffT(G[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || !(B != nil && noarch.PtrdiffT(B[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || xi == nil {
		// check inputs
		return noarch.PtrdiffT((-1))
	}
	n = noarch.PtrdiffT(G[0].n)
	Bp = B[0].p
	Bi = B[0].i
	Gp = G[0].p
	top = n
	for p = Bp[k]; p < Bp[k+noarch.PtrdiffT(1/8)]; p++ {
		if !(Gp[Bi[p]] < noarch.PtrdiffT(0/8)) {
			// start a dfs at unmarked node i
			top = cs_dfs(noarch.PtrdiffT(Bi[p]), G, noarch.PtrdiffT(top), xi, (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&xi[0])) + (uintptr)(int(n))*unsafe.Sizeof(xi[0]))))[:], pinv)
		}
	}
	{
		// restore G
		for p = top; p < n; p++ {
			Gp[xi[p]] = -noarch.PtrdiffT((Gp[xi[p]])) - noarch.PtrdiffT(2/8)
		}
	}
	return noarch.PtrdiffT((top))
}

// cs_scatter - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_scatter.c:3
// x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse
func cs_scatter(A *cs, j int, beta float64, w *int, x []float64, mark int, C *cs, nz int) int {
	var i noarch.PtrdiffT
	var p noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var Ax []float64
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || w == nil || !(C != nil && noarch.PtrdiffT(C[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return noarch.PtrdiffT((-1))
	}
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	Ci = C[0].i
	for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
		// A(i,j) is nonzero
		i = Ai[p]
		if w[i] < mark {
			// i is new entry in column j
			w[i] = mark
			// add i to pattern of C(:,j)
			Ci[func() noarch.PtrdiffT {
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
	return noarch.PtrdiffT((nz))
}

// cs_scc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_scc.c:3
// find the strongly connected components of a square matrix
// matrix A temporarily modified, then restored
func cs_scc(A []cs) []csd {
	var n noarch.PtrdiffT
	var i noarch.PtrdiffT
	var k noarch.PtrdiffT
	var b noarch.PtrdiffT
	var nb noarch.PtrdiffT
	var top noarch.PtrdiffT
	var xi []noarch.PtrdiffT
	var pstack []noarch.PtrdiffT
	var p []noarch.PtrdiffT
	var r []noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var ATp []noarch.PtrdiffT
	var rcopy []noarch.PtrdiffT
	var Blk []noarch.PtrdiffT
	var AT []cs
	var D []csd
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	// allocate result
	D = cs_dalloc(noarch.PtrdiffT(n), 0)
	// AT = A'
	AT = cs_transpose(A, 0)
	// get workspace
	xi = cs_malloc(noarch.PtrdiffT((2*int32(n)+1)/8), uint(0)).([]noarch.PtrdiffT)
	if D == nil || AT == nil || xi == nil {
		return (cs_ddone(D, AT, xi, 0))
	}
	Blk = xi
	pstack = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&xi[0])) + (uintptr)(int(n))*unsafe.Sizeof(xi[0]))))[:]
	rcopy = pstack
	p = D[0].p
	r = D[0].r
	ATp = AT[0].p
	top = n
	{
		// first dfs(A) to find finish times (xi)
		for i = 0; i < n; i++ {
			if !(Ap[i] < noarch.PtrdiffT(0/8)) {
				top = cs_dfs(noarch.PtrdiffT(i), A, noarch.PtrdiffT(top), xi, pstack, nil)
			}
		}
	}
	{
		// restore A; unmark all nodes
		for i = 0; i < n; i++ {
			Ap[i] = -noarch.PtrdiffT((Ap[i])) - noarch.PtrdiffT(2/8)
		}
	}
	top = n
	nb = n
	{
		// dfs(A') to find strongly connnected comp
		for k = 0; k < n; k++ {
			// get i in reverse order of finish times
			i = xi[k]
			if ATp[i] < noarch.PtrdiffT(0/8) {
				// skip node i if already ordered
				continue
			}
			// node i is the start of a component in p
			r[func() noarch.PtrdiffT {
				defer func() {
					nb--
				}()
				return nb
			}()] = top
			top = cs_dfs(noarch.PtrdiffT(i), AT, noarch.PtrdiffT(top), p, pstack, nil)
		}
	}
	// first block starts at zero; shift r up
	r[nb] = 0
	for k = nb; k <= n; k++ {
		r[k-nb] = r[k]
	}
	nb = n - nb
	// nb = # of strongly connected components
	D[0].nb = nb
	{
		// sort each block in natural order
		for b = 0; b < nb; b++ {
			for k = r[b]; k < r[b+noarch.PtrdiffT(1/8)]; k++ {
				Blk[p[k]] = b
			}
		}
	}
	for b = 0; b <= nb; b++ {
		rcopy[b] = r[b]
	}
	for i = 0; i < n; i++ {
		p[func() noarch.PtrdiffT {
			tempVar := &rcopy[Blk[i]]
			defer func() {
				*tempVar++
			}()
			return *tempVar
		}()] = i
	}
	return (cs_ddone(D, AT, xi, 1))
}

// cs_schol - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_schol.c:3
// ordering and symbolic analysis for a Cholesky factorization
func cs_schol(order noarch.PtrdiffT, A []cs) []css {
	var n noarch.PtrdiffT
	var c []noarch.PtrdiffT
	var post []noarch.PtrdiffT
	var P []noarch.PtrdiffT
	var C []cs
	var S []css
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	// allocate result S
	S = cs_calloc(1, uint(0)).([]css)
	if S == nil {
		// out of memory
		return nil
	}
	// P = amd(A+A'), or natural
	P = cs_amd(noarch.PtrdiffT(order), A)
	// find inverse permutation
	S[0].pinv = cs_pinv(P, noarch.PtrdiffT(n))
	cs_free(P)
	if bool(order) && S[0].pinv == nil {
		return (cs_sfree(S))
	}
	// C = spones(triu(A(P,P)))
	C = cs_symperm(A, S[0].pinv, 0)
	// find etree of C
	S[0].parent = cs_etree(C, 0)
	// postorder the etree
	post = cs_post(S[0].parent, noarch.PtrdiffT(n))
	// find column counts of chol(C)
	c = cs_counts(C, S[0].parent, post, 0)
	cs_free(post)
	cs_spfree(C)
	// allocate result S->cp
	S[0].cp = cs_malloc(n+noarch.PtrdiffT(1/8), uint(0)).([]noarch.PtrdiffT)
	S[0].lnz = cs_cumsum(S[0].cp, c, noarch.PtrdiffT(n))
	// find column pointers for L
	S[0].unz = S[0].lnz
	cs_free(c)
	return (func() []css {
		if S[0].lnz >= 0 {
			return S
		}
		return cs_sfree(S)
	}())
}

// cs_spsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_spsolve.c:3
// solve Gx=b(:,k), where G is either upper (lo=0) or lower (lo=1) triangular
func cs_spsolve(G []cs, B []cs, k noarch.PtrdiffT, xi []noarch.PtrdiffT, x []float64, pinv []noarch.PtrdiffT, lo noarch.PtrdiffT) noarch.PtrdiffT {
	var j noarch.PtrdiffT
	var J noarch.PtrdiffT
	var p noarch.PtrdiffT
	var q noarch.PtrdiffT
	var px noarch.PtrdiffT
	var top noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Gp []noarch.PtrdiffT
	var Gi []noarch.PtrdiffT
	var Bp []noarch.PtrdiffT
	var Bi []noarch.PtrdiffT
	var Gx []float64
	var Bx []float64
	if !(G != nil && noarch.PtrdiffT(G[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || !(B != nil && noarch.PtrdiffT(B[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || xi == nil || x == nil {
		return noarch.PtrdiffT((-1))
	}
	Gp = G[0].p
	Gi = G[0].i
	Gx = G[0].x
	n = noarch.PtrdiffT(G[0].n)
	Bp = B[0].p
	Bi = B[0].i
	Bx = B[0].x
	// xi[top..n-1]=Reach(B(:,k))
	top = cs_reach(G, B, noarch.PtrdiffT(k), xi, pinv)
	{
		// clear x
		for p = top; p < n; p++ {
			x[xi[p]] = 0
		}
	}
	{
		// scatter B
		for p = Bp[k]; p < Bp[k+noarch.PtrdiffT(1/8)]; p++ {
			x[Bi[p]] = Bx[p]
		}
	}
	for px = top; px < n; px++ {
		// x(j) is nonzero
		j = xi[px]
		// j maps to col J of G
		J = noarch.PtrdiffT(func() int32 {
			if pinv != nil {
				return int32(noarch.PtrdiffT((pinv[j])))
			}
			return int32(noarch.PtrdiffT(j))
		}() / 8)
		if J < noarch.PtrdiffT(0/8) {
			// column J is empty
			continue
		}
		// x(j) /= G(j,j)
		x[j] /= Gx[func() int32 {
			if bool(noarch.PtrdiffT(lo)) {
				return int32(noarch.PtrdiffT((Gp[J])))
			}
			return (int32(Gp[J+noarch.PtrdiffT(1/8)] - noarch.PtrdiffT(1/8)))
		}()]
		// lo: L(j,j) 1st entry
		p = noarch.PtrdiffT(func() int32 {
			if bool(noarch.PtrdiffT(lo)) {
				return (int32(Gp[J] + noarch.PtrdiffT(1/8)))
			}
			return int32(noarch.PtrdiffT((Gp[J])))
		}() / 8)
		// up: U(j,j) last entry
		q = noarch.PtrdiffT(func() int32 {
			if bool(noarch.PtrdiffT(lo)) {
				return int32(noarch.PtrdiffT((Gp[J+noarch.PtrdiffT(1/8)])))
			}
			return (int32(Gp[J+noarch.PtrdiffT(1/8)] - noarch.PtrdiffT(1/8)))
		}() / 8)
		for ; p < q; p++ {
			// x(i) -= G(i,j) * x(j)
			x[Gi[p]] -= Gx[p] * x[j]
		}
	}
	// return top of stack
	return noarch.PtrdiffT((top))
}

// cs_vcount - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_sqr.c:3
// compute nnz(V) = S->lnz, S->pinv, S->leftmost, S->m2 from A and S->parent
func cs_vcount(A []cs, S []css) noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var k noarch.PtrdiffT
	var p noarch.PtrdiffT
	var pa noarch.PtrdiffT
	var n noarch.PtrdiffT = noarch.PtrdiffT(A[0].n)
	var m noarch.PtrdiffT = noarch.PtrdiffT(A[0].m)
	var Ap []noarch.PtrdiffT = A[0].p
	var Ai []noarch.PtrdiffT = A[0].i
	var next []noarch.PtrdiffT
	var head []noarch.PtrdiffT
	var tail []noarch.PtrdiffT
	var nque []noarch.PtrdiffT
	var pinv []noarch.PtrdiffT
	var leftmost []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var parent []noarch.PtrdiffT = S[0].parent
	pinv = cs_malloc(m+n, uint(0)).([]noarch.PtrdiffT)
	// allocate pinv,
	S[0].pinv = pinv
	leftmost = cs_malloc(noarch.PtrdiffT(m), uint(0)).([]noarch.PtrdiffT)
	// and leftmost
	S[0].leftmost = leftmost
	// get workspace
	w = cs_malloc(m+noarch.PtrdiffT(3*int32(n)/8), uint(0)).([]noarch.PtrdiffT)
	if pinv == nil || w == nil || leftmost == nil {
		// pinv and leftmost freed later
		cs_free(w)
		// out of memory
		return noarch.PtrdiffT((0))
	}
	next = w
	head = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:]
	tail = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&(*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0])) + (uintptr)(int(n))*unsafe.Sizeof((*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0]))))[:]
	nque = (*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&(*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0])) + (uintptr)(int(2*int32(n)))*unsafe.Sizeof((*(*[1000000000]noarch.PtrdiffT)(unsafe.Pointer(uintptr(unsafe.Pointer(&w[0])) + (uintptr)(int(m))*unsafe.Sizeof(w[0]))))[:][0]))))[:]
	{
		// queue k is empty
		for k = 0; k < n; k++ {
			head[k] = noarch.PtrdiffT(-1)
		}
	}
	for k = 0; k < n; k++ {
		tail[k] = noarch.PtrdiffT(-1)
	}
	for k = 0; k < n; k++ {
		nque[k] = 0
	}
	for i = 0; i < m; i++ {
		leftmost[i] = noarch.PtrdiffT(-1)
	}
	for k = n - noarch.PtrdiffT(1/8); k >= 0; k-- {
		for p = Ap[k]; p < Ap[k+noarch.PtrdiffT(1/8)]; p++ {
			// leftmost[i] = min(find(A(i,:)))
			leftmost[Ai[p]] = k
		}
	}
	{
		// scan rows in reverse order
		for i = m - noarch.PtrdiffT(1/8); i >= 0; i-- {
			// row i is not yet ordered
			pinv[i] = noarch.PtrdiffT(-1)
			k = leftmost[i]
			if k == noarch.PtrdiffT(int32(-1)/8) {
				// row i is empty
				continue
			}
			if func() noarch.PtrdiffT {
				tempVar := &nque[k]
				defer func() {
					*tempVar++
				}()
				return *tempVar
			}() == noarch.PtrdiffT(0/8) {
				// first row in queue k
				tail[k] = i
			}
			// put i at head of queue k
			next[i] = head[k]
			head[k] = i
		}
	}
	S[0].lnz = 0
	S[0].m2 = m
	{
		// find row permutation and nnz(V)
		for k = 0; k < n; k++ {
			// remove row i from queue k
			i = head[k]
			// count V(k,k) as nonzero
			S[0].lnz++
			if i < noarch.PtrdiffT(0/8) {
				// add a fictitious row
				i = func() noarch.PtrdiffT {
					tempVar := &S[0].m2
					defer func() {
						*tempVar++
					}()
					return *tempVar
				}()
			}
			// associate row i with V(:,k)
			pinv[i] = k
			if func() noarch.PtrdiffT {
				tempVar := &nque[k]
				*tempVar--
				return *tempVar
			}() <= 0 {
				// skip if V(k+1:m,k) is empty
				continue
			}
			// nque [k] is nnz (V(k+1:m,k))
			S[0].lnz += float64(nque[k])
			if (func() noarch.PtrdiffT {
				pa = parent[k]
				return pa
			}()) != noarch.PtrdiffT(int32(-1)/8) {
				if nque[pa] == noarch.PtrdiffT(0/8) {
					// move all rows to parent of k
					tail[pa] = tail[k]
				}
				next[tail[k]] = head[pa]
				head[pa] = next[i]
				nque[pa] += nque[k]
			}
		}
	}
	for i = 0; i < m; i++ {
		if pinv[i] < noarch.PtrdiffT(0/8) {
			pinv[i] = func() noarch.PtrdiffT {
				defer func() {
					k++
				}()
				return k
			}()
		}
	}
	cs_free(w)
	return noarch.PtrdiffT((1))
}

// cs_sqr - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_sqr.c:60
// symbolic ordering and analysis for QR or LU
func cs_sqr(order noarch.PtrdiffT, A []cs, qr noarch.PtrdiffT) []css {
	var n noarch.PtrdiffT
	var k noarch.PtrdiffT
	var ok noarch.PtrdiffT = 1
	var post []noarch.PtrdiffT
	var S []css
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	// allocate result S
	S = cs_calloc(1, uint(0)).([]css)
	if S == nil {
		// out of memory
		return nil
	}
	// fill-reducing ordering
	S[0].q = cs_amd(noarch.PtrdiffT(order), A)
	if bool(order) && S[0].q == nil {
		return (cs_sfree(S))
	}
	if bool(noarch.PtrdiffT(qr)) {
		var C []cs = func() []cs {
			if bool(noarch.PtrdiffT(order)) {
				return cs_permute(A, nil, S[0].q, 0)
			}
			return (A)
		}()
		// QR symbolic analysis
		// etree of C'*C, where C=A(:,q)
		S[0].parent = cs_etree(C, 1)
		post = cs_post(S[0].parent, noarch.PtrdiffT(n))
		// col counts chol(C'*C)
		S[0].cp = cs_counts(C, S[0].parent, post, 1)
		cs_free(post)
		ok = noarch.PtrdiffT(C != nil && S[0].parent != nil && S[0].cp != nil && bool(cs_vcount(C, S)))
		if bool(noarch.PtrdiffT(ok)) {
			S[0].unz = 0
			k = 0
			for k = 0; k < n; k++ {
				S[0].unz += float64(S[0].cp[k])
			}
		}
		if bool(noarch.PtrdiffT(order)) {
			cs_spfree(C)
		}
	} else {
		// for LU factorization only,
		S[0].unz = float64(4*int32(A[0].p[n]) + int32(n))
		// guess nnz(L) and nnz(U)
		S[0].lnz = S[0].unz
	}
	// return result S
	return (func() []css {
		if bool(noarch.PtrdiffT(ok)) {
			return S
		}
		return cs_sfree(S)
	}())
}

// cs_symperm - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_symperm.c:3
// C = A(p,p) where A and C are symmetric the upper part stored; pinv not p
func cs_symperm(A []cs, pinv []noarch.PtrdiffT, values noarch.PtrdiffT) []cs {
	var i noarch.PtrdiffT
	var j noarch.PtrdiffT
	var p noarch.PtrdiffT
	var q noarch.PtrdiffT
	var i2 noarch.PtrdiffT
	var j2 noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var Cx []float64
	var Ax []float64
	var C []cs
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	// alloc result
	C = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(n), noarch.PtrdiffT(Ap[n]), noarch.PtrdiffT(bool(values) && Ax != nil), 0)
	// get workspace
	w = cs_calloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	if C == nil || w == nil {
		// out of memory
		return (cs_done(C, w, nil, 0))
	}
	Cp = C[0].p
	Ci = C[0].i
	Cx = C[0].x
	{
		// count entries in each column of C
		for j = 0; j < n; j++ {
			// column j of A is column j2 of C
			j2 = noarch.PtrdiffT(func() int32 {
				if pinv != nil {
					return int32(noarch.PtrdiffT(pinv[j]))
				}
				return int32(noarch.PtrdiffT(j))
			}() / 8)
			for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
				i = Ai[p]
				if i > j {
					// skip lower triangular part of A
					continue
				}
				// row i of A is row i2 of C
				i2 = noarch.PtrdiffT(func() int32 {
					if pinv != nil {
						return int32(noarch.PtrdiffT(pinv[i]))
					}
					return int32(noarch.PtrdiffT(i))
				}() / 8)
				// column count of C
				w[func() int32 {
					if i2 > j2 {
						return int32(noarch.PtrdiffT((i2)))
					}
					return int32(noarch.PtrdiffT((j2)))
				}()]++
			}
		}
	}
	// compute column pointers of C
	cs_cumsum(Cp, w, noarch.PtrdiffT(n))
	for j = 0; j < n; j++ {
		// column j of A is column j2 of C
		j2 = noarch.PtrdiffT(func() int32 {
			if pinv != nil {
				return int32(noarch.PtrdiffT(pinv[j]))
			}
			return int32(noarch.PtrdiffT(j))
		}() / 8)
		for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			i = Ai[p]
			if i > j {
				// skip lower triangular part of A
				continue
			}
			// row i of A is row i2 of C
			i2 = noarch.PtrdiffT(func() int32 {
				if pinv != nil {
					return int32(noarch.PtrdiffT(pinv[i]))
				}
				return int32(noarch.PtrdiffT(i))
			}() / 8)
			Ci[(func() noarch.PtrdiffT {
				q = func() noarch.PtrdiffT {
					tempVar := &w[func() int32 {
						if i2 > j2 {
							return int32(noarch.PtrdiffT((i2)))
						}
						return int32(noarch.PtrdiffT((j2)))
					}()]
					defer func() {
						*tempVar++
					}()
					return *tempVar
				}()
				return q
			}())] = noarch.PtrdiffT(func() int32 {
				if i2 < j2 {
					return int32(noarch.PtrdiffT((i2)))
				}
				return int32(noarch.PtrdiffT((j2)))
			}() / 8)
			if Cx != nil {
				Cx[q] = Ax[p]
			}
		}
	}
	// success; free workspace, return C
	return (cs_done(C, w, nil, 1))
}

// cs_tdfs - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_tdfs.c:3
// depth-first search and postorder of a tree rooted at node j
func cs_tdfs(j noarch.PtrdiffT, k noarch.PtrdiffT, head []noarch.PtrdiffT, next []noarch.PtrdiffT, post []noarch.PtrdiffT, stack []noarch.PtrdiffT) noarch.PtrdiffT {
	var i noarch.PtrdiffT
	var p noarch.PtrdiffT
	var top noarch.PtrdiffT
	if head == nil || next == nil || post == nil || stack == nil {
		// check inputs
		return noarch.PtrdiffT((-1))
	}
	// place j on the stack
	stack[0] = j
	for top >= 0 {
		// while (stack is not empty)
		// p = top of stack
		p = stack[top]
		// i = youngest child of p
		i = head[p]
		if i == noarch.PtrdiffT(int32(-1)/8) {
			// p has no unordered children left
			top--
			// node p is the kth postordered node
			post[func() noarch.PtrdiffT {
				defer func() {
					k++
				}()
				return k
			}()] = p
		} else {
			// remove i from children of p
			head[p] = next[i]
			// start dfs on child node i
			stack[func() noarch.PtrdiffT {
				top++
				return top
			}()] = i
		}
	}
	return noarch.PtrdiffT((k))
}

// cs_transpose - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_transpose.c:3
// C = A'
func cs_transpose(A []cs, values noarch.PtrdiffT) []cs {
	var p noarch.PtrdiffT
	var q noarch.PtrdiffT
	var j noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var n noarch.PtrdiffT
	var m noarch.PtrdiffT
	var Ap []noarch.PtrdiffT
	var Ai []noarch.PtrdiffT
	var w []noarch.PtrdiffT
	var Cx []float64
	var Ax []float64
	var C []cs
	if !(A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8)) {
		// check inputs
		return nil
	}
	m = noarch.PtrdiffT(A[0].m)
	n = noarch.PtrdiffT(A[0].n)
	Ap = A[0].p
	Ai = A[0].i
	Ax = A[0].x
	// allocate result
	C = cs_spalloc(noarch.PtrdiffT(n), noarch.PtrdiffT(m), noarch.PtrdiffT(Ap[n]), noarch.PtrdiffT(bool(values) && Ax != nil), 0)
	// get workspace
	w = cs_calloc(noarch.PtrdiffT(m), uint(0)).([]noarch.PtrdiffT)
	if C == nil || w == nil {
		// out of memory
		return (cs_done(C, w, nil, 0))
	}
	Cp = C[0].p
	Ci = C[0].i
	Cx = C[0].x
	{
		// row counts
		for p = 0; p < Ap[n]; p++ {
			w[Ai[p]]++
		}
	}
	// row pointers
	cs_cumsum(Cp, w, noarch.PtrdiffT(m))
	for j = 0; j < n; j++ {
		for p = Ap[j]; p < Ap[j+noarch.PtrdiffT(1/8)]; p++ {
			// place A(i,j) as entry C(j,i)
			Ci[(func() noarch.PtrdiffT {
				q = func() noarch.PtrdiffT {
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
	return (cs_done(C, w, nil, 1))
}

// cs_updown - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_updown.c:3
// sparse Cholesky update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1)
func cs_updown(L []cs, sigma noarch.PtrdiffT, C []cs, parent []noarch.PtrdiffT) noarch.PtrdiffT {
	var n noarch.PtrdiffT
	var p noarch.PtrdiffT
	var f noarch.PtrdiffT
	var j noarch.PtrdiffT
	var Lp []noarch.PtrdiffT
	var Li []noarch.PtrdiffT
	var Cp []noarch.PtrdiffT
	var Ci []noarch.PtrdiffT
	var Lx []float64
	var Cx []float64
	var alpha float64
	var beta float64 = 1
	var delta float64
	var gamma float64
	var w1 float64
	var w2 float64
	var w []float64
	var beta2 float64 = 1
	if !(L != nil && noarch.PtrdiffT(L[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || !(C != nil && noarch.PtrdiffT(C[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || parent == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	Lp = L[0].p
	Li = L[0].i
	Lx = L[0].x
	n = noarch.PtrdiffT(L[0].n)
	Cp = C[0].p
	Ci = C[0].i
	Cx = C[0].x
	if (func() noarch.PtrdiffT {
		p = Cp[0]
		return p
	}()) >= Cp[1] {
		// return if C empty
		return noarch.PtrdiffT((1))
	}
	// get workspace
	w = cs_malloc(noarch.PtrdiffT(n), uint(8)).([]float64)
	if w == nil {
		// out of memory
		return noarch.PtrdiffT((0))
	}
	f = Ci[p]
	for ; p < Cp[1]; p++ {
		// f = min (find (C))
		f = noarch.PtrdiffT(func() int32 {
			if f < Ci[p] {
				return int32(noarch.PtrdiffT((f)))
			}
			return int32(noarch.PtrdiffT((Ci[p])))
		}() / 8)
	}
	{
		// clear workspace w
		for j = f; j != noarch.PtrdiffT(int32(-1)/8); j = parent[j] {
			w[j] = 0
		}
	}
	{
		// w = C
		for p = Cp[0]; p < Cp[1]; p++ {
			w[Ci[p]] = Cx[p]
		}
	}
	{
		// walk path f up to root
		for j = f; j != noarch.PtrdiffT(int32(-1)/8); j = parent[j] {
			p = Lp[j]
			// alpha = w(j) / L(j,j)
			alpha = w[j] / Lx[p]
			beta2 = beta*beta + float64(sigma)*alpha*alpha
			if beta2 <= 0 {
				// not positive definite
				break
			}
			beta2 = math.Sqrt(beta2)
			delta = func() float64 {
				if sigma > noarch.PtrdiffT(0/8) {
					return (beta / beta2)
				}
				return (beta2 / beta)
			}()
			gamma = float64(sigma) * alpha / (beta2 * beta)
			Lx[p] = delta*Lx[p] + func() float64 {
				if sigma > noarch.PtrdiffT(0/8) {
					return (gamma * w[j])
				}
				return 0
			}()
			beta = beta2
			for p += 1; p < Lp[j+noarch.PtrdiffT(1/8)]; p++ {
				w1 = w[Li[p]]
				w2 = w1 - alpha*Lx[p]
				w[Li[p]] = w2
				Lx[p] = delta*Lx[p] + gamma*func() float64 {
					if sigma > noarch.PtrdiffT(0/8) {
						return w1
					}
					return w2
				}()
			}
		}
	}
	cs_free(w)
	return noarch.PtrdiffT((beta2 > 0))
}

// cs_usolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_usolve.c:3
// solve Ux=b where x and b are dense.  x=b on input, solution on output.
func cs_usolve(U []cs, x []float64) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Up []noarch.PtrdiffT
	var Ui []noarch.PtrdiffT
	var Ux []float64
	if !(U != nil && noarch.PtrdiffT(U[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || x == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(U[0].n)
	Up = U[0].p
	Ui = U[0].i
	Ux = U[0].x
	for j = n - noarch.PtrdiffT(1/8); j >= 0; j-- {
		x[j] /= Ux[Up[j+noarch.PtrdiffT(1/8)]-noarch.PtrdiffT(1/8)]
		for p = Up[j]; p < Up[j+noarch.PtrdiffT(1/8)]-noarch.PtrdiffT(1/8); p++ {
			x[Ui[p]] -= Ux[p] * x[j]
		}
	}
	return noarch.PtrdiffT((1))
}

// cs_spalloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:3
// allocate a sparse matrix (triplet form or compressed-column form)
func cs_spalloc(m, n, nzmax int, values bool, triplet int) *cs {
	A := new(cs)
	if A == nil {
		// allocate the cs struct
		// out of memory
		return nil
	}
	// define dimensions and nzmax
	A.m = m
	A.n = n
	nzmax = noarch.PtrdiffT(func() int32 {
		if nzmax > noarch.PtrdiffT(int32(1)/8) {
			return int32(noarch.PtrdiffT((nzmax)))
		}
		return int32((1))
	}() / 8)
	A[0].nzmax = nzmax
	// allocate triplet or comp.col
	A[0].nz = noarch.PtrdiffT(func() int {
		if bool(noarch.PtrdiffT(triplet)) {
			return 0
		}
		return -1
	}())
	A[0].p = cs_malloc(noarch.PtrdiffT(func() int32 {
		if bool(noarch.PtrdiffT(triplet)) {
			return int32(noarch.PtrdiffT(nzmax))
		}
		return int32(n + noarch.PtrdiffT(1/8))
	}()/8), uint(0)).([]noarch.PtrdiffT)
	A[0].i = cs_malloc(noarch.PtrdiffT(nzmax), uint(0)).([]noarch.PtrdiffT)
	A[0].x = func() interface{} {
		if bool(noarch.PtrdiffT(values)) {
			return cs_malloc(noarch.PtrdiffT(nzmax), uint(8))
		}
		return nil
	}().([]float64)
	return (func() []cs {
		if A[0].p == nil || A[0].i == nil || bool(values) && A[0].x == nil {
			return cs_spfree(A)
		}
		return A
	}())
}

// cs_sprealloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:18
// change the max # of entries sparse matrix
func cs_sprealloc(A []cs, nzmax noarch.PtrdiffT) noarch.PtrdiffT {
	var ok noarch.PtrdiffT
	var oki noarch.PtrdiffT
	var okj noarch.PtrdiffT = 1
	var okx noarch.PtrdiffT = 1
	if A == nil {
		return noarch.PtrdiffT((0))
	}
	if nzmax <= 0 {
		nzmax = noarch.PtrdiffT(func() int32 {
			if A != nil && noarch.PtrdiffT(A[0].nz) == noarch.PtrdiffT(int32(-1)/8) {
				return int32(noarch.PtrdiffT((A[0].p[noarch.PtrdiffT(A[0].n)])))
			}
			return int32(noarch.PtrdiffT(A[0].nz))
		}() / 8)
	}
	nzmax = noarch.PtrdiffT(func() int32 {
		if nzmax > noarch.PtrdiffT(int32(1)/8) {
			return int32(noarch.PtrdiffT((nzmax)))
		}
		return int32((1))
	}() / 8)
	A[0].i = cs_realloc(A[0].i, noarch.PtrdiffT(nzmax), uint(0), (*[100000000]noarch.PtrdiffT)(unsafe.Pointer(&oki))[:]).([]noarch.PtrdiffT)
	if A != nil && noarch.PtrdiffT(A[0].nz) >= 0 {
		A[0].p = cs_realloc(A[0].p, noarch.PtrdiffT(nzmax), uint(0), (*[100000000]noarch.PtrdiffT)(unsafe.Pointer(&okj))[:]).([]noarch.PtrdiffT)
	}
	if A[0].x != nil {
		A[0].x = cs_realloc(A[0].x, noarch.PtrdiffT(nzmax), uint(8), (*[100000000]noarch.PtrdiffT)(unsafe.Pointer(&okx))[:]).([]float64)
	}
	ok = noarch.PtrdiffT(bool(oki) && bool(okj) && bool(okx))
	if bool(noarch.PtrdiffT(ok)) {
		A[0].nzmax = nzmax
	}
	return noarch.PtrdiffT((ok))
}

// cs_spfree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:33
// free a sparse matrix
func cs_spfree(A *cs) *cs {
	if A == nil {
		// do nothing if A already NULL
		return nil
	}
	cs_free(A[0].p)
	cs_free(A[0].i)
	cs_free(A[0].x)
	// free the cs struct and return NULL
	return (cs_free(A).([]cs))
}

// cs_nfree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:43
// free a numeric factorization
func cs_nfree(N []csn) []csn {
	if N == nil {
		// do nothing if N already NULL
		return nil
	}
	cs_spfree(N[0].L)
	cs_spfree(N[0].U)
	cs_free(N[0].pinv)
	cs_free(N[0].B)
	// free the csn struct and return NULL
	return (cs_free(N).([]csn))
}

// cs_sfree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:54
// free a symbolic factorization
func cs_sfree(S []css) []css {
	if S == nil {
		// do nothing if S already NULL
		return nil
	}
	cs_free(S[0].pinv)
	cs_free(S[0].q)
	cs_free(S[0].parent)
	cs_free(S[0].cp)
	cs_free(S[0].leftmost)
	// free the css struct and return NULL
	return (cs_free(S).([]css))
}

// cs_dalloc - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:66
// allocate a cs_dmperm or cs_scc result
func cs_dalloc(m noarch.PtrdiffT, n noarch.PtrdiffT) []csd {
	var D []csd
	D = cs_calloc(1, uint(0)).([]csd)
	if D == nil {
		return nil
	}
	D[0].p = cs_malloc(noarch.PtrdiffT(m), uint(0)).([]noarch.PtrdiffT)
	D[0].r = cs_malloc(m+noarch.PtrdiffT(6/8), uint(0)).([]noarch.PtrdiffT)
	D[0].q = cs_malloc(noarch.PtrdiffT(n), uint(0)).([]noarch.PtrdiffT)
	D[0].s = cs_malloc(n+noarch.PtrdiffT(6/8), uint(0)).([]noarch.PtrdiffT)
	return (func() []csd {
		if D[0].p == nil || D[0].r == nil || D[0].q == nil || D[0].s == nil {
			return cs_dfree(D)
		}
		return D
	}())
}

// cs_dfree - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:79
// free a cs_dmperm or cs_scc result
func cs_dfree(D []csd) []csd {
	if D == nil {
		// do nothing if D already NULL
		return nil
	}
	cs_free(D[0].p)
	cs_free(D[0].q)
	cs_free(D[0].r)
	cs_free(D[0].s)
	// free the csd struct and return NULL
	return (cs_free(D).([]csd))
}

// cs_done - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:90
// free workspace and return a sparse matrix result
func cs_done(C *cs, w *int, x []float64, ok bool) *cs {
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
func cs_idone(p []noarch.PtrdiffT, C []cs, w interface{}, ok noarch.PtrdiffT) []noarch.PtrdiffT {
	// free temporary matrix
	cs_spfree(C)
	// free workspace
	cs_free(w)
	// return result, or free it
	return (func() []noarch.PtrdiffT {
		if bool(noarch.PtrdiffT(ok)) {
			return p
		}
		return cs_free(p).([]noarch.PtrdiffT)
	}())
}

// cs_ndone - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:106
// free workspace and return a numeric factorization (Cholesky, LU, or QR)
func cs_ndone(N []csn, C []cs, w interface{}, x interface{}, ok noarch.PtrdiffT) []csn {
	// free temporary matrix
	cs_spfree(C)
	// free workspace
	cs_free(w)
	cs_free(x)
	// return result if OK, else free it
	return (func() []csn {
		if bool(noarch.PtrdiffT(ok)) {
			return N
		}
		return cs_nfree(N)
	}())
}

// cs_ddone - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_util.c:115
// free workspace and return a csd result
func cs_ddone(D []csd, C []cs, w interface{}, ok noarch.PtrdiffT) []csd {
	// free temporary matrix
	cs_spfree(C)
	// free workspace
	cs_free(w)
	// return result if OK, else free it
	return (func() []csd {
		if bool(noarch.PtrdiffT(ok)) {
			return D
		}
		return cs_dfree(D)
	}())
}

// cs_utsolve - transpiled function from  $GOPATH/src/github.com/Konstantin8105/sparse/CSparse/Source/cs_utsolve.c:3
// solve U'x=b where x and b are dense.  x=b on input, solution on output.
func cs_utsolve(U []cs, x []float64) noarch.PtrdiffT {
	var p noarch.PtrdiffT
	var j noarch.PtrdiffT
	var n noarch.PtrdiffT
	var Up []noarch.PtrdiffT
	var Ui []noarch.PtrdiffT
	var Ux []float64
	if !(U != nil && noarch.PtrdiffT(U[0].nz) == noarch.PtrdiffT(int32(-1)/8)) || x == nil {
		// check inputs
		return noarch.PtrdiffT((0))
	}
	n = noarch.PtrdiffT(U[0].n)
	Up = U[0].p
	Ui = U[0].i
	Ux = U[0].x
	for j = 0; j < n; j++ {
		for p = Up[j]; p < Up[j+noarch.PtrdiffT(1/8)]-noarch.PtrdiffT(1/8); p++ {
			x[j] -= Ux[p] * x[Ui[p]]
		}
		x[j] /= Ux[Up[j+noarch.PtrdiffT(1/8)]-noarch.PtrdiffT(1/8)]
	}
	return noarch.PtrdiffT((1))
}
