from itertools import *

## Definitions
B = 1000
pr = list(primes(B))

print("Computing hilbert for primes...")
prime_hilbert = [ [ hilbert_symbol(p,q,p) == -1 for q in pr ] for p in pr ]

def squarefrees( n ):
	return filter( is_squarefree, range(n) )

def products_less_than( L ):
	M = [(1,[])]
	for i, y in enumerate(L):
		n = len(M)
		for j in range(n):
			x,f = y*M[j][0], M[j][1] + [i]
			if x < B:
				M.append((x,f))
	return dict(M)

def sections( A, B ):
	i, j = 0, 0
	n, m = len(A), len(B)
	L, M, R = [], [], []
	while i < n and j < m:
		if A[i] < B[j]:
			L.append( A[i] )
			i += 1
		elif A[i] > B[j]:
			R.append( B[j] )
			j += 1
		else:
			M.append( A[i] )
			i += 1
			j += 1
	while i < n:
		L.append( A[i] )
		i += 1
	while j < m:
		R.append( B[j] )
		j += 1
	return M, L, R

def global_hilbert( fa, fb ):
	_, a = fa
	_, b = fb
	M, A, B = sections( a, b )
	for i in A:
		if i == 0:
			continue
		s = 0
		for j in b:
			s += prime_hilbert[i][j]
		if s % 2 == 1:
			return 0
	for i in B:
		if i == 0:
			continue
		s = 0
		for j in a:
			s += prime_hilbert[i][j]
		if s % 2 == 1:
			return 0
	for i in M:
		if i == 0:
			continue
		s = prime_hilbert[i][i]
		for j in chain(a,b):
			s += prime_hilbert[i][j]
		if s % 2 == 1:
			return 0
	return 1

def global_hilbert_check( fa, fb ):
	a, _ = fa
	b, _ = fb
	for p in pr:
		if hilbert_symbol( a, b, p ) == -1:
			return 0
	return 1

## Computation
print("Computing squarefree integers...")
sf = products_less_than(pr)

print("Computing hilbert for pairs...")
hilbert_pairs = [ [ global_hilbert(x,y) for y in sf.items()] for x in sf.items() ] 
adj = [ [ j for j,y in enumerate(sf.keys()) if hilbert_pairs[i][j] ] for i,x in enumerate(sf.keys()) ]

print("Computing hilbert triples...")
c = 0
for x in range(len(sf.keys())):
	for y in adj[x]:
		if y > x:
			break
		for z in adj[x]:
			if z > y:
				break
			if hilbert_pairs[y][z]:
				#print(x,y,z)
				c += 1
print(c)