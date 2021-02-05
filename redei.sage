from itertools import *

def get_sol( a, b ):
	try:
		M = DiagonalQuadraticForm(QQ, [a, b, -1] )
		(x,y,z) = M.solve()
		return (abs(x),abs(y),abs(z))
	except:
		return None

def valid_norm( a, b ):
	if b == 1 or b == -1:
		return False
	return get_sol( a, b ) != None

def gen_norms( a ):
	A = list({ x*x-a*y*y for x in range( 1, 15 ) for y in range( 1, 15 ) })
	B = [ squarefree_part(b) for b in A ]
	# C = [ b for b in B if valid_norm( a, b ) ]
	return B


def get_candidate():
	#for a in primes(100):
	while True:
		a = squarefree_part( randint(2,1000) )
		da = disc(a)
		for b in gen_norms( a ):
			db = disc(b)
			for c in gen_norms( a ):
				if gcd([da,db,disc(c)]) == 1:
					if get_sol( b, c ) != None:
						return (a,b,c)

def str_solution( a, b ):
	(x,y,z) = get_sol( a, b )
	return "{a}*{x}^2 + {b}*{y}^2 = {z}^2".format( a=a, b=b, x=x, y=y, z=z )

def disc( a ):
	if a % 4 == 1:
		return a
	return 4*a

def triple_solutions( a, b, c ):
	print( a, b, c )
	print( disc(a), disc(b), disc(c) )
	print( str_solution( a, b ) )
	print( str_solution( b, c ) )
	print( str_solution( a, c ) )

def prime_divisors( a ):
	return [ x for x, _ in a.factor() ]

def redei( a, b, c, sol = None, verbose = True ):
	# 0 is error
	if verbose:
		print(a,b,c)

	# assumptions
	a = squarefree_part(a)
	b = squarefree_part(b)
	c = squarefree_part(c)
	if 1 in [a,b,c]:
		return 1
	if a == b:
		a, b, c = b, c, a
	if disc(b) % 2 == 0:
		a, b = b, a
	if gcd([disc(a),disc(b),disc(c)]) != 1:
		return 0

	if verbose:
		print(a,b,c)

	# compute twist
	(x,y,z) = (0,0,0)
	if sol != None:
		(x,y,z) = sol
		if a*x*x+b*y*y != z*z:
			raise Exception( "Not a solution" )
	else:
		try:
			(x,y,z) = get_sol( a, b )
		except:
			return 0
	if verbose:
		print(x,y,z)
	P.<X> = QQ[]
	Ka.<sa> = QQ.extension( X^2-a )
	# beta = (z/y)+(x/y)*sa
	beta = z+x*sa
	da = disc(a)
	if verbose:
		print( "beta:", beta )
	K.<sab> = QQ.extension( X^2 - a*b )
	t = 1
	if disc(b) % 2 == 1:
		if da % 2 == 1 or disc(b) % 8 == 1:
			if verbose:
				print("Case 1 & 2:",disc(b))
			for t in [-2,-1,1,2]:
				Ft.<sbetat> = K.extension( X^4 - t*beta.trace()*X^2 + t*t*beta.norm() )
				norm = Ft.relative_discriminant().norm() 
				if verbose:
					print("norm:",norm)
				if norm % 2 == 1:
					break
			else:
				return 0
		else:
			if verbose: # BAD CASE
				print("Case 3")
			p2 = Ka.prime_above(2)
			I = p2**3
			#d = 1
			#for n in range(1,9):
			#	if I.reduce( y**n-1 ) == 0:
			#		d = y**n
			#		break
			#else:
			#	return 1337
			#betad = beta*d
			#print( betad )
			tau = (1+sa)^2/2
			for t, te in [(-2,-tau),(-1,-1),(1,1),(2,tau)]:
				if verbose:
					print( I.reduce( te*beta-1 ) )
				if I.reduce( te*beta -1 ) == 0:
					break
			else:
				return 0
		if verbose:
			print( "Twisting parameter:", t )
		beta = beta*t
	elif verbose:
		print("Case 4")

	if verbose:
		F.<sbeta> = K.extension( X^4 - beta.trace()*X^2 + beta.norm() )
		print( "disc: ", F.relative_discriminant() )

	# local reduction
	Kasbeta.<sbeta> = Ka.extension(X^2-beta)
	global_symbol = 1
	for p in prime_divisors(c):
		if verbose:
			print("Checking prime", p)
		Kap = None
		if p == 2:
			for Kap in Ka.primes_above(2):
				Kasbetap = Kasbeta.prime_above( Kap )
				e = Kasbetap.relative_ramification_index()
				if verbose:
					print("ramification index:",e)
				if e == 1:
					break
			else:
				return 0
		else:
			for Kap in Ka.primes_above( p ):
				if Ka.valuation( Kap )( beta ) % 2 == 0:
					break
			else:
				return 0
		pi = Ka.uniformizer( Kap )
		hs = Ka.hilbert_symbol( beta, pi, Kap )
		if verbose:
			print( "hilbert:", hs )
		global_symbol *= hs

	# prime at infinity
	if c < 0 and a > 0 and b > 0:
		if verbose:
			print("Checking prime infinity")
		hs = Ka.hilbert_symbol( beta, -1, Ka.real_embeddings()[0] )
		if verbose:
			print( "hilbert:",hs)
		global_symbol *= hs

	return global_symbol