# The mod n period
def fib_per(b,c,n):
    R = Integers(n)
    f0 = R(0); f1 = R(1)
    A1 = R(0); A2 = R(1)
    i = 0
    while True:
	      i += 1
        x = R(b*A2 + c*A1)
        A1 = A2
        A2 = x
        if A1 == f0 and A2 == f1:
	        return i
    

def fib_recur(b,c,u0,u1,n):
    if n == 0:
        return u0
    elif  n ==1:
        return u1
    else :
        a0 = u0
        a1 = u1
        for i in xrange(n-1):
            x = b*a1 + c*a0
            a0 = a1
            a1 = x
        return x

# Finds primes q = 1 mod p such that FibPeriod(q) | M.
# The second condition lets us use q without increasing M.
def find_primes(p,b,c,M,B):
	L = []
	for q in prime_range(B):
		if (q % p == 1) and (M % fib_per(b,c,q) == 0) :
			L.append(q)
	return L

# Is "a" a pth power mod q?
def is_p_power(a,q,p):
	return GF(q)(a)^((q-1)/p) in (0,1)

def find_cong(p,b,c,M, cong=[0], Modulus=1,B=10000):
	L = find_primes(p,b,c,M,B)
	for ell in L:
		modu = fib_per(b,c,ell)
		
		F = GF(ell)
		u0 = 0; u1 = 1
		a0 = F(u0); a1 = F(u1)
		
		# Note that 0,1 are pth powers for all p
		# This is a list of all congruences we find from working mod ell
		cong1 = [0,1]
		
		# We look at our sequence mod ell
		for m in xrange(2, modu):
			x = F(b*a1 + c*a0)
  		a0 = a1
  		a1 = x
  		
			if is_p_power(x,ell,p):
				cong1.append(m)
				
	  # We "intersect" the two "sets" of congruences.
	  # CRT fails for numbers which are "in" one, but not the other.
		congnew = set([])
		for i in cong1:
			for j in cong:
				try:
					cc = crt(i,j,modu,Modulus)
				except ValueError:
					continue
				congnew.add(cc)					
		cong = list(congnew)		
		Modulus = lcm(modu, Modulus)
  
  # Note Modulus <= M, as FibPeriod(ell) divides M by choice.
  # If Modulus did not reach M, we must increase the bound B on primes ell.
	if Modulus != M:
		return find_cong(p,b,c,M, cong, Modulus=Modulus, B=2*B)
	return cong, Modulus

def is_pth_power(n, p):
	return bool(n**(1/p) == int(n**(1/p)))


# Main Method:
# b,c - recurrence relation
# p - the prime power to eliminate
def check_prime(b,c,p):
  # TODO: This bound on the index should be derived from the bound on p, Thm 2.
	Bound = p^(10*p)
	
	# We make a list of the small powers in our sequence, in order to ignore their residue classes later
	powers = [0,1]
	for i in xrange(2,20):
		if is_pth_power(fib_recur(b,c,0,1,i),p):
			powers.append(i)
		
	M = 2^5*3^3*5^2*7*11	
	cong, Modulus = [0], 1
	while True:
		cong, Modulus = find_cong(p,b,c,M,cong,Modulus)
		cong.sort()
		if (len(cong) > len(powers) and cong[len(powers)] > Bound)
		  break
		elsif M > Bound:
			break
		M *= next_prime(M.prime_factors()[-1])
	return powers






