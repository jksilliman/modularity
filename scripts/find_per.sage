def fib_per(b,c, n):
    R = Integers(n)
    f0 = 0; f1 = 1
    A1 = R(0); A2 = R(1)
    found = False
    i = 0
    while not found:
	i += 1
        x = R(b*A2 + c*A1)
        A1 = R(A2)
        A2 = R(x)
        if i>=1 and R(A1) == R(f0) and R(A2) == R(f1):
	    found = True
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

def fib_recur_mod(b,c,u0,u1,p,n):
    F = GF(p)
    if n == 0:
        return F(u0)
    elif  n ==1:
        return F(u1)
    else :
        a0 = F(u0)
        a1 = F(u1)
        for i in xrange(n-1):
            x = F(b*a1 + c*a0)
            a0 = F(a1)
            a1 = F(x)
        return x

def find_primes(p,b,c,M,B):
	L = []
	for q in prime_range(B):
		#print q, fib_per(b,c,q)
		if (q % p == 1) and (M % fib_per(b,c,q) == 0) :
			L.append(q)
	return L

def is_p_power(a,q,p):
	return GF(q)(a)^((q-1)/p) in (0,1)

def find_cong(p,b,c,M, cong=[0], Modulus=1, B=10000):
	L = find_primes(p,b,c,M,B)
	print L
	for ell in L:
		print "Working on", ell
		modu = fib_per(b,c,ell)
		F = GF(ell)
		u0 = 0; u1 = 1
		a0 = F(u0); a1 = F(u1)
		cong1 = [0,1]
		for m in xrange(2, fib_per(b,c,ell)):
			x = F(b*a1 + c*a0)
            		a0 = F(a1)
            		a1 = F(x)
			if is_p_power(x,ell,p):
				cong1.append(m)
		print "About to CRT"
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
		print Modulus
	if Modulus != M:
		return find_cong(p,b,c,M, cong, Modulus, 2*B)
	cong.sort()
	return cong, Modulus

def is_pth_power(n, p):
	return bool(n**(1/p) == int(n**(1/p)))

def check_prime(b,c,p):
	Bound = p^(10*p)
	powers = [0,1]
	for i in xrange(2,20):
		if is_pth_power(fib_recur(b,c,0,1,i),p):
			powers.append(i)
	M = 2^5*3^3*5^2*7*11
	qq = M.prime_factors()[-1]
	cong, Modulus = find_cong(p,b,c,M)
	Done  = False
	if len(cong) > len(powers) and cong[len(powers)] > Bound:
		Done = True
	elif Modulus > Bound:
		Done = True
	while not Done :
		M *= next_prime(qq)
		qq = M.prime_factors()[-1]
		cong, Modulus = find_cong(p,b,c,M,cong, Modulus)

		if len(cong) > len(powers) and cong[len(powers)] > Bound:
			Done = True
		elif Modulus > Bound:
			Done = True
	return powers






