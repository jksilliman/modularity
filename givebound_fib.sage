def check_congruence(f, add_reduction, depth=29, torsion_divisor=2, mult_only = False, hasse_only = False):
    q=f.q_expansion(depth+1)
    g = 0
    
    irrational = True
    if f.base_ring().degree() == 1:
        irrational = False
        
    n = 2
    for n in range(2,depth+1):
        bad = False
        for p in add_reduction:
            if n==p:
                bad = True
        if bad:
            continue
        
        if is_prime(n) and (f.level())%n != 0:
            if mult_only == True:
            	b = ((n+1)**2 - (q[n])**2).norm()
            if hasse_only == True:
                b = 1
                for m in xrange(-floor(sqrt(n)*2),floor(sqrt(n)*2)+1):
                    if torsion_divisor == 1 or (m-n-1)%torsion_divisor == 0:                         
                        c = (m - q[n]).norm()
                        b *= c
            else:
                b = ((n+1)**2 - (q[n])**2).norm()
                for m in xrange(-floor(sqrt(n)*2),floor(sqrt(n)*2)+1):
                    if torsion_divisor == 1 or (m-n-1)%torsion_divisor == 0:                         
                        c = (m - q[n]).norm()
                        b *= c
            if irrational:
                b *= n # needed if f irrational?
            if b != 0:
                if g == 0:
                    g = b
                g = gcd(g,b)
    return g
 
 
 
def dim_ratnl_newform(level):
    n = 0
    forms = Newforms(level,2, names='a')
    for f in forms:
        r = f.base_ring()
        if r.degree() == 1:
            n += 1
    return n
    
def isogeny_classes_level(level):
    curves = []
    for n in range(1, dim_ratnl_newform(level)+1):
        letter = chr(96 + n)
        curves.append(EllipticCurve(str(level) + letter))
    return curves
def largest_prime(num):
    if abs(num) > 1:
        return factor(num)[-1][0]
    return num
def is_rad(num):
    for n, exp in num.factor():
        if exp > 1:
            return false
    return true

def try_factor(num):
  if abs(num) < 2:
    return num
  else:
    return factor(num)

def all_factors(num):
    f = num.factor()
    fs = []
    pw = [0] * len(f)
    pos = 0
    while pos < len(f):
        if pw[pos] < f[pos][1]:
            pw[pos] += 1
            pos = 0
            factor = 1
            for i in range(0, len(f)):
                factor *= pow(f[i][0],pw[i])
            fs.append(factor)
        else: 
            pw[pos] = 0
            pos += 1
    return fs
def rad(num):
    res = 1
    for n, exp in num.factor():
        res *= n
    return res

def square_factors(num):
    squared = []
    for (f,e) in num.factor():
        if e >= 2:
            squared.append(f)
    return squared
    
q_precision = 42
trace_depth = 29

def analyze_level(level, bad_prime_size=7):
    bad_abvars = []
    
    add_reduction = square_factors(level)
    
    # Check congruences of ModForms, as in Bennett & Skinner Section 4.
    forms = Newforms(level,2, names='a')
    bad_q_expn = []
    for f in forms:
        #print f
        r = f.base_ring()
        div = check_congruence(f, add_reduction,trace_depth)
        largest_p = largest_prime(div)
        if largest_p >= bad_prime_size or abs(largest_p) < 2:
            if r.degree() == 1:
                bad_q_expn.append((f.q_expansion(q_precision), div))
            else:
                bad_abvars.append([f, div])  
    # Fetch bad elliptic curves
    bad_curves = []
    curves = isogeny_classes_level(level)
    for c in curves:
        q = c.modular_form().q_expansion(q_precision)
        for (bad_q, div) in bad_q_expn:
            if q == bad_q:
                bad_curves.append([c, div])
    return (bad_curves, bad_abvars)

def describe_level(n):
    (curves, abvars) = analyze_level(n)
    if len(curves) + len(abvars) > 0:
        return False
    return True

def upbound(N):
	B = 1
	for p,e in list(N.factor()):
		B *= p^(e-1)*(p-1)*(p+1)
	return B


def is_com_divisible(b,c,u0,u1,N):
	R = Integers(N)
	v0 = 2*u1-b*u0
	v1 = b*u1+2*c*u0
	if R(v0)==0 or R(v1)==0:
		return True
	a0 = R(v0); a1 = R(v1)
	divi = False
	per = False
	i = 1
	while not divi and not per and i < upbound(N) + 2:
		i += 1
		x = R(b*a1 + c*a0)
		a0 = a1
		a1 = x
		if a1 == 0:
			divi = True
			return True
		if a1 == R(v1) and a0 == R(v0):
			per = True
			return False
	return False

def is_fib_divisible(b,c,u0,u1,N):
	R = Integers(N)
	a0 = R(u0); a1 = R(u1)
	divi = False
	per = False
	i = 1
	while not divi and not per and i < upbound(N) + 2:
		i += 1
		x = R(b*a1 + c*a0)
		a0 = a1
		a1 = x
		if a1 == 0:
			divi = True
			return True
		if a1 == R(u1) and a0 == R(u0):
			per = True
			return False
	return False

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

def give_bound_fib(b,c):
	if gcd(b,c) != 1:
		return 'not relatively prime'
	Bad_primes = []	
	D = b^2+4*c
	print "D is", D.factor()
	if D%4 == 1:
		if c%2 != 0:
			if not is_com_divisible(b,c,0,1,2) :
				if Integers(4)(c) == Integers(4)(-1):
					poss = [2^3*rad(D*c)]
				else :
					poss = [2^2*rad(D*c), 2^3*rad(D*c)]
			else :
				if Integers(4)(c) == Integers(4)(-1):
					poss = [2^3*rad(D*c)]
				else :
					poss = [2^2*rad(D*c), 2^3*rad(D*c)]
				poss.append(2*rad(c*D))
		if c%2 ==0:
			poss = [rad(c*D)]
	else :
		k = D.valuation(2)
		A = D/2^k
		if k == 2:
			poss = [2^5*rad(c*A)]
		if k == 3:
			poss = [2^7*rad(c*A)]
		if k == 4:
			if Integers(4)(A) == Integers(4)(-1):
				poss = [2^2*rad(c*A)]
			if Integers(4)(A) == Integers(4)(1):
				poss = [2^3*rad(c*A)]
		if k == 5:
			poss = [2^2*rad(c*A)]
		if k == 6 or k == 7:
			poss = [2^3*rad(c*A)]
		if k == 8:
			poss = [rad(c*A)]
		if k > 8:
			poss = [2*rad(c*A)]
		if is_fib_divisible(b,c,0,1,2) and not (k > 8) :
			poss.append(2*rad(c*A))
	print "poss is", poss	
	Bound = 0
	for lev in poss:
		#print "Hey, I'm checking a level called", lev, ":)"
		add_reduction = square_factors(lev)
		#print add_reduction
		forms = Newforms(lev,2,names='a')
		#print "I've found the forms!"
		num_rational = 0
		for f in forms:
			#print "Hey, I'm checking a form"
			R = f.base_ring()
			d = R.degree()
			if d > 1:
				div = check_congruence(f, add_reduction)
				#print "yo, I'm at the point where I find the gcd", div
				div = ZZ(div)
				Bad_primes.extend(div.prime_factors())
			if d == 1:
				num_rational += 1
		if num_rational != 0:
			for p in primes(2,17+1):
				Bad_primes.append(p)
		P = lev.prime_factors()[-1]
		possi_bound = max(30, P)
		Bound = max(Bound, possi_bound)
	#print 'yo, Bad_primes is', Bad_primes
	return list(set(Bad_primes)), Bound

def rule_out(b,c,n):
	no_pow = True
	pow = []
	for i in xrange(2,n+2):
		if fib_recur(b,c,0,1,i).is_perfect_power() and fib_recur(b,c,0,1,i) != 0 and fib_recur(b,c,0,1,i) != 1 :
			no_pow = False
			Un = fib_recur(b,c,0,1,i)
			p = 0
			for q,e in list(Un.factor()):
				p = gcd(p,e)
			pow.append((p,i,Un))
	return no_pow, pow



