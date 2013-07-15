#This is code to sieve for perfect pth powers in Lucas sequences


#This function determines if the period of the sequence (b,c) mod q divides M 
#using the fact that the period mod q is the order of the matrix F = [[0,1],[c,b]]
#in GL_2(F_q).  As the order divides (q-1) -- because b^2+4c is a square mod q -- the 
#order of F divides M iff it divides M mod (q-1)
def per_divide(b,c,q,M):
	F = matrix(GF(q), [[0,1],[c,b]])
	I = matrix(GF(q), [[1,0],[0,1]])
	return (F^(M % (q - 1)) == I)

#This function returns some prime q>ell which is 1 mod p AND 
#AND b^2+4c is a square mod ell, AND period(q) divides M2 but 
#period(q) does NOT divide M1.
def next_good_prime(p,b,c,M1,M2,ell):
	found = False
	while not found:
		ell  = next_prime(ell)
		if ell % p == 1:	#only computes further 1/pth of the time
			if legendre_symbol(b^2 + 4*c, ell) == 1: # makes the order divide (ell-1) similar to Siksek
				if per_divide(b, c, ell, M2):
					if not per_divide(b, c, ell, M1):
						found  = True
						return ell 

#Naively calculates the nth term in an arbitrary recurrence sequence
#defined by (b,c) with starting conditions u0 and u1 (note, could be 
#optimized using matrices with a double and add algorithm)
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

#This function checks if a is a pth power mod q
def is_p_power_mod(a,q,p):
	return GF(q)(a)^((q-1)/p) in (0,1)


#This function finds the list of possible congruences for index n mod M,
#for the first time (optimized for starting)
def find_cong_start(p,b,c,M):
	done = False
	cong = [0]
	Modulus = 1
	ell = 1
	tries = 0
	while not done:
		ell = next_good_prime(p,b,c,1,M,ell)  	#take one prime at a time
		print "Working on", ell
		F = GF(ell)
		u0 = 0; u1 = 1
		a0 = F(u0); a1 = F(u1)
		#cong1 keeps track of congruences mod period(ell) for the sequence mod ell
		cong1 = [0,1]		#u_0=0 and u_1 = 1 will be pth powers mod any prime for any p
		found_period = False
		n = 1
		while not found_period:
			n += 1		#n is the index of the a1 variable
			x = F(b*a1 + c*a0)
            		a0 = F(a1)
            		a1 = F(x)
			if a0 == F(u0) and a1 == F(u1):		#found the period
				modu = n - 1		#the period is the index of the a0 variable not the a1 variable
				del cong1[-1]		#we have already added a second a0=0 to our congruence list, so we delete
				found_period = True	
			elif is_p_power_mod(x,ell,p):
				cong1.append(n)		#if haven't found the period, add a possible congruence to the list

		#use the CRT to intersect the congruences of cong1 and the 
		congnew = set([])
		for i in cong1:
			for j in cong:
				try:
					cc = crt(i,j,modu,Modulus)
				except ValueError:
					continue
				congnew.add(cc)
		congnew = list(congnew)
		Modulus = lcm(modu, Modulus)
		#print Modulus
		
		print Modulus, M, len(congnew)

		if Modulus == M:
			if len(cong) == len(congnew):
				tries += 1
				if tries == 2:		#try twice to get less congruences when Modulus = M
					done = True
		cong = congnew
	cong.sort()
	return cong, Modulus

#This is the main function that finds the list of possible congruences for index n mod M2,
#given a list of congruences in cong already mod Modulus (optimized for later steps)
#note that M1 = Modulus
def find_cong(p,b,c,M1,M2, cong, Modulus):
	done = False
	ell = 1
	tries = 0
	CongNew = []		#makes a new list from cong that is now mod M2 instead of M1
	for k in xrange(M2/M1):
		for i in cong:
			CongNew.append(k*M1+i)
	cong = set(CongNew)
	while not done:
		ell = next_good_prime(p,b,c,M1,M2,ell)  	#take one prime at a time
		print "Working on", ell
		#modu = fib_per(b,c,ell)
		F = GF(ell)
		u0 = F(0); u1 = F(1)
		a0 = u0; a1 = u1
		bf, cf = F(b), F(c)
		#cong1 keeps track of congruences mod period(ell) for the sequence mod ell
		cong1 = [0,1]		#u_0=0 and u_1 = 1 will be pth powers mod any prime for any p
		found_period = False
		n = 1
		ppowers = set([])		#make the set (F_ell)^p
		for i in F:
			ppowers.add(i^p)
		while not found_period:
			n += 1		#n is the index of the a1 variable
            		a0, a1 = a1, bf*a1 + cf*a0
			if a0 == u0 and a1 == u1:		#found the period
				modu = n - 1		#the period is the index of the a0 variable not the a1 variable
				del cong1[-1]		#but we have already added a second 0 to our congruence list, so we delete
				found_period = True	
			elif a1 in ppowers:
				cong1.append(n)		#if haven't found the period, add a possible congruence to the list
		cong1 = set(cong1)

		killed_something = False		#keeps track of when cong1 can rule out a congruence in cong
		
		#CRT by hand to gain speed
		for i in list(cong):
			if not (i % modu in cong1):		#congruence in cong is inconsistent with any in cong1
				cong.remove(i)			#remove that congruence
				killed_something = True
		
		print M2, cong

		if not killed_something:		
			tries += 1
			if tries == 2:			#try twice to rule out congruences
				done = True
	cong = list(cong)
	cong.sort()
	return cong, M2

#check for when a is a perfect pth power
def is_pth_power(a, p):
	return bool(a^(1/p) == int(a^(1/p)))

#This is the main function that rules out pth powers given a bound Bound for n in terms of p
def check_prime(b,c,p,Bound):
	powers = [0,1]
	for i in xrange(2,20):			#find all perfect powers for n small, this is our "conjecturally complete set"
		if is_pth_power(fib_recur(b,c,0,1,i),p):
			powers.append(i)
	M2 = 2^5*3^3*5^2*7*11		#starting modulus
	if not (M2 % p == 0):		#make sure p divides the starting modulus
		M2 *= p
	qq = 1
	cong, Modulus = find_cong_start(p,b,c,M2)		#find the first set of congruences for our starting modulus
	Done  = False

	#This section check for a contradiction
	if len(cong) > len(powers) and cong[len(powers)] > Bound:	#contradiction if smallest congruence that grows with Modulus > Bound
	elif Modulus > Bound:		#if there are no such congruences, then if Modulus > Bound
		Done = True

	j = 1			#variable that keeps track of how many times we had to add prime powers to Modulus
	while not Done :
		j += 1
		#print "adding the", j, "prime power"
		M1 = M2
		M2 = lcm(M2, qq)		#step up prime powers from initially qq = 1
		qq = next_prime_power(qq)
		if M1 == M2:			#such that the new modulus is actually distinct
			continue
		cong, Modulus = find_cong(p, b, c,M1, M2, cong, Modulus)	#using this new modulus, find the new cong

		#check for a contraduction again
		if len(cong) > len(powers) and cong[len(powers)] > Bound:
			Done = True
		elif Modulus > Bound:
			Done = True

	#print "Added", j, "new factors"
	return powers, cong, Modulus






