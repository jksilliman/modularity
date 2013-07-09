attach givebound_fib.sage

import sys

N = int(sys.argv[1])

b = (N % 10) + 1
c = int(N/10) + 1

if gcd(b,c) ==1 and (b,c) != (4,3):
	print b, ",", c, "\n"
    	D = b^2+4*c
	if not D.is_square():
		d = QuadraticField(D).gen()
	else :
		d = sqrt(D)
	alpha = (b + d)/2
	beta = (b - d)/2
	if D<=0:
		print 'D was 0'
	else:
		P,n = give_bound_fib(b,c)
		no_pow, pow = rule_out(b,c,n)
		if len(P) == 0:
			print "\t","wow! good level"
		if no_pow:
			print "\t There are no perfect pth powers except perhaps for p=", P, "\n"
		else :
			print "\t All of the perfect powers are", pow, "except perhaps for p=", P, "\n"
