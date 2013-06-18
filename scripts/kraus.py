def check_congruence(f, frey_add_reduction, depth=29, torsion_divisor=2):
    q=f.q_expansion(depth+1)
    g = 0
    
    irrational = True
    if f.base_ring().degree() == 1:
        irrational = False
        
    for n in range(3,depth+1):
        if is_prime(n) and (f.level())%n != 0 and frey_add_reduction%n != 0:
            b = ((n+1)**2 - (q[n])**2).norm()
            for m in xrange(-floor(sqrt(n)*2),floor(sqrt(n)*2)+1):
                if (m-n-1)%torsion_divisor == 0:                         
                    c = (m - q[n]).norm()
                    b *= c
            if irrational:
                b *= n # needed if f irrational?
            if b != 0:
                if g == 0:
                    g = b
                g = gcd(g,b)
            #print b
            #print (n,b.factor())
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
    
q_precision = 42
trace_depth = 29

def analyze_level(level, bad_prime_size=7):
    bad_abvars = []
    
    # Check congruences of ModForms, as in Bennett & Skinner Section 4.
    forms = Newforms(level,2, names='a')
    bad_q_expn = []
    for f in forms:
        r = f.base_ring()
        div = check_congruence(f, trace_depth)
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
        print "Level", Integer(n).factor()
        for (curve, div) in curves:
            print "\tEC", (div.factor(), curve.modular_form().q_expansion(5))
        for (f, div) in abvars:
            print "\tAV", (div.factor(), f.q_expansion(5))
        return False
    return True
   
   
   
   
   
 
 
def find_period_fib(b,c):
    R = Integers(4)
    f0 = 0; f1 = 1
    A1 = 0; A2 = 1
    for i in xrange(20):
        x = b*A2 + c*A1
        A1 = A2
        A2 = x
        if (i>= 1) and R(A1) == R(f0) and R(A2) == R(f1):
            return i+1
            
def find_period_luc(b,c):
    R = Integers(4)
    l0 = 2; l1 = b
    A1 = 2; A2 = b
    for i in xrange(20):
        x = b*A2 + c*A1
        A1 = A2
        A2 = x
        if (i>= 0) and R(A1) == R(l0) and R(A2) == R(l1):
            return i+1
            
def  get_congruences(b,c):
    R = Integers(4)
    f0 = 0; f1 = 1
    l0 = 2; l1 = b
    fib_per = find_period_fib(b,c)
    luc_per = find_period_luc(b,c)
    m = lcm(fib_per, luc_per)
    F = [R(0), R(1)]
    L = [R(2), R(b)]
    for i in xrange(m-2):
        x = b*f1 + c*f0
        f0 = f1
        f1 = x
        y = b*l1 + c*l0
        l0 = l1
        l1 = y
        F.append(R(x))
        L.append(R(y))
    return F,L
    
def find_levels(b,c):
    R = Integers(4)
    F,L = get_congruences(b,c)
    m = len(F)
    D = b**2+4*c
    for n in xrange(m):
        f = F[n]
        l = L[n]
        if R(f) == R(2):
            print "n =", n, "mod", m
            print "\tF =",f,"mod 4"
            print "\tF is not a perfect power"
            continue
        elif GF(2)(f) == GF(2)(0) and R(l) == R(2):
            j = (radical(D)).valuation(2)
            radD = radical(D)/2**j
            N = 2*radD
        elif GF(2)(f) == GF(2)(1) and GF(2)(l) == GF(2)(1) \
        and GF(2)(D) == GF(2)(1):
            if GF(2)(n) == GF(2)(1):
                N = 2**2*radical(D)
            if GF(2)(n) == GF(2)(0):
                N = 2**3*radical(D)
        elif GF(2)(f) == GF(2)(1) and GF(2)(l) == GF(2)(0) \
        and GF(2)(D) == GF(2)(0):
            k = D.valuation(2)
            A = D/2**k
            if k == 2:
                N = 2**5*radical(A)
            if k == 3:
                N = 2**7*radical(A)
            if k == 4:
                if R(-A) == R(f):
                    N = 2**2*radical(A)
                if R(A) == R(f):
                    N = 2**3*radical(A)
            if k == 5:
                N = 2**5*radical(A)
            if k == 6 or k == 7:
                N = 2**3*radical(A)
            if k == 8:
                N = radical(A)
            if k >8:
                N = 2*radical(A)
        print "n =", n, "mod", m
        print "\tF =", f, "mod 4"
        print "\tL =", l, "mod 4"
        print "\tDescent level:", N
   
   
   
   
   
   
   
   
   
   
   
   
# Kraus Method
# (b**2 + 4c)Fn**2 + 4(-c)**n = Ln**2
# c = +- 1, Ln, Fn odd.

def print_unconditional_success(y,n):
    print "Success! Fn = y**p has no solutions for y >=",y, "n =", n, "mod 2"
def print_success(y,n,p_range,try_global):
    if try_global:
        print "Success! Fn = y**p has no solutions for y >=",y, "n =", n, "mod 2", "p in",
        print_range(p_range) 
    else:
        print "Conditional success. If:"
        print "\tFn = y**p has a solutions for y >=",y, "n =", n, "mod 2", "p in",
        print_range(p_range)
        print "then:"
        print"\t y**2p != 1 mod l for one of our allowed l"
def print_range(r):
    print "["+str(r[0])+", " + str(r[-1]) +"]" 
    
def kraus(b,c,p_range,n_parity=1,max_k=1000, try_global=True, max_global_attempt=50):
    print "b =",b
    print "c =",c
    print "p ",
    print_range(p_range)
    print "n =",n_parity, "mod 2"
    
    
    d = b**2 + 4*c
    # Roots of x**2 - bx - c
    #omega = (b + sqrt(d))/2
    #tau = (b + sqrt(d))/2
    
    # Fix p >= 7 prime. Fix parity of n. We attempt to prove that Fn = y**p has no solutions.
    # Frey Curve E: Y**2 = X**3 + Hn X**2 + (-c)**nX, 
    # where Hn = +- Ln such that Hn = -(-c)**n mod 4.
    # N(E) = 2**alpha * rad(d)
    n = n_parity
    
    if ((-c)**n+1)%4 == 0:
        alpha = 2
    elif ((-c)**n-1)%4 == 0:
        alpha = 3
    level = 2**alpha * rad(d)
    add_red = 2
    print "Level: ", level
    (curves, abvars) = analyze_level(level,add_red)
    if len(abvars) > 0:
        print "WARNING: Unruly Abelian Varieties Found"
        print abvars 
    if len(curves) == 0:
        print_unconditional_success(2,n)
        return
    overall_success = True
    for (F, bad_p) in curves:
        print F        
        for p in p_range:
            if not p.is_prime():
                continue
            
            success = False
            attempt_counter = 0
            
            k = 0
            while k < max_k: 
                k += 1
                if try_global and attempt_counter > max_global_attempt:
                    try_global = False
                    print "WARNING: We stop trying global solution at p =",p                   
                    k = 0
                    continue
                if (k%1000)==0: 
                    print "\tk = ", k
                    
                l = 2*p*k + 1
        
                
                
                if not l.is_prime():
                    continue
        
                Fl = GF(l)
                g = Fl.multiplicative_generator()**(2*p)     
                # Does l divide the conductor? If so, skip. THEORY-PROBLEM: This works for y >= |b| + |c|.
                if (2*d)%l == 0:
                    continue
                # We can't check unless b**2 + 4c a qr.
                if Fl(d).is_square():
                    omega = (Fl(b) + Fl(d).sqrt())/Fl(2)
                    if omega**(2*k) == Fl(1):
                        continue
                else:
                    continue
                
                
                #print "\tPossible l:", l
                # We conclude that either l does not divide y, or y < |b| + |c|. We assume l does not divide y, hence not the conductor of Frey Curve.
                
                attempt_counter += 1
        
                # A = []
                Es = []
                
                z = g
                
                if try_global: # If we want to allow y**2p = 1 mod l 
                    top = k+1
                else:
                    top = k
                for i in range(1,top): #z is of order k, we go up to z**k=1
                    m = d*z + 4*(-c)**n
                    if m.is_square():
                        # A.append(z)
                        d_z = sqrt(m)
                        Es.append( EllipticCurve(GF(l), [0,d_z,0,(-c)**n,0]) )
                    z *= g
                    
                
                success = True
                for Ez in Es:
                    if Ez.trace_of_frobenius() == F.change_ring(Fl).trace_of_frobenius():
                        # This l was no good. Move along.     
                        success = False
                        break
                if success:
                    #print "\tSuccessful l:", l
                    break
            if not success:
                print "No success at p =", p, ". Improve Kraus' method, or find a soln to Fn = y**p. Curve:", F  
                overall_success = False      
    if overall_success:
        print_success(abs(b)+abs(c),n,p_range, try_global)
    else:
      print "There was failure at some p."
