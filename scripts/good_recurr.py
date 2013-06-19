c = 1
for b in range(1,100000,2):
  d = b^2+4*c
  if not is_prime(d):
    continue
  good = True
  F = 1
  Fprev = 0
  L = b
  Lprev = 2 
  n = 1
  while n <= 12:
    temp = Fprev
    Fprev = F
    F = temp + b*F

    temp = Lprev
    Lprev = L
    L = temp + b*L

    n += 1
    #print(n, F)
    if n==9:
      if Integer(F/(b**2+1)).is_perfect_power():
        good = False
        print "b =", b, "F/b =", factor(F/b) 
    if n==12:
      if Integer(L/2).is_perfect_power():
        good = False
        print "b =", b, "L/2 =", factor(L/2)
      if Integer(F).is_perfect_power():
        print "For b = ", b, "F(12) =", factor(F), " is the only perfect power of even index"
    #if n==d^2:
    #  if Integer(F/d^2).is_perfect_power():
    #    good == False
  if good:
    pass
    #print b, "good"
  else:
    print b, "bad"
