for b in xrange(1,72,2):
  good = True
  n = 2
  u = 1
  uprev = 0
  while not (u == 0 and uprev == 1):
    temp = uprev
    uprev = u
    u = (temp + b*u)%9
    #print u
    if u==0 and n%3==0 and n%2==1:
      good = False
    n += 1
  if good:
    print b, "good", n-1
  else:
    print b, "bad", n-1
