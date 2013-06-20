consts = {10: 32.3, 12: 29.9, 14: 28.2, 16: 26.9, 18: 26.0, 20: 25.2, 22: 24.5, 24: 24.0, 26: 23.5, 28: 23.1, 30: 22.8}
D = 5

for (m, C1) in consts.iteritems():
  min_p = ceil(exp(m/2 - 0.21 + 2 * log(10)))
  p = min_p
  last_diff = 10**10
  while True:
    val = C1 * 4 * log(10) * (log(p) + 0.21 - log(10))**2
    print p, val.n()
    if val.n() > p:
      print "For m =", m, "p <", p
      break


    diff = p - val
    if diff > last_diff:
      print "No success for m =", m
      break
    last_diff = diff

    p += 1
