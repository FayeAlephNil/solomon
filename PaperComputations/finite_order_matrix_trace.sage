load("lib.sage")

C4 = matrix([[0, -1], [1,0]])
C3 = matrix([[0,-1],[1,-1]])
C6 = matrix([[1,-1],[1,0]])
var('p q')
Tv = matrix_from_scc(matrix([[p],[q]]))
f4 = (Tv*C4*Tv*C4**(-1)).trace()
f3 = (Tv * C3 * Tv * C3**(-1) * C3**2 * Tv * C3**(-2)).trace()
f6 = (Tv * C6 * Tv * C6**(-1) * C6**2 * Tv * C6**(-2)).trace()

assert (f3 - f3_cand).full_simplify() == 0
assert (f6 - f6_cand).full_simplify() == 0
assert (f4 - f4_cand).full_simplify() == 0
