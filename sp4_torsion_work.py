from sage.matrix.special import block_matrix
from sage.matrix.constructor import matrix
from sage.matrix.constructor import identity_matrix
from operator import mul
from functools import reduce
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.arith.misc import euler_phi
from sage.combinat.partition import Partitions
import collections
import itertools

# Much of this is based on Qingjie Yang's thesis, hereby abbreviated as QY.
# Here is a link to this thesis.
# https://open.library.ubc.ca/soa/cIRcle/collections/ubctheses/831/items/1.0079840

U = matrix([[1, 0], [1, -1]])

Y = matrix([[0, 1, 0, 0], [0, 0, -1, 0], [0, 0, -1, 1], [1, 1, -1, 0]])

Z = matrix([[0, -1, 1, 0], [-1, 0, 1, 1], [-1, 1, 0, 0], [0, -1, 0, 0]])


def V(x, y, z):
    return matrix(
        [[2 * x, 0, -y, x], [0, -2 * x, -x, -z], [z, x, -x, z], [-x, y, y, x]]
    )


def ident(n):
    return identity_matrix(n)


def block_it(X):
    n = int(X.nrows() / 2)
    A = X.matrix_from_rows_and_columns(range(0, n), range(0, n))
    B = X.matrix_from_rows_and_columns(range(n, 2 * n), range(0, n))
    C = X.matrix_from_rows_and_columns(range(0, n), range(n, 2 * n))
    D = X.matrix_from_rows_and_columns(range(n, 2 * n), range(n, 2 * n))
    return (A, B, C, D)


# See (1.3) of QY, also standard
def direct_sum(X, Y):
    return block_matrix([[X, 0], [0, Y]], subdivide=False)


# See (1.7) of QY
def symp_direct_sum(X, Y):
    (A1, B1, C1, D1) = block_it(X)
    (A2, B2, C2, D2) = block_it(Y)
    return block_matrix(
        [[A1, 0, B1, 0], [0, A2, 0, B2], [C1, 0, D1, 0], [0, C2, 0, D2]],
        subdivide=False,
    )


def J(n):
    k = int(n / 2)
    return block_matrix([[0, ident(k)], [-ident(k), 0]], subdivide=False)


def W(lam):
    return matrix([[0, -1], [1, -lam]])

sl2_torsion = {
    2: [-ident(2)],
    3: [matrix([[0,-1],[1,-1]])],
    4: [matrix([[0, -1], [1,0]])],
    6: [matrix([[1,-1],[1,0]])]
}

# Torsion in Sp(2g; Z) up to conjugacy which is realized by
# a cyclic action on a Riemann surface.
# (See proposition 6.7 of QY)
sp4_torsion_2 = [-ident(4), direct_sum(U, U.transpose())]
sp4_torsion_3 = [symp_direct_sum(W(1), W(1).transpose())]
sp4_torsion_4 = [symp_direct_sum(-J(2), J(2))]
sp4_torsion_5 = [Y, Y**2, Y**3, Y**4]
sp4_torsion_6 = [-sp4_torsion_3[0], V(0, 1, 1), V(0, -1, -1)]
sp4_torsion_8 = [Z, -Z]
sp4_torsion_10 = [-Y, -(Y**2), -(Y**3), -(Y**4)]
sp4_torsion = {
    2: sp4_torsion_2,
    3: sp4_torsion_3,
    4: sp4_torsion_4,
    5: sp4_torsion_5,
    6: sp4_torsion_6,
    8: sp4_torsion_8,
    10: sp4_torsion_10,
}

def basis_vector(i,dim):
    return matrix(ident(dim).column(i)).transpose()

def symp(x, y, JJ=None):
    if JJ is None:
        JJ = J(x.nrows())
    return (x.transpose() * JJ * y)[0, 0]


def dehn_twist(v, x, JJ=None):
    return x + symp(x, v, JJ=JJ) * v


def dehn_twist_mat(v, JJ=None):
    if JJ is None:
        JJ = J(v.nrows())
    II = ident(v.nrows())
    return II + v * (v.transpose() * JJ.transpose())


def primitive_factor_for_conj(conj, period, R = None, dim=None):
    if dim is None:
        dim = conj.nrows()
    if R is None:
        R = PolynomialRing(QQ, "x", dim)
    x = R.gens()
    v = matrix([[x[i]] for i in range(0, dim)])
    ts = [dehn_twist_mat(conj**k * v) for k in range(0, period)]
    bigT = reduce(mul, ts, ident(v.nrows()))
    return bigT

def calc_char_poly_conj(conj,period,dim=None):
    return primitive_factor_for_conj(conj,period,dim).charpoly('t')


def possible_cycs(dim):
    R = ZZ["t"]
    factors_idx = [
        (i, euler_phi(i)) for i in range(1, 2 * (dim**2) + 1) if euler_phi(i) <= dim
    ]
    factors = [(j, R.cyclotomic_polynomial(i)) for (i, j) in factors_idx]
    factor_dict = collections.defaultdict(list)
    for k, v in factors:
        factor_dict[k].append(v)
    return factor_dict


def possible_torsion_char(dim):
    factor_dict = possible_cycs(dim)
    possible = []
    for part_lst in Partitions(dim).list():
        factor_iter = itertools.product(*[factor_dict[k] for k in part_lst])
        new_possible = [reduce(mul, factor, 1) for factor in factor_iter]
        possible += new_possible
    return set(possible)


# For the analogue in QY, see Proposition 6.2 and equations (6.2)-(6.17)
def possible_symp_torsion_char(dim):
    assert dim % 2 == 0, "Symplectic Space must have even dimension"
    return [
        p
        for p in possible_torsion_char(dim)
        if p.coefficients() == p.coefficients()[::-1]
    ]

-x1^6 + 3*x1^4*y1^2 - 3*x1^2*y1^4 + y1^6 - 3*x1^5*x2 + 6*x1^3*y1^2*x2 - 3*x1*y1^4*x2 - 6*x1^4*x2^2 + 9*x1^2*y1^2*x2^2 - 3*y1^4*x2^2 - 7*x1^3*x2^3 + 6*x1*y1^2*x2^3 - 6*x1^2*x2^4 + 3*y1^2*x2^4 - 3*x1*x2^5 - x2^6 - 3*x1^4*y1*y2 + 6*x1^2*y1^3*y2 - 3*y1^5*y2 - 6*x1^3*y1*x2*y2 + 6*x1*y1^3*x2*y2 - 9*x1^2*y1*x2^2*y2 + 6*y1^3*x2^2*y2 - 6*x1*y1*x2^3*y2 - 3*y1*x2^4*y2 + 3*x1^4*y2^2 - 9*x1^2*y1^2*y2^2 + 6*y1^4*y2^2 + 6*x1^3*x2*y2^2 - 9*x1*y1^2*x2*y2^2 + 9*x1^2*x2^2*y2^2 - 9*y1^2*x2^2*y2^2 + 6*x1*x2^3*y2^2 + 3*x2^4*y2^2 + 6*x1^2*y1*y2^3 - 7*y1^3*y2^3 + 6*x1*y1*x2*y2^3 + 6*y1*x2^2*y2^3 -3*x1^2*y2^4 + 6*y1^2*y2^4 - 3*x1*x2*y2^4 - 3*x2^2*y2^4 - 3*y1*y2^5 + y2^6 + 3*x1^4- 6*x1^2*y1^2 + 3*y1^4 + 6*x1^3*x2 - 6*x1*y1^2*x2 + 9*x1^2*x2^2 - 6*y1^2*x2^2 + 6*x1*x2^3 + 3*x2^4 + 6*x1^2*y1*y2 - 6*y1^3*y2 + 6*x1*y1*x2*y2 + 6*y1*x2^2*y2 - 6*x1^2*y2^2 + 9*y1^2*y2^2 - 6*x1*x2*y2^2 - 6*x2^2*y2^2 - 6*y1*y2^3 + 3*y2^4 - 4
