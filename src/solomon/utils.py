import flint
import numpy as np

def convert_flint_to_np(flint_mat):
    cols = flint_mat.ncols()
    rows = flint_mat.nrows()
    entries = [[int(flint_mat[i,j]) for j in range(0,cols)] for i in range(0,rows)]
    return np.matrix(entries)

def convert_np_to_flint(np_mat,N):
    entries = [round(x) for x in np_mat.flatten().tolist()[0]]
    return flint.nmod_mat(np_mat.shape[0], np_mat.shape[1], entries ,N)

def primes(n):
    out = list()
    sieve = [True] * (n+1)
    for p in range(2, n+1):
        if (sieve[p]):
            out.append(p)
            for i in range(p, n+1, p):
                sieve[i] = False
    return out

def cofactors(A):
    U,sigma,Vt = np.linalg.svd(A)
    N = len(sigma)
    g = np.tile(sigma,N)
    g[::(N+1)] = 1
    G = np.diag((-1)**N * np.prod(np.reshape(g,(N,N)),1))
    return U @ G @ Vt
