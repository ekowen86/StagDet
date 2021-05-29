# stag_det_su2.py

import sys
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu


# default parameters
N = 6
m = 0.5
d = 4

# get parameters from command line if provided
if len(sys.argv) > 1:
	N = int(sys.argv[1]) # lattice size
if len(sys.argv) > 2:
	m = float(sys.argv[2]) # fermion mass
if len(sys.argv) > 3:
	d = int(sys.argv[3]) # number of dimensions

print("N: %d" % (N))
print("m: %.3f" % (m))
print("d: %d" % (d))

V = N**d  # lattice volume


# generate V random SU(2) variables
def random_links(V):
	a = 1.0 - 2.0 * np.random.rand(V)
	b = np.random.rand(V)
	c = np.random.rand(V)
	g = np.sqrt(1.0 - a * a)
	theta = np.arccos(2.0 * b - 1.0)
	phi = 2.0 * np.pi * c
	a0 = a
	a1 = g * np.sin(theta) * np.cos(phi)
	a2 = g * np.sin(theta) * np.sin(phi)
	a3 = g * np.cos(theta)
	return np.stack((a0, a1, a2, a3))


# generate a random field configuration
# indices order is mu, color, x
field = np.full((d,4,V), random_links(V))

# generate the nonzero elements of the staggered Dirac operator
row = []
col = []
data = []

def append_element(i_r, i_c, value):
	row.append(i_r)
	col.append(i_c)
	data.append(value)

for i in range(0, V):
	# get coordinates from site index
	r = i  # remainder
	x = np.empty(d, dtype=np.int32)
	for n in range(0, d):
		x[n] = r % N
		r = r // N
	
	for mu in range(0, d):
		# get forward and backward site positions
		xp1 = x.copy()
		xp1[mu] = (xp1[mu] + 1) % N
		xm1 = x.copy()
		xm1[mu] = (xm1[mu] - 1 + N) % N

		# get forward and backward site indices
		i_xp1 = 0
		i_xm1 = 0
		for n in range(0, d):
			i_xp1 += xp1[n] * pow(N, n)
			i_xm1 += xm1[n] * pow(N, n)

		# get the eta factor
		eta = 1.0
		for n in range(0, mu):
			if (x[n] & 1):
				eta *= -1.0

		# U term
		c_xp1 = 0.5 * eta * field[mu,:,i]
		append_element(i * 2 + 0, i_xp1 * 2 + 0, c_xp1[0])
		append_element(i * 2 + 0, i_xp1 * 2 + 1, c_xp1[1])
		append_element(i * 2 + 1, i_xp1 * 2 + 0, c_xp1[2])
		append_element(i * 2 + 1, i_xp1 * 2 + 1, c_xp1[3])

		# U^dagger term (conjugate, swap 1 and 2 to get transpose)
		c_xm1 = -0.5 * eta * np.conj(field[mu,:,i_xm1])
		c_xm1[1], c_xm1[2] = c_xm1[2], c_xm1[1]
		if (mu == (d - 1) and x[mu] == 0):
			# antiperiodic boundary conditions in t direction
			c_xm1 *= -1.0
		append_element(i * 2 + 0, i_xm1 * 2 + 0, c_xm1[0])
		append_element(i * 2 + 0, i_xm1 * 2 + 1, c_xm1[1])
		append_element(i * 2 + 1, i_xm1 * 2 + 0, c_xm1[2])
		append_element(i * 2 + 1, i_xm1 * 2 + 1, c_xm1[3])

	# mass term
	append_element(i * 2 + 0, i * 2 + 0, m)
	append_element(i * 2 + 1, i * 2 + 1, m)

start_time = time.time()

# create a sparse matrix and find the log determinant
D = csc_matrix((data, (row, col)), dtype=np.complex128)

# sparse version
LU = splu(D)
diagU = LU.U.diagonal().astype(np.complex128)
logdet = np.log(diagU).sum()  # det(L) = 1 by construction

# dense version (slower)
# w, v = np.linalg.eig(D.todense())
# logdet = np.sum(np.log(w))

end_time = time.time()

print("logdet: %.12f + %.12f * i * pi" % (np.real(logdet), np.imag(logdet) / np.pi))
print("time: %.12fs" % (end_time - start_time))
