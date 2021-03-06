# stag_det.py

import sys
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

# default parameters
N = 64
m = 0.5

# get parameters from command line if provided
if len(sys.argv) == 3:
	N = int(sys.argv[1]) # lattice size
	m = float(sys.argv[2]) # fermion mass

print("N: %d" % (N))
print("m: %.3f" % (m))

# generate a random configuration
field1 = np.exp(2j * np.pi * np.random.rand(N * N))
field2 = np.exp(2j * np.pi * np.random.rand(N * N))

# generate the nonzero elements of the staggered Dirac operator
row = []
col = []
data = []
for i in range(0, len(field1)):
	x = i % N
	y = i // N

	xp1 = (x + 1) % N
	xm1 = (x - 1 + N) % N
	yp1 = (y + 1) % N
	ym1 = (y - 1 + N) % N

	i_xp1 = xp1 + y * N
	i_xm1 = xm1 + y * N
	i_yp1 = x + yp1 * N
	i_ym1 = x + ym1 * N

	c_xp1 = 0.5 * field1[i]
	row.append(i)
	col.append(i_xp1)
	data.append(c_xp1)

	c_xm1 = -0.5 * np.conj(field1[i_xm1])
	row.append(i)
	col.append(i_xm1)
	data.append(c_xm1)

	c_yp1 = 0.5 * field2[i]
	if (x & 1 != 0):
		c_yp1 *= -1.0  # eta factor
	if (yp1 == 0):
		c_yp1 *= -1.0  # antiperiodic boundary conditions in y direction
	row.append(i)
	col.append(i_yp1)
	data.append(c_yp1)

	c_ym1 = -0.5 * np.conj(field2[i_ym1])
	if (x & 1 != 0):
		c_ym1 *= -1.0  # eta factor
	if (y == 0):
		c_ym1 *= -1.0  # antiperiodic boundary conditions in y direction
	row.append(i)
	col.append(i_ym1)
	data.append(c_ym1)

	# mass term
	row.append(i)
	col.append(i)
	data.append(m)

start_time = time.time()

# create a sparse matrix and find the log determinant
D = csc_matrix((data, (row, col)), dtype=np.complex128)
LU = splu(D)
diagU = LU.U.diagonal().astype(np.complex128)
logdet = np.log(diagU).sum()  # det(L) = 1 by construction

# measure chiral condensate (single point source at origin)
src = np.zeros(N * N)
src[0] = 1.0
cc = np.real(LU.solve(src)[0])
print("chiral cond.: %.12f (single)" % cc)

# measure chiral condensate (entire volume, slow)
# cc_sum = 0.0
# for i in range(0, N * N):
# 	src = np.zeros(N * N)
# 	src[i] = 1.0
# 	cc_sum += np.real(LU.solve(src)[i])
# print("chiral cond.: %.12f (all)" % (cc_sum / (N * N)))

end_time = time.time()

print("logdet: %.12f + %.12f * i * pi" % (np.real(logdet), np.imag(logdet) / np.pi))
print("time: %.12fs" % (end_time - start_time))
