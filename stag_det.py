# stag_det.py

import sys
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

N = int(sys.argv[1]) # lattice size
print("N: %d" % (N))

m = float(sys.argv[2]) # fermion mass
print("m: %.3f" % (m))

# generate a random configuration
field1 = np.exp(2j * np.pi * np.random.rand(N * N))
field2 = np.exp(2j * np.pi * np.random.rand(N * N))

start_time = time.time()

# generate the nonzero elements of the staggered Dirac operator
row = []
col = []
data = []
for i in range(0, len(field1)):
	x = i % N
	y = i / N

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
		c_yp1 *= -1.0;  # eta factor
	row.append(i)
	col.append(i_yp1)
	data.append(c_yp1)

	c_ym1 = -0.5 * np.conj(field2[i_ym1])
	if (x & 1 != 0):
		c_ym1 *= -1.0;  # eta factor
	row.append(i)
	col.append(i_ym1)
	data.append(c_ym1)

	# mass term
	row.append(i)
	col.append(i)
	data.append(m)

# create a sparse matrix and find the log determinant
D = csc_matrix((data, (row, col)), dtype=np.complex128)
LU = splu(D)
diagL = LU.L.diagonal().astype(np.complex128)
diagU = LU.U.diagonal().astype(np.complex128)
logdet = np.log(diagL).sum() + np.log(diagU).sum()

end_time = time.time()

print("logdet: %.15f" % np.real(logdet))
print("time: %.12fs" % (end_time - start_time))
