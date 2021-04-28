# StagDet

This program checks that the determinant of the staggered Dirac operator is
invariant under translations and rotations on a 2D square lattice with a U(1)
gauge field and periodic boundary conditions.

First it generates a random gauge field configuration (base case), then it
generates several additional configurations by transforming the base case in a
few different ways:

- x+1: translate by 1 in the positive x-direction
- x+2: translate by 2 in the positive x-direction
- y+1: translate by 1 in the positive y-direction
- y+2: translate by 2 in the positive y-direction
- xflip: reflect in the x-direction (about y=N/2)
- yflip: reflect in the y-direction (about x=N/2)
- xyswap: swaps the x and y indices
- r90: rotate CCW by 90 degrees
- r180: rotate CCW by 180 degrees
- r270: rotate CCW by 270 degrees

Next, it computes the average plaquette for each configuration. This is done to
check that the transformation was applied correctly. Finally, it uses the
Eigen/Sparse library to compute the log determinant of the staggered Dirac
operator for each configuration.

## Usage

    stag_det N m
    N: lattice size
    m: fermion mass

Note that with staggered fermions the determinant is not invariant if N is odd.

## Typical Output

    ./stag_det 60 0.4
    N: 60
    m: 0.400
    generating random U(1) gauge field configuration...

    measuring average plaquette...
    0.012204462119 (base case)
    0.012204462119 (x+1)
    0.012204462119 (x+2)
    0.012204462119 (y+1)
    0.012204462119 (y+2)
    0.012204462119 (xflip)
    0.012204462119 (yflip)
    0.012204462119 (xyswap)
    0.012204462119 (r90)
    0.012204462119 (r180)
    0.012204462119 (r270)

    calculating log(det D)...
    -377.426047102154484 (base case)
    -377.426047102155337 (x+1)
    -377.426047102158577 (x+2)
    -377.426047102155678 (y+1)
    -377.426047102157156 (y+2)
    -377.426047102155508 (xflip)
    -377.426047102156815 (yflip)
    -377.426047102156019 (xyswap)
    -377.426047102156076 (r90)
    -377.426047102157725 (r180)
    -377.426047102155678 (r270)
