//
// stag_det.cc
//
// checks that the determinant of the staggered Dirac operator is invariant
// under translations and rotations on a 2D square lattice with U(1) a gauge
// field and periodic boundary conditions.

#include <cstdio>
#include <complex>
#include <vector>
#include <random>
#include <Eigen/Sparse>

typedef std::complex<double> Complex;
typedef std::vector<Complex> Field;
typedef Eigen::SparseMatrix<Complex> Matrix;
typedef Eigen::Triplet<Complex> Triplet;

int N;  // lattice size
double m;  // fermion mass

double CalcLogDet(Field field1, Field field2);
double CalcPlaq(Field field1, Field field2);

int main(int argc, char* argv[]) {

  N = std::atoi(argv[1]); printf("N: %d\n", N);
  m = std::stod(argv[2]); printf("m: %.3f\n", m);

  // generate a random configuration
  printf("generating random U(1) gauge field configuration...\n");
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(-M_PI, M_PI);
  Field field1(N * N);
  Field field2(N * N);
  for (int i = 0; i < field1.size(); i++) {
    field1[i] = std::polar(1.0, dist(gen));
    field2[i] = std::polar(1.0, dist(gen));
  }

  // generate shifted and flipped fields
  Field field1_x1(field1.size());
  Field field2_x1(field2.size());
  Field field1_x2(field1.size());
  Field field2_x2(field2.size());
  Field field1_y1(field1.size());
  Field field2_y1(field2.size());
  Field field1_y2(field1.size());
  Field field2_y2(field2.size());
  Field field1_xflip(field1.size());
  Field field2_xflip(field2.size());
  Field field1_yflip(field1.size());
  Field field2_yflip(field2.size());
  Field field1_xyswap(field1.size());
  Field field2_xyswap(field2.size());
  Field field1_r90(field1.size());
  Field field2_r90(field2.size());
  Field field1_r180(field1.size());
  Field field2_r180(field2.size());
  Field field1_r270(field1.size());
  Field field2_r270(field2.size());

  for (int i = 0; i < field1.size(); i++) {
    int x = i % N;
    int y = i / N;

    int xp1 = (x + 1) % N;
    int yp1 = (y + 1) % N;
    int xp2 = (x + 2) % N;
    int yp2 = (y + 2) % N;
    int xflip = N - x - 1;
    int yflip = N - y - 1;
    int xflipm1 = (xflip - 1 + N) % N;
    int yflipm1 = (yflip - 1 + N) % N;

    field1_x1[i] = field1[xp1 + y * N];
    field2_x1[i] = field2[xp1 + y * N];
    field1_x2[i] = field1[xp2 + y * N];
    field2_x2[i] = field2[xp2 + y * N];
    field1_y1[i] = field1[x + yp1 * N];
    field2_y1[i] = field2[x + yp1 * N];
    field1_y2[i] = field1[x + yp2 * N];
    field2_y2[i] = field2[x + yp2 * N];
    field1_xflip[i] = conj(field1[xflipm1 + y * N]);
    field2_xflip[i] = field2[xflip + y * N];
    field1_yflip[i] = field1[x + yflip * N];
    field2_yflip[i] = conj(field2[x + yflipm1 * N]);
    field1_xyswap[i] = field2[y + x * N];
    field2_xyswap[i] = field1[y + x * N];
    field1_r90[i] = conj(field2[y + xflipm1 * N]);
    field2_r90[i] = field1[y + xflip * N];
    field1_r180[i] = conj(field1[xflipm1 + yflip * N]);
    field2_r180[i] = conj(field2[xflip + yflipm1 * N]);
    field1_r270[i] = field2[yflip + x * N];
    field2_r270[i] = conj(field1[yflipm1 + x * N]);
  }

  printf("\nmeasuring average plaquette...\n");
  printf("%.12f (base case)\n", CalcPlaq(field1, field2));
  printf("%.12f (x+1)\n", CalcPlaq(field1_x1, field2_x1));
  printf("%.12f (x+2)\n", CalcPlaq(field1_x2, field2_x2));
  printf("%.12f (y+1)\n", CalcPlaq(field1_y1, field2_y1));
  printf("%.12f (y+2)\n", CalcPlaq(field1_y2, field2_y2));
  printf("%.12f (xflip)\n", CalcPlaq(field1_xflip, field2_xflip));
  printf("%.12f (yflip)\n", CalcPlaq(field1_yflip, field2_yflip));
  printf("%.12f (xyswap)\n", CalcPlaq(field1_xyswap, field2_xyswap));
  printf("%.12f (r90)\n", CalcPlaq(field1_r90, field2_r90));
  printf("%.12f (r180)\n", CalcPlaq(field1_r180, field2_r180));
  printf("%.12f (r270)\n", CalcPlaq(field1_r270, field2_r270));

  // calculate abs log det of dirac operator
  printf("\ncalculating log(det D)...\n");
  printf("%.15f (base case)\n", CalcLogDet(field1, field2));
  printf("%.15f (x+1)\n", CalcLogDet(field1_x1, field2_x1));
  printf("%.15f (x+2)\n", CalcLogDet(field1_x2, field2_x2));
  printf("%.15f (y+1)\n", CalcLogDet(field1_y1, field2_y1));
  printf("%.15f (y+2)\n", CalcLogDet(field1_y2, field2_y2));
  printf("%.15f (xflip)\n", CalcLogDet(field1_xflip, field2_xflip));
  printf("%.15f (yflip)\n", CalcLogDet(field1_yflip, field2_yflip));
  printf("%.15f (xyswap)\n", CalcLogDet(field1_xyswap, field2_xyswap));
  printf("%.15f (r90)\n", CalcLogDet(field1_r90, field2_r90));
  printf("%.15f (r180)\n", CalcLogDet(field1_r180, field2_r180));
  printf("%.15f (r270)\n", CalcLogDet(field1_r270, field2_r270));

  return 0;
}

double CalcLogDet(Field field1, Field field2) {

  // generate the nonzero elements of the staggered Dirac operator
  std::vector<Triplet> dirac_elements;
  for (int i = 0; i < field1.size(); i++) {

    int x = i % N;
    int y = i / N;

    int xp1 = (x + 1) % N;
    int xm1 = (x - 1 + N) % N;
    int yp1 = (y + 1) % N;
    int ym1 = (y - 1 + N) % N;

    int i_xp1 = xp1 + y * N;
    int i_xm1 = xm1 + y * N;
    int i_yp1 = x + yp1 * N;
    int i_ym1 = x + ym1 * N;

    Complex c_xp1 = 0.5 * field1[i];
    dirac_elements.push_back(Triplet(i, i_xp1, c_xp1));

    Complex c_xm1 = -0.5 * conj(field1[i_xm1]);
    dirac_elements.push_back(Triplet(i, i_xm1, c_xm1));

    Complex c_yp1 = 0.5 * field2[i];
    if (x & 1) c_yp1 *= -1.0;  // eta factor
    dirac_elements.push_back(Triplet(i, i_yp1, c_yp1));

    Complex c_ym1 = -0.5 * conj(field2[i_ym1]);
    if (x & 1) c_ym1 *= -1.0;  // eta factor
    dirac_elements.push_back(Triplet(i, i_ym1, c_ym1));

    // mass term
    dirac_elements.push_back(Triplet(i, i, m));
  }

  Matrix D(N*N,N*N);
  D.setFromTriplets(dirac_elements.begin(), dirac_elements.end());

  Eigen::SparseLU<Matrix> solver;
  solver.compute(D);
  if(solver.info() != Eigen::Success) {
    printf("compute failed!\n");
    return 0.0;
  }
  solver.factorize(D);
  if(solver.info() != Eigen::Success) {
    printf("factorize failed!\n");
    return 0.0;
  }
  Complex log_det = solver.logAbsDeterminant();
  if(solver.info() != Eigen::Success) {
    printf("logAbsDeterminant failed!\n");
    return 0.0;
  }
  return real(log_det);
}

double CalcPlaq(Field field1, Field field2) {
  Complex plaq = 0.0;
  for (int i = 0; i < field1.size(); i++) {
    int x = i % N;
    int y = i / N;

    int xp1 = (x + 1) % N;
    int yp1 = (y + 1) % N;

    Complex u1 = field1[i];
    Complex u2 = field2[xp1 + y * N];
    Complex u3 = conj(field1[x + yp1 * N]);
    Complex u4 = conj(field2[i]);

    plaq += u1 * u2 * u3 * u4;
  }
  return real(plaq) / double(field1.size());
}
