
#include <algorithm>
#include <cmath>
#include <limits>

#include <boost/log/trivial.hpp>

#include <Problem.hpp>

namespace lp {

std::ostream& operator<<(std::ostream& out, const Problem& problem) {
  using namespace std;

  out << endl << "##################### Problem Start" << endl;
  out << "Value of c: " << endl << problem.c.rows() << endl;
  out << "Value of G: " << endl << "Rows: " << problem.G.rows()
      << " Cols: " << problem.G.cols() << endl;
  out << "Value of h: " << endl << problem.h.rows() << endl;
  out << "##################### Problem End" << endl;

  return out;
}

NewtonDirection operator*(const double lhs, const NewtonDirection& rhs) {
  NewtonDirection result;

  result.x = lhs * rhs.x;
  result.z = lhs * rhs.z;

  return result;
}

std::ostream& operator<<(std::ostream& out, const NewtonDirection& direction) {
  using namespace std;

  out << endl << "##################### NewtonDirection Start" << endl;
  out << "Value of x: " << endl << direction.x << endl;
  out << "Value of z: " << endl << direction.z << endl;
  out << "##################### NewtonDirection End" << endl;

  return out;
}

/**
 * Scaling structure to hold D and E
 */
class Scalings {
 public:
  Scalings(const Problem& problem) : D(problem.G.rows()), E(problem.G.cols()) {}

  Eigen::VectorXd D;
  Eigen::VectorXd E;

  double rho;
  double sigma;
};

/**
 * MF => Matrix Free Scalings taken from "ALGORITHMS FOR THE EQUILIBRATION OF
 * MATRICES AND THEIR APPLICATION TO LIMITED-MEMORY QUASI-NEWTON METHODS"
 *
 *
 */
Scalings getMFScalings(const Problem& problem, double iterationPercentage) {
  using namespace std;

  const double minScale = 1e-3;
  const double maxScale = 1e3;

  Scalings scalings(problem);

  // As it is percentage we have divide by 100
  int maxIterations = static_cast<int>((iterationPercentage / 100) *
                                       max(problem.G.rows(), problem.G.cols()));

  Eigen::ArrayXd r(problem.G.rows());
  r.setOnes();

  Eigen::ArrayXd c(problem.G.cols());
  c.setOnes();

  for (int k = 0; k < maxIterations; ++k) {
    double omega = pow(2, -max(min(floor(log2(k)) - 1.0, 4.0), 1.0));

    // Computation of D
    Eigen::ArrayXd rRand(problem.G.cols());
    rRand.setRandom();
    rRand = rRand.abs() / c.sqrt();

    Eigen::ArrayXd rRhs = (problem.G * rRand.matrix()).array();
    rRhs = rRhs.square();

    r = ((1 - omega) * r / r.sum()) + (omega * rRhs / rRhs.sum());

    // Computation of E
    Eigen::ArrayXd cRand(problem.G.rows());
    cRand.setRandom();
    cRand = cRand.abs() / r.sqrt();

    Eigen::ArrayXd cRhs = (problem.G.transpose() * cRand.matrix()).array();
    cRhs = cRhs.square();

    c = ((1 - omega) * c / c.sum()) + (omega * cRhs / cRhs.sum());
  }

  scalings.D = r.pow(-0.5).matrix();
  scalings.E = c.pow(-0.5).matrix();

  return scalings;
}

/**
 * Row mean norm
 */
double getRho(const Problem& problem) {
  double rowNorm = 0;
  for (int col = 0; col < problem.G.cols(); ++col) {
    double rowNormSquare = 0;
    for (int row = 0; row < problem.G.rows(); ++row) {
      double value = problem.G.coeff(row, col);

      rowNormSquare += value * value;
    }

    rowNorm += std::sqrt(rowNormSquare);
  }

  return rowNorm / problem.G.rows();
}

// FIXME Wrong row and col reversed
double getSigma(const Problem& problem) {
  double colNorm = 0;
  for (int row = 0; row < problem.G.rows(); ++row) {
    double colNormSquare = 0;
    for (int col = 0; col < problem.G.cols(); ++col) {
      double value = problem.G.coeff(row, col);

      colNormSquare += value * value;
    }

    colNorm += std::sqrt(colNormSquare);
  }

  return colNorm / problem.G.cols();
}

Eigen::ArrayXd getRow2Norm(const SparseMatrix& sp) {

  Eigen::ArrayXd row2Norms(sp.rows());

  for (int row = 0; row < sp.rows(); ++row) {
    double rowNormSquare = 0;
    for (int col = 0; col < sp.cols(); ++col) {
      double value = sp.coeff(row, col);

      rowNormSquare += value * value;
    }

    row2Norms(row) = std::sqrt(rowNormSquare);
  }

  return row2Norms;
}

Eigen::ArrayXd getCol2Norm(const SparseMatrix& sp) {

  Eigen::ArrayXd col2Norms(sp.cols());

  for (int col = 0; col < sp.cols(); ++col) {
    double colNormSquare = 0;
    for (int row = 0; row < sp.rows(); ++row) {
      double value = sp.coeff(row, col);

      colNormSquare += value * value;
    }

    col2Norms(col) = std::sqrt(colNormSquare);
  }

  return col2Norms;
}

Scalings getRuizScalings(const Problem& problem) {

  Eigen::ArrayXd D(problem.G.rows());
  D.setOnes();

  Eigen::ArrayXd E(problem.G.cols());
  E.setOnes();

  SparseMatrix GHat = problem.G;
  Eigen::ArrayXd row2Norm = getRow2Norm(GHat);
  Eigen::ArrayXd col2Norm = getCol2Norm(GHat);

  double r1 = row2Norm.maxCoeff() / row2Norm.minCoeff();
  double r2 = col2Norm.maxCoeff() / col2Norm.minCoeff();

  const double mByn = std::pow((GHat.rows() / GHat.cols()), (1 / 4));

  double tolerance = 1e-2;

  while (r1 > tolerance || r2 > tolerance) {

    D = D * row2Norm.pow(-0.5);

    BOOST_LOG_TRIVIAL(info) << "D " << D;

    E = E * mByn * col2Norm.pow(-0.5);

    GHat = D.matrix().asDiagonal() * problem.G * E.matrix().asDiagonal();

    row2Norm = getRow2Norm(GHat);
    col2Norm = getCol2Norm(GHat);

    r1 = row2Norm.maxCoeff() / row2Norm.minCoeff();
    r2 = col2Norm.maxCoeff() / col2Norm.minCoeff();

    BOOST_LOG_TRIVIAL(info) << "r1 " << r1;
    BOOST_LOG_TRIVIAL(info) << "r2 " << r2;
  }

  Scalings scalings(problem);

  scalings.D = D.matrix();
  scalings.E = E.matrix();

  return scalings;
}

Eigen::ArrayXd getRowInfNorm(const SparseMatrix& sp) {

  Eigen::ArrayXd rowInfNorms(sp.rows());

  for (int row = 0; row < sp.rows(); ++row) {
    double rowNorm = 0;
    for (int col = 0; col < sp.cols(); ++col) {
      double value = std::abs(sp.coeff(row, col));
      rowNorm = rowNorm < value ? value : rowNorm;
    }

    rowInfNorms(row) = rowNorm == 0 ? 1 : rowNorm;
  }

  return rowInfNorms;
}

Eigen::ArrayXd getColInfNorm(const SparseMatrix& sp) {

  Eigen::ArrayXd colInfNorms(sp.cols());

  for (int col = 0; col < sp.cols(); ++col) {
    double colNorm = 0;
    for (int row = 0; row < sp.rows(); ++row) {
      double value = std::abs(sp.coeff(row, col));
      colNorm = colNorm < value ? value : colNorm;
    }

    colInfNorms(col) = colNorm == 0 ? 1 : colNorm;
  }

  return colInfNorms;
}

Scalings getRuizInfNormScalings(const Problem& problem) {

  SparseMatrix GHat = problem.G;

  Eigen::ArrayXd D(GHat.rows());
  D.setOnes();
  Eigen::ArrayXd E(GHat.cols());
  E.setOnes();

  const double tolerance = 1e-3;

  while (true) {
    Eigen::ArrayXd rowInfNorm = getRowInfNorm(GHat);
    Eigen::ArrayXd colInfNorm = getColInfNorm(GHat);

    double r = (1 - rowInfNorm).abs().maxCoeff();
    double c = (1 - colInfNorm).abs().maxCoeff();

    if (r <= tolerance && c <= tolerance) {
      break;
    }

    rowInfNorm = rowInfNorm.pow(-0.5);
    colInfNorm = colInfNorm.pow(-0.5);

    D = D * rowInfNorm;
    E = E * colInfNorm;

    GHat = rowInfNorm.matrix().asDiagonal() * GHat *
           colInfNorm.matrix().asDiagonal();
  }

  Scalings scalings(problem);

  scalings.D = D.matrix();
  scalings.E = E.matrix();

  return scalings;
}

void Problem::normalize() {

  Scalings scalings = getRuizInfNormScalings(*this);

  this->G = scalings.D.asDiagonal() * this->G * scalings.E.asDiagonal();

  this->c = scalings.E.asDiagonal() * this->c;

  this->h = scalings.D.asDiagonal() * this->h;
}
}