#ifndef ADMM_CHOLESKY_LDLT_HPP
#define ADMM_CHOLESKY_LDLT_HPP

#include <vector>

#include <boost/log/trivial.hpp>

#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include <Core/Timer.hpp>
#include <Core/Point.hpp>

#include <Problem.hpp>

namespace lp {

class ADMMCholeskyLDLT {
 public:
  ADMMCholeskyLDLT(const Problem& problem)
      : problem(problem), rhsEta(problem.c.rows() + problem.h.rows()) {

    // AMD and Metis // FIXME Check this setting!!!
    //     choleskySolver.cholmod().nmethods = 2;

    timer.start("Create MHat");
    SparseMatrix mHat = createMHat();
    timer.end("Create MHat");

    //     BOOST_LOG_TRIVIAL(info) << std::endl <<mHat;

    // TODO May be use std::move for mHat as we do not use it anymore
    timer.start("Cholesky Factorize");
    choleskySolver.compute(mHat);
    timer.end("Cholesky Factorize");

    if (choleskySolver.info() != Eigen::Success) {
      BOOST_LOG_TRIVIAL(error)
          << "Factorization failed...Matrix is not positive definite";
      // FIXME As of now I do not know what to do, if matrix is not positive
      // definite
      throw std::runtime_error(
          "Factorization failed...Matrix is not positive definite");
    }

    // rhsEta =  [c; h]
    rhsEta << problem.c, problem.h;

    mInverseEta = solve(rhsEta);

    denominator = 1 + rhsEta.transpose() * mInverseEta;
  }

  void doSubspaceProjection(Point& point) {
    // TODO These two can be parallelized, but is it worth
    Eigen::VectorXd xTauC = point.x * problem.admmRho - point.tau * problem.c;
    Eigen::VectorXd zTauH = point.z - point.tau * problem.h;

    Eigen::VectorXd rhs(problem.c.rows() + problem.h.rows());
    rhs << xTauC, zTauH;

    BOOST_LOG_TRIVIAL(info) << "rhs" << rhs;

    Eigen::VectorXd mInverseRhs = solve(rhs);

    // (mil) Matrix Inverse Lemma solution
    Eigen::VectorXd milSolution =
        mInverseRhs -
        (mInverseEta * (rhsEta.transpose() * mInverseRhs)) / denominator;

    point.x = milSolution.head(problem.c.rows());
    point.z = milSolution.tail(problem.h.rows());
    point.tau = point.tau + problem.c.transpose() * point.x +
                problem.h.transpose() * point.z;
  }

 private:
  Eigen::VectorXd solve(const Eigen::VectorXd& rhs) {
    Eigen::VectorXd lhs = choleskySolver.solve(rhs);
    lhs.tail(problem.h.rows()) = -lhs.tail(problem.h.rows());
    return lhs;
  }
  /**
   * Find a better way to insert elements into Eigen SparseMatrix instead of
   * using Triplets, direct insert for some reason is slower than
   * setFromTriplets approach
   */
  SparseMatrix createMHat() {
    std::vector<Eigen::Triplet<double>> mHatTriplets;
    // Number of non zeros in final resultant SparseMatrix
    // [ I*rho  -G']
    // [ -G      I ]
    // is nz in I*rho + nz in I + nz in G
    mHatTriplets.reserve(problem.G.cols() + problem.G.rows() +
                         problem.G.nonZeros());

    // First fill diagonal entries i.e. I and I*rho
    // Final Matrix size is (m+n) X (m+n)
    // So there are m+n diagonal entries, I*rho dimension is n X n and I
    // dimension is m X m
    for (int i = 0; i < problem.G.cols(); ++i) {
      mHatTriplets.emplace_back(i, i, problem.admmRho);
    }

    for (int i = problem.G.cols(); i < problem.G.cols() + problem.G.rows();
         ++i) {
      mHatTriplets.emplace_back(i, i, -1);
    }

    // Lower bottom corner which is -G
    for (int col = 0; col < problem.G.outerSize(); ++col) {
      for (SparseMatrix::InnerIterator it(problem.G, col); it; ++it) {
        mHatTriplets.emplace_back(it.row() + problem.G.cols(), it.col(),
                                  -it.value());
      }
    }

    SparseMatrix mHat(problem.G.cols() + problem.G.rows(),
                      problem.G.cols() + problem.G.rows());

    mHat.setFromTriplets(mHatTriplets.begin(), mHatTriplets.end());

    return mHat;
  }

  const Problem& problem;
  Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double, Eigen::Lower>>
      choleskySolver;

  Eigen::VectorXd rhsEta;
  // M{-1} * eta
  Eigen::VectorXd mInverseEta;
  // 1 + h' * M{-1} * h
  double denominator;

  Timer& timer = Timer::getADMMInstance();
};
}

#endif  // ADMM_CHOLESKY_LDLT_HPP