#ifndef ADMM_CHOLESKY_LDLT_HPP
#define ADMM_CHOLESKY_LDLT_HPP

#include <boost/log/trivial.hpp>

#include <Eigen/PaStiXSupport>

#include <Core/Point.hpp>
#include <Problem.hpp>

namespace lp {

/**
 * When solve method is called Point is returned
 *
 */
template <typename CholeskySolverLDLT>
class ADMMCholeskyLDLT {

 public:
  // Call factor, which is done only once
  // Build lower/upper matrice (quasi definite matrice) to factor into LDLT
  ADMMCholeskyLDLT(const Problem& problem)
      : problem(problem),
        rhsF(problem.c.rows() + problem.h.rows()),
        MInverseF(rhsF.rows()) {

    linearSolver.compute(createMHat());
    if (linearSolver.info() != Eigen::Success) {
      BOOST_LOG_TRIVIAL(error)
          << "Factorization failed...Matrix is not positive definite";
      // FIXME As of now I do not know what to do, if matrix is not positive
      // definite
      throw std::runtime_error(
          "Factorization failed...Matrix is not positive definite");
    }

    rhsF << problem.c, problem.h;

    MInverseF = linearSolver.solve(rhsF);
    // Negate y part as we have normalized M to MHat to confirm with spd
    MInverseF.tail(problem.h.rows()) = -MInverseF.tail(problem.h.rows());

    // Denominator of Matrix inverse lemma
    denominator = 1 + rhsF.transpose() * MInverseF;
  }

  /**
   *
   */
  Point doProjection(const Point& point) const {

    Eigen::VectorXd deltaXTauC =
        point.x * xScalingParameter - point.tau * problem.c;
    Eigen::VectorXd deltaYTauB = point.z - point.tau * problem.h;
    // RHS of Matrix Inverse Lemma (mil)
    Eigen::VectorXd milRhs(problem.c.rows() + problem.h.rows());
    milRhs << deltaXTauC, deltaYTauB;

    Eigen::VectorXd MInverseMilRhs = linearSolver.solve(milRhs);

    // Negate y part as we have normalized M to MHat to confirm with spd
    MInverseMilRhs.tail(problem.h.rows()) =
        -MInverseMilRhs.tail(problem.h.rows());

    // (mil) Matrix Inverse Lemma solution
    // FIXME MInverseF * rhsF.transpose() can be cached, but I am getting
    // operator= no found for the arguments when I do in constructor
    Eigen::VectorXd milSolution =
        MInverseMilRhs -
        (MInverseF * rhsF.transpose() * MInverseMilRhs) / denominator;

    Point projectedPoint;

    projectedPoint.x = milSolution.head(problem.G.cols());
    projectedPoint.z = milSolution.tail(problem.G.rows());
    // tau + c' * x + h' * z
    projectedPoint.tau = point.tau + rhsF.transpose() * milSolution;

    return projectedPoint;
  }

 private:
  /**
   * MHat = [ I  -A' ] [x] = [rhsX]
   *        [-A  -I  ] [y] = [rhsY]
   * We are going to build only lower matrix
   */
  SparseMatrix createMHat() const {

    int GRows = problem.G.rows();
    int Gcols = problem.G.cols();
    int matrixSize = GRows + Gcols;

    SparseMatrix lowerMHat(matrixSize, matrixSize);
    // FIXME Check the performace of this loops
    // FIXME Fix the loops for better code readability
    // As storage schem is column compressed, we fill data column wise to
    // achieve O(1) complexity. And we know it is symmetric so we fill only
    // lower traingular matrix
    for (int colIndex = 0; colIndex < matrixSize; ++colIndex) {
      // Change the rowIndex as we need to traverse only lower triangle, thats
      // the reason for rowIndex = colIndex
      for (int rowIndex = colIndex; rowIndex < matrixSize; ++rowIndex) {
        // fill top left blocks diagonal
        if (colIndex == rowIndex && colIndex < Gcols && rowIndex < Gcols) {
          lowerMHat.insert(rowIndex, colIndex) = xScalingParameter;
        }
        // fill bottom right block diagonal
        if (colIndex == rowIndex && colIndex >= Gcols && rowIndex >= Gcols) {
          lowerMHat.insert(rowIndex, colIndex) = -1;
        }
        // fill bottom left block using G
        if (rowIndex >= Gcols && colIndex < Gcols) {
          lowerMHat.insert(rowIndex, colIndex) =
              -problem.G.coeff(rowIndex - Gcols, colIndex);
        }
      }
    }

    return lowerMHat;
  }

  const double xScalingParameter = 1e-3;
  const Problem& problem;
  Eigen::VectorXd rhsF;
  // Denominator of Matrix Inverse Lemma
  // 1 + f' * M{-1} * f
  // where f = [c, h]'
  double denominator;
  // M{-1}*f
  Eigen::VectorXd MInverseF;
  // Caching another constant which is reused during solve phase
  // M{-1} * h * h'
  //   SparseMatrix MInverseFFTranspose;

  CholeskySolverLDLT linearSolver;
};
}
#endif  // ADMM_CHOLESKY_LDLT_HPP