#ifndef SUITESPARSE_CHOLESKYLLT_HPP
#define SUITESPARSE_CHOLESKYLLT_HPP

#include <algorithm>
#include <cstdlib>
#include <cmath>

#include <boost/log/trivial.hpp>
#include <boost/mpl/int.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include <Core/Timer.hpp>
#include <Problem.hpp>

namespace lp {

// TODO Define new namespace
// FIXME As of now A (inequality constraints) is not considered
// rely on presolve to normalize equality to inequality constraints

/**
 *SuiteSparse cholmod
 *TODO Move below comment to IPM algorithm class
 *Q. Why call factor, why not do it in constructor of Linear solver provided
 *A. Hmmm... Both approaches have their own advantages and disadvantages, 1.
 *factor in constructor is much clean and straight but user (creator of the
 *linear solver) will not have
 *any context of previous iteration as object is destroyed at the end of each
 *iteration in IPM. But advantage is user doesn't have to delete any objects at
 *each iteration that he
 *has created to find solution, he can rely on destructor to do proper cleaning
 *in face of any disater 2. Seperate factor method gives user flexibiltiy to
 *store context information and use it later to improve performace like reuse
 *analyze, update/downdate. But user have to maintain objects created carefully.
 *As performace is utmost importance second option is chosen
 *
 */
template <typename Scalings, typename CholeskySolver>
class IPMCholeskyLLT {

 public:
  IPMCholeskyLLT(const lp::Problem& problem) : problem(problem) {
    // TODO In Linear programming non zero pattern of SPD does not change, so it
    // done only once, but what about other cones?
    Timer& timer = Timer::getIPMInstance();
    timer.start("Cholesky Analyze");
    solver.analyzePattern(problem.G.transpose() * problem.G);
    timer.end("Cholesky Analyze");
  }

  // Factor matrix omega = G' * W{-1} * W{-T} * G
  // If scalingMatrix W is Identity factor G' * G
  // If A (equality constaints) is not empty then factor A * omega{-1} * A' to
  // find y
  // First find y, then x and finally z
  // Diagonal Scaling matrix is only good for Linear programming
  // TODO Input is not anymore DiagonalMatrix but Scalings object adapt!!
  template <lp::SolveFor solveFor>
  void factor(const Scalings& scalings) {

    Timer& timer = Timer::getIPMInstance();
    // FIXME Check the possibility of avoiding muliplication, check the cholmod
    // doc
    timer.start("PSD Matrix");
    // First compute G * W{-T} and store it
    // Used also to compute omegaTilde
    omega = getOmega(scalings, boost::mpl::int_<static_cast<int>(solveFor)>());

    // Used in solve method, cached as it is used multiple times
    omegaTilde =
        getOmegaTilde(scalings, boost::mpl::int_<static_cast<int>(solveFor)>());

    // G' * W{-1} * W{-T} * G
    SparseMatrix M = omega.transpose() * omega;

    timer.end("PSD Matrix");

    // SparseMatrix MHat = scaleSymmetricMatrix(M);

    timer.start("Cholesky Factorize");
    //     solver.compute(MHat);
    solver.compute(M);
    timer.end("Cholesky Factorize");

    if (solver.info() != Eigen::Success) {
      BOOST_LOG_TRIVIAL(error)
          << "Factorization failed...Matrix is not positive definite";
      // FIXME As of now I do not know what to do, if matrix is not positive
      // definite
      throw std::runtime_error(
          "Factorization failed...Matrix is not positive definite");
    }
  }

  // rhsX, rhsY, rhsZ are not primal/dual variables but symbolic names of right
  // hand side of three equations
  // When A(equality constaints) is empty rhsY is not required/considered
  // When problem is in standard for without inequality constraints rhs has only
  // two vectors
  // First calculate y from rhs
  // (A * L{-T} * L{-1} * A') * y = A * L{-T} * L{-1} * (rhsX + G' * W{-1} *
  // W{-T} * rhsZ + A' * rhsY) - rhsY
  // then calculate x from rhs
  // (L * L') * x = rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
  // finally calculate z
  // z = W{-1} * W{-T} * (G*x - rhsZ)
  // TODO Error message to check if compute/factor is called before solve
  template <lp::SolveFor solveFor>
  lp::NewtonDirection solve(const Eigen::VectorXd& rhsX,
                            const Eigen::VectorXd& rhsZ,
                            const Scalings& scalings) const {

    lp::NewtonDirection direction;

    {
      Timer& timer = Timer::getIPMInstance();
      timer.start("Cholesky Solve");
      // rhsX + G' * W{-1} * W{-T} * rhsZ
      Eigen::VectorXd rhs = /*diagonalMatrix **/(rhsX + omegaTilde * rhsZ);

      direction.x = /*diagonalMatrix **/ solver.solve(rhs);
      timer.end("Cholesky Solve");
    }
    {
      direction.z = computeWZ(scalings, direction, rhsZ,
                              boost::mpl::int_<static_cast<int>(solveFor)>());
    }

    return direction;
  }

 private:
  // G
  SparseMatrix getOmega(
      const Scalings& scalings,
      boost::mpl::int_<static_cast<int>(lp::SolveFor::Initial)>) const {
    // TODO For initial, G Transpose is done twice unnecesarly
    // but I think its OK, as it is only one time, for the convinience of code
    return problem.G;
  }

  // W{-T} * G
  SparseMatrix getOmega(
      const Scalings& scalings,
      boost::mpl::int_<static_cast<int>(lp::SolveFor::StepDirection)>) const {
    return scalings.NNOInverse.asDiagonal() * problem.G;
  }

  // G' * W{-1} * W{-T}
  // For initial its just G'
  SparseMatrix getOmegaTilde(
      const Scalings& scalings,
      boost::mpl::int_<static_cast<int>(lp::SolveFor::Initial)>) const {
    return problem.G.transpose();
  }

  // G' * W{-1} * W{-T}
  // Transpose of diagonal matrix alters nothing
  SparseMatrix getOmegaTilde(
      const Scalings& scalings,
      boost::mpl::int_<static_cast<int>(lp::SolveFor::StepDirection)>) const {
    return omega.transpose() * scalings.NNOInverse.asDiagonal();
  }

  // Wz = W{-T} * (G*x - rhsz)
  Eigen::VectorXd computeWZ(
      const Scalings& scalings, const lp::NewtonDirection& direction,
      const Eigen::VectorXd& rhsZ,
      boost::mpl::int_<static_cast<int>(lp::SolveFor::StepDirection)>) const {
    return scalings.NNOInverse.asDiagonal() * (problem.G * direction.x - rhsZ);
  }

  // Wz = W{-T} * (G*x - rhsz)
  Eigen::VectorXd computeWZ(
      const Scalings& scalings, const lp::NewtonDirection& direction,
      const Eigen::VectorXd& rhsZ,
      boost::mpl::int_<static_cast<int>(lp::SolveFor::Initial)>) const {
    return problem.G * direction.x - rhsZ;
  }

  /**
   *  Logic is given in paper http://amestoy.perso.enseeiht.fr/doc/adru.pdf
   *
   */
  SparseMatrix scaleSymmetricMatrix(const SparseMatrix& symmetricMatrix) {

    SparseMatrix dm(symmetricMatrix.rows(), symmetricMatrix.rows());
    dm.setIdentity();

    SparseMatrix MHat = symmetricMatrix;

    while (true) {
      Eigen::VectorXd scalingVector = getScalingDiagonalVector(MHat);

      if (isScalingConverged(scalingVector)) {
        break;
      }

      dm = dm * scalingVector.asDiagonal();

      MHat = dm * symmetricMatrix * dm;
    }

    diagonalMatrix = dm;
    return MHat;
  }

  bool isScalingConverged(const Eigen::VectorXd& scalingVector) {

    double infNorm = std::pow(scalingVector.minCoeff(), 2);

    if ((1 - infNorm) <= problem.scalingTolerance) {
      return true;
    }

    return false;
  }

  Eigen::VectorXd getScalingDiagonalVector(
      const SparseMatrix& symmetricMatrix) {
    Eigen::VectorXd diagonalVector(symmetricMatrix.rows());

    Eigen::SparseVector<double> column;

    for (int i = 0; i < symmetricMatrix.cols(); ++i) {
      column = symmetricMatrix.col(i);

      double maxElem = 0;
      // Find infinity norm in a vector and assign it to diagonalVector
      column = column.unaryViewExpr([&](double coeff)->double {
        maxElem = std::max<double>(maxElem, std::abs(coeff));
        return coeff;
      });

      if (maxElem == 0) {
        diagonalVector(i) = 1;
      } else {
        diagonalVector(i) = 1 / std::sqrt(maxElem);
      }
    }

    return diagonalVector;
  }

  const lp::Problem& problem;
  CholeskySolver solver;
  // Custom Solver
  //   Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double, Eigen::Lower>>
  // solver;
  //   Eigen::PastixLDLT<Eigen::SparseMatrix<double>, Eigen::Lower> solver;
  // All below variables are used to store intermediate results to avoid
  // repeated calculations or copies (Product computaions are coslty)

  // G' * W{-1} * W{-T} Should be used in multiple solves to solve
  SparseMatrix omegaTilde;  // rhsZCoefficient;

  // G' * W{-1} Lets call it Omega
  SparseMatrix omega;
  // Used in scaling matrix for numerical stability
  SparseMatrix diagonalMatrix;
};

/**
 * Extension of Eigen cholmod wrapper to do symbolic factorize only once at the
 * start solver of iteration
 *
 * AnalyzePatter should be called by the user, before calling compute
 *
 */
template <typename _MatrixType, int _UpLo = Eigen::Lower>
class CholmodSupernodalLLT
    : public Eigen::CholmodBase<_MatrixType, _UpLo,
                                CholmodSupernodalLLT<_MatrixType, _UpLo>> {

  typedef Eigen::CholmodBase<_MatrixType, _UpLo, CholmodSupernodalLLT> Base;

  using Base::m_cholmod;

 public:
  typedef _MatrixType MatrixType;

  CholmodSupernodalLLT() : Base() { init(); }

  void compute(const MatrixType& matrix) { Base::factorize(matrix); }

  ~CholmodSupernodalLLT() {}
  /*
    void analyzePattern(const MatrixType& matrix) {
      if (this->m_cholmodFactor) {
        cholmod_free_factor(&this->m_cholmodFactor, &this->m_cholmod);
        this->m_cholmodFactor = 0;
      }
      cholmod_sparse A = viewAsCholmod(matrix);
      A.stype = -1;
      this->m_cholmodFactor = cholmod_analyze(&A, &m_cholmod);

      this->m_isInitialized = true;
      this->m_info = Eigen::Success;
      this->m_analysisIsOk = true;
      this->m_factorizationIsOk = false;
    }

    void factorize(const MatrixType& matrix) {
      eigen_assert(Base::m_analysisIsOk && "You must first call
    analyzePattern()");
      cholmod_sparse A = viewAsCholmod(matrix);
      A.stype = -1;
      cholmod_factorize_p(&A, this->m_shiftOffset, 0, 0, this->m_cholmodFactor,
    &m_cholmod);

      // If the factorization failed, minor is the column at which it did. On
      // success minor == n.
      this->m_info =
          (this->m_cholmodFactor->minor == this->m_cholmodFactor->n ?
    Eigen::Success
                                                        :
    Eigen::NumericalIssue);
      this->m_factorizationIsOk = true;
    }*/

 protected:
  void init() {
    m_cholmod.final_asis = 1;
    m_cholmod.supernodal = CHOLMOD_SUPERNODAL;
    // FIXME Set to metis and test
    m_cholmod.nmethods = 2;
  }
};
}

#endif  // SUITESPARSE_CHOLESKYLLT_HPP