#ifndef SUITESPARSE_CHOLESKYLLT_HPP
#define SUITESPARSE_CHOLESKYLLT_HPP

#include <boost/log/trivial.hpp>
#include <boost/mpl/int.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

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
template <typename Scalings>
class SuiteSparseCholeskyLLT {

 public:
  SuiteSparseCholeskyLLT(const lp::Problem& problem) : problem(problem) {
    BOOST_LOG_TRIVIAL(info) << "Constructor of Linear Solver";
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
    BOOST_LOG_TRIVIAL(info) << "Step directions Factorization";

    // First compute G' * W{-1} and store it
    omega = getOmega(scalings, boost::mpl::int_<solveFor>());

    // Used in solve method, cached as it is used multiple times
    omegaTilde = getOmegaTilde(scalings, boost::mpl::int_<solveFor>());

    solver.compute(omega * omega.transpose());

    if (solver.info() != Eigen::Success) {
      // TODO Raise exception maybe Matrix is singular
      BOOST_LOG_TRIVIAL(error) << "Factorization failed";
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
  template <lp::SolveFor solveFor>
  lp::NewtonDirection solve(const Eigen::VectorXd& rhsX,
                            const Eigen::VectorXd& rhsY,
                            const Eigen::VectorXd& rhsZ,
                            const Scalings& scalings) const {

    lp::NewtonDirection direction;

    {
      // rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
      Eigen::VectorXd rhs = rhsX + omegaTilde * rhsZ;

      direction.x = solver.solve(rhs);
    }
    {
      direction.z =
          computeWZ(scalings, direction, rhsZ, boost::mpl::int_<solveFor>());
    }

    return direction;
  }

 private:
  // G'
  Eigen::SparseMatrix<double> getOmega(
      const Scalings& scalings, boost::mpl::int_<lp::SolveFor::Initial>) const {
    // TODO For initial, G Transpose is done twice unnecesarly
    // but I think its OK, as it is only one time, for the convinience of code
    return problem.G.transpose();
  }

  // G' * W{-1}
  Eigen::SparseMatrix<double> getOmega(
      const Scalings& scalings,
      boost::mpl::int_<lp::SolveFor::StepDirection>) const {
    return problem.G.transpose() * scalings.NNOInverse.asDiagonal();
  }

  // G' * W{-1} * W{-T}
  // For initial its just G'
  Eigen::SparseMatrix<double> getOmegaTilde(
      const Scalings& scalings, boost::mpl::int_<lp::SolveFor::Initial>) const {
    return problem.G.transpose();
  }

  // G' * W{-1} * W{-T}
  // Transpose of diagonal matrix alters nothing
  Eigen::SparseMatrix<double> getOmegaTilde(
      const Scalings& scalings,
      boost::mpl::int_<lp::SolveFor::StepDirection>) const {
    return omega * scalings.NNOInverse.asDiagonal();
  }

  // Wz = W{-T} * (G*x - rhsz)
  Eigen::VectorXd computeWZ(
      const Scalings& scalings, const lp::NewtonDirection& direction,
      const Eigen::VectorXd& rhsZ,
      boost::mpl::int_<lp::SolveFor::StepDirection>) const {
    return scalings.NNOInverse.asDiagonal() * problem.G * direction.x - rhsZ;
  }

  // Wz = W{-T} * (G*x - rhsz)
  Eigen::VectorXd computeWZ(const Scalings& scalings,
                            const lp::NewtonDirection& direction,
                            const Eigen::VectorXd& rhsZ,
                            boost::mpl::int_<lp::SolveFor::Initial>) const {
    return problem.G * direction.x - rhsZ;
  }

  const lp::Problem& problem;
  // Custom Solver
  Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double, Eigen::Lower>> solver;
  // All below variables are used to store intermediate results to avoid
  // repeated calculations or copies (Product computaions are coslty)

  // G' * W{-1} * W{-T} Should be used in multiple solves to solve
  Eigen::SparseMatrix<double> omegaTilde;  // rhsZCoefficient;

  // G' * W{-1} Lets call it Omega
  Eigen::SparseMatrix<double> omega;
};
}

#endif  // SUITESPARSE_CHOLESKYLLT_HPP