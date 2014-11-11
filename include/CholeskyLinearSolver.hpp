#ifndef CHOLESKY_LINEAR_SOLVER_HPP
#define CHOLESKY_LINEAR_SOLVER_HPP

#include <type_traits>

#include <boost/utility.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Problem.hpp>
// TODO Using namespace lp everywhere temporarly
// maybe use better name
namespace lp {

/**
 *All of the IPM related work is done in this base class
 *Users who extend this Base object should solve matrices
 *provided to them
 *
 *Implement following methods
 *TODO
 *
 *Q. What this class does
 *A. Calculate all computations related to factoriazation and solving rhs which
 *involves IPM alogithm
 *
 *
 */
template <class Derived, class Factor>
class CholeskyLinearSolver : boost::noncopyable {

 public:
  CholeskyLiearSolver(Problem& problem) : problem(problem) {}

  // Method is final as I do not know at present this would not bring any good,
  // I have to think about performance though
  // TODO Check if making method final is a good choice
  // Lets call scaling matrix as omega
  // Template can be only two things integer with value 1 or SparseMatrix
  //
  // Factor matrix G' * W{-1} * W{-T} * G
  // If scalingMatrix W is Identity factor G' * G (for initial values)
  // If A (equality constaints) is not empty then factor A * omega{-1} * A' to
  // find y
  // First find y, then x and finally z
  // Diagonal Scaling matrix is only good for Linear programming
  template <class ScalingMatrix>
  void factor(ScalingMatrix& omega) final {
    // This is more of documentation than error message
    // Tells what are supported type as scaling matrix
    static_assert(std::is_integral<ScalingMatrix>::value ||
                      std::is_class<ScalingMatrix>::value,
                  "Not Correct Scaling Matrix");

    // Find G' * W{-1} * W{-T} * G
    // TODO Singualar case is not considered (In next releases)
    firstFactor = child()->compute(problem.G.transpose() *
                                   getCOmega<ScalingMatrix>(omega) * problem.G);

    if(problem.A.size() != 0){
	Eigen::SparseMatrix<double> firstFactorAsEigen = child()->getFactorAsEigen(firstFactor);

    }
  }

  // TODO Check above method for reason for being final
  // SolveRhs why not just solve, as solve is left for dervied objects to
  // implement
  void solveRhs() final {}

 private:
  const Problem& problem;

  Factor const* firstFactor;
  Factor const* secondFactor;

  // helper
  Derived& child() { return *static_cast<Derived*>(this); }

  // Get W{-1} * W{-T} call it cOmega (Capital Omega)
  template <class ScalingMatrix>
  ScalingMatrix getCOmega(ScalingMatrix& omega) {}

  // For Initial return integer, its done this way to avoid duplication of
  // expression
  // TODO Check if its over engineering
  template <>
  int getCOmega<int>(int& omega) {
    return 1;
  }
};
}

#endif  // CHOLESKY_LINEAR_SOLVER_HPP