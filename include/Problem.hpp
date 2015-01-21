
#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <cstddef>

#include <Eigen/Sparse>
#include <Eigen/Dense>

// lp means lineat programming
// TODO Change to relavent namespace later
namespace lp {

// TODO Is this correct place to have solver typedef? It might make sense as
// users also use this header to represent their matrices
typedef Eigen::SparseMatrix<double, Eigen::ColMajor, std::ptrdiff_t>
    SparseMatrix;

// Used in LinearSolver
// Template specialization for solve method (instead of boolean flag)
// May be its too much but its clean
enum class SolveFor {
  Initial,
  StepDirection
};

// Solver State
enum class SolverState {
  Feasible,
  PrimalInfeasible,
  DualInfeasible,
  InProgress,
  MaximumIterations
};

/**
 *Defines the problem and options to solver
 *
 *Solves primal problem
 * 	minimize     c'x
 *    subject to   Gx + s = h or Gx <= h (Inequality constraints)
 * 		   Given Inequality constraints should be of the form Gx <= h
 * 		   (s is added to simplify algorithm)
 * 		   Ax = b (Equality constraints)
 *		   s >= 0, s is subset of Cone
 * 	Dimensions c = n x 1
 * 		   x = n x 1
 * 		   G = m x n
 *		   h = m x 1 (Slack variable s has same dimension as h)
 *		   s = m x 1 (Slack variable is not passed from user)
 *		   A = me x n (me is m rows for equality constraints)
 * 		   b = me x 1
 *
 *And its dual problem
 * 	maximize    -h'z - b'y
 * 	subject to  G'z + A'y + c = 0
 *		    z >= 0, z is subset of Cone
 *
 *Cone is cartesian product of muliple cones ( Postive orthant, Second order
 *cone, SDP)
 *XXX Currently solver supports only for positive orthant cone i.e. Linear
 *programming
 */
class Problem {
 public:
  // Problem constructor forces to enter dimensions to have valid entries for
  // all related matrices

  Problem(int inequalityRows, int equalityRows, int cols,
          int ipmMaxIterations_ = 200, double ipmTolerance_ = 1e-7,
          double scalingTolerance_ = 1e-3, int admmMaxIterations_ = 25000,
          double admmTolerance_ = 1e-3, double admmRho_ = 1e-3,
          double admmRelaxationParameter_ = 1.5)
      : G(inequalityRows, cols),
        h(inequalityRows),
        c(cols),
        ipmMaxIterations(ipmMaxIterations_),
        ipmTolerance(ipmTolerance_),
        scalingTolerance(scalingTolerance_),
        admmMaxIterations(admmMaxIterations_),
        admmTolerance(admmTolerance_),
        admmRho(admmRho_),
        admmRelaxationParameter(admmRelaxationParameter_) {}
  // Objective to minimize
  Eigen::VectorXd c;
  // Inequality constraints
  SparseMatrix G;
  Eigen::VectorXd h;

  // Maximum number of iterations after which algorithm terminates if feasible
  // solution not found
  const int ipmMaxIterations;
  // Duality gap relative to primal or dual objective
  // which is gap/pobj or gap/dobj
  const double ipmTolerance;
  // To attain numerical stability // TODO Is this used?
  const double scalingTolerance;
  // ADMM Options
  const int admmMaxIterations;
  const double admmTolerance;
  // TODO Find More about this option
  const double admmRho;  // Scaling parameter
  // Relaxation parameter alpha
  const double admmRelaxationParameter;

  // TODO Some more Options

  /**
   * iterationPercentage is number of iterations performed to find Diagonal
   * scalings for given matrix, percentage of size of matrix, where size =
   * max(m, n)
   *
   */
  void normalize();
};

std::ostream& operator<<(std::ostream& out, const Problem& problem);

// TODO Check if we can move this class to different header
// Solution after newton step computed
// Helper class to transfer NewtonDirection to IPM algorithm
class NewtonDirection {
 public:
  // Variable names of solution are given in generic names such as x, y, z
  // As solution is for 3 X 3 block matrix we have three dense vectors denoting
  // that
  Eigen::VectorXd x;
  Eigen::VectorXd z;
};

NewtonDirection operator*(const double lhs, const NewtonDirection& rhs);

std::ostream& operator<<(std::ostream& out, const NewtonDirection& direction);
}
#endif  // PROBLEM_HPP