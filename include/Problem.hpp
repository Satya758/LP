
#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>

// lp means lineat programming
// TODO Change to relavent namespace later
namespace lp {

// Template specialization for solve method (instead of boolean flag)
// May be its too much
enum SolveFor {
  Initial,
  StepDirection
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
          int maxIterations = 100, double relativeGapTolerance = 1e-6,
          double gapTolerance = 1e-7, double residualTolerance = 1e-7)
      : G(inequalityRows, cols),
        h(inequalityRows),
        A(equalityRows, cols),
        b(equalityRows),
        c(cols),
        maxIterations(maxIterations),
        relativeGapTolerance(relativeGapTolerance),
        gapTolerance(gapTolerance),
        residualTolerance(residualTolerance) {}
  // Objective to minimize
  Eigen::VectorXd c;
  // Inequality constraints
  Eigen::SparseMatrix<double> G;
  Eigen::VectorXd h;
  // Equality constraints
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;

  // Maximum number of iterations after which algorithm terminates if feasible
  // solution not found
  const int maxIterations;
  // Duality gap relative to primal or dual objective
  // which is gap/pobj or gap/dobj
  const double relativeGapTolerance;
  // Duality gap
  const double gapTolerance;
  // Residual in primal and dual variables after newton step
  const double residualTolerance;

  // TODO Some more Options
};

// TODO Check if we can move this class to different header
// Solution after newton step computed
// Helper class to transfer NewtonDirection to IPM algorithm
class NewtonDirection {
 public:
  // Variable names of solution are given in generic names such as x, y, z
  // As solution is for 3 X 3 block matrix we have three dense vectors denoting
  // that
  // Case when equality constraints are not present then y is not used
  // computations
  // TODO Check when y is not used
  Eigen::VectorXd x;
  Eigen::VectorXd y;
  Eigen::VectorXd z;
};
}
#endif  // PROBLEM_HPP