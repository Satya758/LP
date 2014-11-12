#ifndef SOLVER_HPP
#define SOLVER_HPP

// TODO What shall we name this type?
//#define LinearSolver class

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>

#include <Problem.hpp>

namespace lp {
namespace internal {

/**
 *Solves Convex problems
 *Linear solver is provided as template argument
 */
template <class LinearSolver>
class Solver {
 public:
  Solver(const Problem& problem) : problem(problem), linearSolver(problem) {

    Eigen::DiagonalMatrix<double, Eigen::Dynamic> dummy;

    // Initial values
    linearSolver.factor(dummy, dummy, lp::SolveFor::Initial);

    Eigen::VectorXd rhsX(problem.c.rows());

    rhsX.setZero();

    linearSolver.solve(rhsX, problem.b, problem.h, dummy,
                       lp::SolveFor::Initial);
  }

 private:
  const Problem& problem;
  LinearSolver linearSolver;
};

}  // internal
}  // lp

#endif  // SOLVER_HPP
