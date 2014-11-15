#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <Eigen/Dense>

#include <Problem.hpp>

namespace lp {

// Solver State
enum class SolverState {
  Feasible,
  Infeasible,
  InProgress
};

/**
 * Contains Residual and final Point (Point can be made internal as point is
 * copied(moved) to solution)
 */
class Solution {};
}

#endif  // SOLUTION_HPP
