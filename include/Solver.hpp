#ifndef SOLVER_HPP
#define SOLVER_HPP

// TODO What shall we name this type?
//#define LinearSolver class

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>

#include <Problem.hpp>

#include <Solution.hpp>

namespace lp {
namespace internal {

/**
 *Solves Convex problems
 *Linear solver is provided as template argument
 */
template <typename LinearSolver>
class Solver {
 public:
  Solver(const Problem& problem) : problem(problem), linearSolver(problem) {}

  /**
   *
   */
  void solve() {
    Point initialPoint(this->problem, linearSolver);

    BOOST_LOG_TRIVIAL(info) << initialPoint;
  }

 private:
  const Problem& problem;
  LinearSolver linearSolver;
};

}  // internal
}  // lp

#endif  // SOLVER_HPP
