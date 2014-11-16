#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <algorithm>

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>

#include <Core/Point.hpp>
#include <Core/Residuals.hpp>
#include <Problem.hpp>
#include <Solution.hpp>
// FIXME If there is chance of using mathematical symbols in text editor use it
// Copy symbols from http://math.typeit.org/
namespace lp {

/**
 *Solves Convex problems
 *Linear solver is provided as template argument
 */
template <typename LinearSolver, typename Scalings>
class Solver {
 public:
  Solver(const Problem& problem) : problem(problem), linearSolver(problem) {}

  /**
   *
   */
  void solve() {

    const double residualX0{std::max(1.0, problem.c.norm())};
    const double residualY0{std::max(1.0, problem.b.norm())};
    const double residualZ0{std::max(1.0, problem.h.norm())};

    Point point = getInitialPoint();

    BOOST_LOG_TRIVIAL(info) << "Initial Points";
    BOOST_LOG_TRIVIAL(info) << point;

    Residual tolerantResidual(problem);

    // TODO Either start from 1 or use just lessthan instead of lessthan or eual
    // to operator, or maybe I can use these indications to find start/end of
    // loop
    for (int i = 0; i <= problem.maxIterations; ++i) {

      Residual residual(problem, point, i, residualX0, residualY0, residualZ0);
      SolverState solverState = getSolverState(residual, tolerantResidual);

      BOOST_LOG_TRIVIAL(info) << "Iteration: " << i;
      // TODO Build solution object to return
      if (solverState == SolverState::Feasible) {
        BOOST_LOG_TRIVIAL(info) << "Solution found ";
        break;  // Maybe return
      } else if (solverState == SolverState::Infeasible) {
        BOOST_LOG_TRIVIAL(info) << "Solution found, its Infeasible ";
        break;  // Maybe return
      }

      // Compute scalings for given point
      Scalings scalings(point);

      BOOST_LOG_TRIVIAL(info) << scalings;

      NewtonDirection subSolution = getSubSolution(scalings);

      BOOST_LOG_TRIVIAL(info) << subSolution;

      break;
    }
  }

 private:
  const Problem& problem;
  LinearSolver linearSolver;

  // As name implies its dummy object as Initial point does not need any
  // scaling point
  // TODO Check if creating and passing dummy object can be problamatic for
  // perforamnce
  // TODO Check what is the cost of creating dummy object, if it is high had
  // to find another way
  // TODO May be not much of a problem as it is used only once
  // TODO As Interface dicatates I have to send Identity Matrix, as I know the
  // implementation inside I have taken shortcut but its wrong
  Point getInitialPoint() {

    Point point(problem);
    Scalings scalings;

    linearSolver.template factor<lp::SolveFor::Initial>(scalings);

    {
      // Get Primal initial points
      Eigen::VectorXd rhsX(problem.c.rows());
      rhsX.setZero();

      NewtonDirection direction =
          linearSolver.template solve<lp::SolveFor::Initial>(
              rhsX, problem.b, problem.h, scalings);
      // TODO FIXME Can we invoke move constructor instead of copying data?
      point.x = direction.x;
      point.s = -1 * direction.z;

      adjustInitialPoint(point.s);
    }
    {
      // Get Dual initial points
      Eigen::VectorXd rhsX(-1 * problem.c);
      Eigen::VectorXd rhsY(problem.A.rows());
      rhsY.setZero();
      Eigen::VectorXd rhsZ(problem.G.rows());
      rhsZ.setZero();

      NewtonDirection direction =
          linearSolver.template solve<lp::SolveFor::Initial>(rhsX, rhsY, rhsZ,
                                                             scalings);
      // FIXME Move constructor?
      point.y = direction.y;
      point.z = direction.z;

      adjustInitialPoint(point.z);
    }
    {
      point.kappa = 1;
      point.tau = 1;
    }
    return point;
  }

  /**
   * vector = vector | alpha < 0
   * vector = vector + (1 + alpha) | otherwise
   *
   */
  void adjustInitialPoint(Eigen::VectorXd& vector) {
    double alpha = getMaxStep(vector);

    if (alpha >= -1e-8 * std::max(vector.norm(), 1.0)) {
      vector = (vector.array() + 1 + alpha).matrix();
    }
  }

  /**
   * Finds solution to sub KKT system
   *
   * [0   A'   G'] [x]    [-c]
   * [A   0    0 ] [y] =  [b]
   * [G   0   -W'W][z]	  [h]
   */
  NewtonDirection getSubSolution(const Scalings& scalings) {
    // Only one factorization in loop is done here
    linearSolver.template factor<lp::SolveFor::StepDirection>(scalings);

    Eigen::VectorXd rhsX(-1 * problem.c);

    return linearSolver.template solve<lp::SolveFor::StepDirection>(
        rhsX, problem.b, problem.h, scalings);
  }

  /**
   *
   * o\ is inverse operation for o (diamond is used in cvxopt document)
   * return value of z is scaled, actual value is W * z
   * FIXME Why Scaled z
   *
   * [0   A'   G'] [x]    [dx]
   * [A   0    0 ] [y] =  [-dy]
   * [G   0   -W'W][z]	  [W{T}(lambda o\ ds) - dz]
   *
   * Where
   *
   */
  Point getAffineDirection(const Scalings scalings, const Residual& residual,
                           const NewtonDirection& subSolution) {}
};

}  // lp

#endif  // SOLVER_HPP
