#ifndef SOLVER_HPP
#define SOLVER_HPP

// TODO What shall we name this type?
//#define LinearSolver class

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>

#include <Core/Point.hpp>
#include <Core/Residuals.hpp>
#include <Problem.hpp>
#include <Solution.hpp>

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

    Point point = getInitialPoint();

    BOOST_LOG_TRIVIAL(info) << "Initial Points";
    BOOST_LOG_TRIVIAL(info) << point;

    Residual tolerantResidual(problem);

    // TODO Either start from 1 or use just lessthan instead of lessthan or eual
    // to operator, or maybe I can use these indications to find start/end of
    // loop
    for (int i = 0; i <= problem.maxIterations; ++i) {

      Residual residual(problem, point, i);
      SolverState solverState = getSolverState(residual, tolerantResidual);

      // TODO Build solution object to return
      if(solverState == SolverState::Feasible){
	BOOST_LOG_TRIVIAL(info) << "Solution found ";
	break; // Maybe return
      } else if(solverState == SolverState::Infeasible) {
	BOOST_LOG_TRIVIAL(info) << "Solution found, its Infeasible ";
	break; // Maybe return
      }




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
    }
    {
      point.kappa = 1;
      point.tau = 1;
    }
    return point;
  }

};

}  // lp

#endif  // SOLVER_HPP
