#ifndef ADMM_SOLVER_HPP
#define ADMM_SOLVER_HPP

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Core/Point.hpp>
#include <Core/Residuals.hpp>
#include <Problem.hpp>
#include <Solution.hpp>

namespace lp {

template <typename CholeskySolver>
class ADMMSolver {
 public:
  ADMMSolver(const Problem& problem) : problem(problem) {}

  Solution solve() {

    const double residualX0{std::max(1.0, problem.c.norm())};
    const double residualZ0{std::max(1.0, problem.h.norm())};

    Point currentPoint = getInitialPoint();
    Residual tolerantResidual(problem);

    CholeskySolver solver(problem);

    for (int i = 0; i <= problem.maxIterations; ++i) {

      BOOST_LOG_TRIVIAL(trace) << currentPoint;

      Residual residual(problem, currentPoint, i, residualX0, residualZ0);
      SolverState solverState = getSolverState(residual, tolerantResidual);

      if (solverState == SolverState::Feasible) {
        return Solution(residual, currentPoint, solverState);
      } else if (solverState == SolverState::PrimalInfeasible ||
                 solverState == SolverState::DualInfeasible) {
        return Solution(residual, currentPoint, solverState);
      } else if (i == problem.maxIterations) {
        return Solution(residual, currentPoint, SolverState::MaximumIterations);
      }

      BOOST_LOG_TRIVIAL(info) << residual;

      // Projection onto Feasible space
      Point omegaPrimal = solver.doProjection(getOmega(currentPoint));
      relaxOmegaPrimal(omegaPrimal);

      BOOST_LOG_TRIVIAL(trace) << "After Space proj: " << omegaPrimal;
      // Conic projection
      // currentPoint has dual variables
      Point psiPrimal = getPsi(omegaPrimal, currentPoint);
      doConicProjection(psiPrimal);

      BOOST_LOG_TRIVIAL(trace) << "After Cone proj: " << psiPrimal;
      // Update dual variables
      updateDual(omegaPrimal, psiPrimal, currentPoint);

      currentPoint = omegaPrimal;
    }
  }

 private:
  const Problem& problem;
  // alpha in documentation of scs
  const double relaxationParameter = 1.5;

  // Omega = priaml + dual variables
  // TODO  Notice here there is no dual residual r as in doc
  // Primal varaible is not touched as r is not considered by me
  Point getOmega(const Point& point) const {
    Point omega;
    omega.x = point.x + point.r;
    omega.z = point.z + point.s;
    omega.tau = point.tau + point.kappa;

    return omega;
  }

  // Primal varaible x is not changed
  void relaxOmegaPrimal(Point& point) const {
    point.z =
        relaxationParameter * point.z + (1 - relaxationParameter) * point.z;
    point.tau =
        relaxationParameter * point.tau + (1 - relaxationParameter) * point.tau;
  }

  // omegaPrimal - dual Varaibles
  // TODO dual Residual r is considered by me
  Point getPsi(const Point& point, const Point& dualPoint) const {
    // Problem is passed for dimensions
    Point psi(problem);

    psi.x = point.x - dualPoint.r;
    psi.z = point.z - dualPoint.s;
    psi.tau = point.tau - dualPoint.kappa;

    return psi;
  }

  // Only Positive orthant is implemented for now
  // x is not considered here as well
  void doConicProjection(Point& point) const {

    point.z = point.z.unaryExpr([](const double & element)->double {
      if (element < 0) {
        return 0;
      } else {
        return element;
      }
    });

    if (point.tau < 0) {
      point.tau = 0;
    }
  }

  // v{k=1} = v{k} - uHat{k+1} + u{k+1}
  // Again r is ignored here
  void updateDual(Point& omegaPrimal, const Point& psi,
                  const Point& dualPoint) const {
    omegaPrimal.r = dualPoint.r;
    omegaPrimal.s = psi.z + (dualPoint.s - omegaPrimal.z);
    omegaPrimal.kappa = psi.tau + (dualPoint.kappa - omegaPrimal.tau);
  }

  Point getInitialPoint() {
    Point point;

    point.x = Eigen::VectorXd::Zero(problem.c.rows());
    point.z = Eigen::VectorXd::Zero(problem.h.rows());
    point.s = Eigen::VectorXd::Zero(problem.h.rows());
    point.r = Eigen::VectorXd::Zero(problem.c.rows());
    point.tau = point.kappa =
        std::sqrt(problem.c.rows() + problem.h.rows() + 1);

    return point;
  }
};
}
#endif  // ADMM_SOLVER_HPP