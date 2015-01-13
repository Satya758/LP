#ifndef ADMM_SOLVER_HPP
#define ADMM_SOLVER_HPP

#include <Core/Residuals.hpp>
#include <Core/ADMMProjections.hpp>
#include <Core/Point.hpp>
#include <Core/Timer.hpp>

#include <Ext/ADMMCholeskyLDLT.hpp>

#include <Problem.hpp>
#include <Solution.hpp>

namespace lp {
class ADMMSolver {
 public:
  ADMMSolver(const Problem& problem)
      : problem(problem), timer(Timer::getADMMInstance()) {}

  Solution solve() {

    Residual tolerantResidual(problem);
    ADMMProjections<ADMMCholeskyLDLT> admmProjections(problem);

    const double residualX0{1.0 + problem.c.norm()};
    const double residualZ0{1.0 + problem.h.norm()};

    Point currentPoint = this->getInitialPoint();

    for (int i = 0; i < problem.admmMaxIterations; ++i) {

      timer.start("Residual Computation");
      Residual residual(problem, currentPoint, i, residualX0, residualZ0);
      timer.end("Residual Computation");

      BOOST_LOG_TRIVIAL(info) << residual;

      timer.start("SolverState Computation");
      SolverState solverState = getSolverState(residual, tolerantResidual);
      timer.end("SolverState Computation");

      if (solverState == SolverState::Feasible) {
        return Solution(residual, currentPoint, solverState);
      } else if (solverState == SolverState::PrimalInfeasible ||
                 solverState == SolverState::DualInfeasible) {
        BOOST_LOG_TRIVIAL(info) << "Problem is Infeasible";
        return Solution(residual, currentPoint, solverState);
      } else if (i == problem.admmMaxIterations) {
        return Solution(residual, currentPoint, SolverState::MaximumIterations);
      }

      // Start of ADMM algorithm
      Point subspaceProjectedPoint =
          admmProjections.doSubspaceProjection(currentPoint);

      Point coneProjectedPoint = admmProjections.doConeProjection(
          subspaceProjectedPoint, currentPoint);

      admmProjections.doDualUpdates(subspaceProjectedPoint, currentPoint,
                                    coneProjectedPoint);

      currentPoint = coneProjectedPoint;
    }
  }

 private:
  const Problem& problem;
  Timer& timer;

  Point getInitialPoint() {
    Point point(problem);

    point.x.setZero();
    point.z.setZero();

    point.r.setZero();
    point.s.setZero();

    point.kappa = point.tau =
        std::sqrt(problem.c.rows() + problem.h.rows() + 1);

    return point;
  }
};
}
#endif  // ADMM_SOLVER_HPP