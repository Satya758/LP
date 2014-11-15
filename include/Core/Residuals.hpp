#ifndef RESIDUALS_HPP
#define RESIDUALS_HPP

#include <iostream>
#include <climits>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Core/Point.hpp>
#include <Problem.hpp>
#include <Solution.hpp>

namespace lp {

/**
 * Calculates residual, tolerant values
 */
class Residual {
 public:
  // Compute residual for given point
  Residual(const Problem& problem, const Point& point, int iterations)
      : residualX0(std::max(1.0, problem.c.norm())),
        residualY0(std::max(1.0, problem.b.norm())),
        residualZ0(std::max(1.0, problem.h.norm())),
        subResidualY(problem.A * point.x),
        subResidualZ(point.s + problem.G * point.x),
        subResidualX(-problem.A.transpose() * point.y -
                     problem.G.transpose() * point.z),
        gap(point.s.dot(point.z)),
        primalObjective(problem.c.dot(point.x) / point.tau),
        dualObjective(-(problem.b.dot(point.y) + problem.h.dot(point.z)) /
                      point.tau),
        relativeGap(getRelativeGap()),
        primalResidual(getPrimalResidual(problem, point)),
        dualResidual(getDualResidual(problem, point)),
        primalInfeasibility(getPrimalInfeasibility(problem, point)),
        dualInfeasibility(getDualInfeasibility(problem, point)),
        iterations(iterations),
        primalSlack(-point.s.minCoeff()),
        dualSlack(-point.z.minCoeff()) {}

  // Compute Accepted Tolerence residual
  Residual(const Problem& problem)
      : gap(problem.gapTolerance),
        relativeGap(problem.relativeGapTolerance),
        primalResidual(problem.residualTolerance),
        dualResidual(problem.residualTolerance),
        iterations(problem.maxIterations),
        primalSlack(0),
        dualSlack(0),
        primalInfeasibility(problem.residualTolerance),
        dualInfeasibility(problem.residualTolerance),
        primalObjective(0),
        dualObjective(0) {}

  const double gap;
  const double relativeGap;

  const double primalSlack;
  const double dualSlack;

  const double primalObjective;
  const double dualObjective;

  const double primalResidual;
  const double dualResidual;

  const double primalInfeasibility;
  const double dualInfeasibility;

  const int iterations;

 private:
  // A*x
  Eigen::VectorXd subResidualY;
  // s + G*x
  Eigen::VectorXd subResidualZ;
  // -A'*y - G'*z
  Eigen::VectorXd subResidualX;

  double residualX0;
  double residualY0;
  double residualZ0;

  double getRelativeGap() const {
    if (primalObjective < 0) {
      return gap / std::abs(primalObjective);
    } else if (dualObjective > 0) {
      return gap / std::abs(dualObjective);
    } else {
      // Some large number, I think its good enough
      return 100;
    }
  }

  double getPrimalResidual(const Problem& problem, const Point& point) const {
    //  A*x - b*tau
    double resY = ((subResidualY - problem.b * point.tau).norm()) / point.tau;
    // s + G*x - h*tau
    double resZ = ((subResidualZ - problem.h * point.tau).norm()) / point.tau;

    return std::max(resY / residualY0, resZ / residualZ0);
  }

  double getDualResidual(const Problem& problem, const Point& point) const {
    // -A'*y - G'*z - c*tau
    double resX = ((subResidualX - problem.c * point.tau).norm()) / point.tau;

    return resX / residualX0;
  }

  double getPrimalInfeasibility(const Problem& problem,
                                const Point& point) const {
    // hz + by
    // TODO Check, tau is muliplied to dualObjective
    if (dualObjective > 0) {
      return (subResidualX.norm() * dualObjective) / residualX0;
    } else {
      return 0;
    }
  }

  double getDualInfeasibility(const Problem& problem,
                              const Point& point) const {
    // cx
    // TODO Check, tau is muliplied to primalObjective
    if (primalObjective < 0) {
      return std::max(subResidualY.norm() / residualY0,
                      subResidualZ.norm() / residualZ0) /
             std::abs(primalObjective);
    } else {
      return 0;
    }
  }
};

std::ostream& operator<<(std::ostream& out, const Residual& residual) {
  using namespace std;

  out << endl << "##################### Residual Start" << endl;
  out << "Duality Gap: " << residual.gap << endl;
  out << "Relative Gap: " << residual.relativeGap << endl;
  out << "Primal Slack: " << residual.primalSlack << endl;
  out << "Dual Slack: " << residual.dualSlack << endl;
  out << "Primal Objective: " << residual.primalObjective << endl;
  out << "Dual Objective: " << residual.dualObjective << endl;
  out << "Primal Residual: " << residual.primalResidual << endl;
  out << "Dual Residual: " << residual.dualResidual << endl;
  out << "Primal Infeasibility: " << residual.primalInfeasibility << endl;
  out << "Dual Infeasibility: " << residual.dualInfeasibility << endl;
  out << "Iterations: " << residual.iterations << endl;
  out << "##################### Residual End" << endl;
}

// TODO Only One operator is implmented and only make sense to use it in end
// scenarios
// Very restrictive use case
bool operator<=(const Residual& lhs, const Residual& rhs) {
  return (lhs.primalResidual <= rhs.primalResidual &&
              lhs.dualResidual <= rhs.dualResidual &&
              (lhs.gap <= rhs.gap || lhs.relativeGap <= rhs.relativeGap) ||
          lhs.iterations <= rhs.iterations);
}

// FIXME When problem is in Infeasible, values (Point & Residual) are changed
// which is not done
// here, Its not soo important as of now
// TODO And also when problem is Feasible and optimal values of point should be
// divided by tau is not done here
// TODO I am planning to change the point when building final solution lets see
SolverState getSolverState(const Residual& residual,
                           const Residual& tolerantResidual) {
  if (residual <= tolerantResidual) {
    return SolverState::Feasible;
  } else if (residual.primalInfeasibility <=
                 tolerantResidual.primalInfeasibility ||
             residual.dualInfeasibility <= tolerantResidual.dualInfeasibility) {
    return SolverState::Infeasible;
  } else {
    return SolverState::InProgress;
  }
}
}
#endif  // RESIDUALS_HPP