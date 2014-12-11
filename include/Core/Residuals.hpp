#ifndef RESIDUALS_HPP
#define RESIDUALS_HPP

#include <iostream>
#include <climits>
#include <algorithm>

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Core/Point.hpp>
#include <Problem.hpp>

namespace lp {

/**
 * find alpha = inf{alpha| vector + alpha * e >= 0 }
 * e is vector of ones of same dimension as vector
 * >= is defined for each element in vector
 * alpha can be computed by finding minimum element, which when added to vector
 *land at zero
 *
 * So max step that can be taken before breaching condition
 *
 * For SOCP and SDP there is much more to do (Only works for LP) as vector is
 *normalized in SOCP and SDP
 */
double getMaxStep(const Eigen::VectorXd& vector) { return -vector.minCoeff(); }

/**
 * Calculates residual, tolerant values
 */
class Residual {
  // Declared in the same order as required initilization, see constructor
  // initilization
  // Here residualX... and subResidualX... are required before computing public
  // variables
 private:
  double residualX0;
  double residualZ0;

  // s + G*x
  Eigen::VectorXd subResidualZ;
  // -A'*y - G'*z
  Eigen::VectorXd subResidualX;

 public:
  // Compute residual for given point
  Residual(const Problem& problem, const Point& point, int iterations,
           const double residualX0, const double residualZ0)
      : residualX0(residualX0),
        residualZ0(residualZ0),
        subResidualZ(point.s + problem.G * point.x),
        subResidualX(-problem.G.transpose() * point.z),
        residualX(getResidualX(problem, point)),
        residualZ(getResidualZ(problem, point)),
        residualTau(getResidualTau(problem, point)),
        gap(point.s.dot(point.z)),
        primalObjective(problem.c.dot(point.x) / point.tau),
        dualObjective(-(problem.h.dot(point.z)) / point.tau),
        relativeGap(getRelativeGap()),
        primalResidual(getPrimalResidual(problem, point)),
        dualResidual(getDualResidual(problem, point)),
        primalInfeasibility(getPrimalInfeasibility(problem, point)),
        dualInfeasibility(getDualInfeasibility(problem, point)),
        iterations(iterations),
        primalSlack(getMaxStep(point.s)),
        dualSlack(getMaxStep(point.z)) {}

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
        dualObjective(0),
        // Initializing irrevelant vectors in this scenario with least int to
        // avoid performance hit
        residualX(1),
        residualZ(1),
        residualTau(1) {}

  const Eigen::VectorXd residualX;
  const Eigen::VectorXd residualZ;
  const double residualTau;

  const double gap;

  const double primalObjective;
  const double dualObjective;

  const double relativeGap;

  const double primalResidual;
  const double dualResidual;

  const double primalInfeasibility;
  const double dualInfeasibility;

  const double primalSlack;
  const double dualSlack;

  const int iterations;

 private:
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

  // -G'*z - c*tau
  Eigen::VectorXd getResidualX(const Problem& problem, const Point& point) {
    return subResidualX - problem.c * point.tau;
  }

  // s + G*x - h*tau
  Eigen::VectorXd getResidualZ(const Problem& problem, const Point& point) {
    return subResidualZ - problem.h * point.tau;
  }
  // kappa + c'*x + h'*z
  double getResidualTau(const Problem& problem, const Point& point) {
    // FIXME c'*x and h'*z are computed again for primalObjective and
    // dualObjective
    return point.kappa + problem.c.dot(point.x) + problem.h.dot(point.z);
  }

  double getPrimalResidual(const Problem& problem, const Point& point) const {

    // s + G*x - h*tau
    double resZ = (residualZ.norm()) / point.tau;

    return resZ / residualZ0;
  }

  double getDualResidual(const Problem& problem, const Point& point) const {
    // -A'*y - G'*z - c*tau
    double resX = (residualX.norm()) / point.tau;

    return resX / residualX0;
  }

  double getPrimalInfeasibility(const Problem& problem,
                                const Point& point) const {
    // hz + by
    // TODO Check, tau is muliplied to dualObjective
    if (dualObjective > 0) {
      return (subResidualX.norm() * dualObjective) / residualX0;
    } else {
      // TODO Some large number, is it correct way!
      return 100;
    }
  }

  double getDualInfeasibility(const Problem& problem,
                              const Point& point) const {
    // cx
    // TODO Check, tau is muliplied to primalObjective
    if (primalObjective < 0) {
      return (subResidualZ.norm() / residualZ0) / std::abs(primalObjective);
    } else {
      // TODO Some large number, is it correct way!
      return 100;
    }
  }
};

// std::ostream& operator<<(std::ostream& out, const Residual& residual) {
//   using namespace std;
//
//   out << endl << "##################### Residual Start" << endl;
//   out << "Residual X : " << endl << residual.residualX << endl;
//   out << "Residual Y : " << endl << residual.residualY << endl;
//   out << "Residual Z : " << endl << residual.residualZ << endl;
//   out << "Residual Tau : " << endl << residual.residualTau << endl;
//   out << "Duality Gap: " << residual.gap << endl;
//   out << "Relative Gap: " << residual.relativeGap << endl;
//   out << "Primal Slack: " << residual.primalSlack << endl;
//   out << "Dual Slack: " << residual.dualSlack << endl;
//   out << "Primal Objective: " << residual.primalObjective << endl;
//   out << "Dual Objective: " << residual.dualObjective << endl;
//   out << "Primal Residual: " << residual.primalResidual << endl;
//   out << "Dual Residual: " << residual.dualResidual << endl;
//   out << "Primal Infeasibility: " << residual.primalInfeasibility << endl;
//   out << "Dual Infeasibility: " << residual.dualInfeasibility << endl;
//   out << "Iterations: " << residual.iterations << endl;
//   out << "##################### Residual End" << endl;
//
//   return out;
// }

std::ostream& operator<<(std::ostream& out, const Residual& residual) {
  using namespace std;

  if (residual.iterations == 0) {
    out << "Iteration \t\t"
        << "PCost \t\t"
        << "DCost \t\t"
        << "Gap \t\t"
        << "PResidual \t\t"
        << "Dresidual \t\t" << endl;
  }

  out << residual.iterations << "\t\t" << residual.primalObjective << "\t\t"
      << residual.dualObjective << "\t\t" << residual.gap << "\t\t"
      << residual.primalResidual << "\t\t" << residual.dualResidual << "\t\t";

  return out;
}

// TODO Only One operator is implmented and only make sense to use it in end
// scenarios
// Very restrictive use case
bool operator<=(const Residual& lhs, const Residual& rhs) {
  return (lhs.primalResidual <= rhs.primalResidual &&
          lhs.dualResidual <= rhs.dualResidual &&
          (lhs.gap <= rhs.gap || lhs.relativeGap <= rhs.relativeGap));
  // TODO Had to consider maximumIterations reached
  //|| lhs.iterations <= rhs.iterations);
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
             tolerantResidual.primalInfeasibility) {
    return SolverState::PrimalInfeasible;
  } else if (residual.dualInfeasibility <= tolerantResidual.dualInfeasibility) {
    return SolverState::DualInfeasible;
  } else {
    return SolverState::InProgress;
  }
}

// FIXME Remove this method after confirming that i do not need different way to
// initial points
SolverState getSolverStateForInitialPoint(const Residual& residual,
                                          const Residual& tolerantResidual) {
  if (residual.primalSlack <= 0 && residual.dualSlack <= 0 &&
      (residual.gap <= tolerantResidual.gap ||
       residual.relativeGap <= tolerantResidual.relativeGap)) {
    return SolverState::Feasible;
  } else {
    return SolverState::InProgress;
  }
}
}
#endif  // RESIDUALS_HPP