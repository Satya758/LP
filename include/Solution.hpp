#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <iostream>
#include <climits>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Problem.hpp>

namespace lp {
namespace internal {}

// TODO Where should I place this object? In what namespace internal or lp, does
// solution header make sense? If not here what is interface for the user to get
// solution

// TODO Make Point internal and friend of Residual and Solution object
class Point {
 public:
  // Defines sizes, This is default
  Point(const Problem& problem)
      : x(problem.c.rows()),
        s(problem.G.rows()),
        y(problem.A.rows()),
        z(problem.G.rows()) {}

  // Find Initial Point, it needs LinearSolver
  template <typename LinearSolver>
  Point(const Problem& problem, LinearSolver& linearSolver)
      : x(problem.c.rows()),
        s(problem.G.rows()),
        y(problem.A.rows()),
        z(problem.G.rows()) {
    // As name implies its dummy object as Initial point does not need any
    // scaling point
    // TODO Check if creating and passing dummy object can be problamatic for
    // perforamnce
    // TODO Check what is the cost of creating dummy object, if it is high had
    // to find another way
    // TODO May be not much of a problem as it is used only once
    // TODO As Interface dicatates I have to send Identity Matrix, as I know the
    // implementation inside I have taken shortcut but its wrong
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> dummy;

    linearSolver.template factor<lp::SolveFor::Initial>(dummy, dummy);

    {
      // Get Primal initial points
      Eigen::VectorXd rhsX(problem.c.rows());
      rhsX.setZero();

      NewtonDirection direction =
          linearSolver.template solve<lp::SolveFor::Initial>(rhsX, problem.b,
                                                             problem.h, dummy);
      // TODO FIXME Can we invoke move constructor instead of copying data?
      this->x = direction.x;
      this->s = -1 * direction.z;
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
                                                             dummy);
      // FIXME Move constructor?
      this->y = direction.y;
      this->z = direction.z;
    }
    {
      this->kappa = 1;
      this->tau = 1;
    }
  }

  // I dont have to declare this as friend as all member variables are public
  // anyway
  // friend std::ostream& operator<< (std::ostream &out, const Point& point);

  // Primal Variables
  Eigen::VectorXd x;
  Eigen::VectorXd s;
  // Dual Variables
  Eigen::VectorXd y;
  Eigen::VectorXd z;
  // Homogenizing variables
  double kappa;
  double tau;
};

// TODO Point members are not constant be carefull
std::ostream& operator<<(std::ostream& out, const Point& point) {
  using namespace std;

  out << endl << "##################### Point Start" << endl;
  out << "Primal Variable x:" << endl << point.x << endl;
  out << "Primal Variable s:" << endl << point.s << endl;
  out << "Dual Variable y:" << endl << point.y << endl;
  out << "Dual Variable z:" << endl << point.z << endl;
  out << "Homogenizing Variable kappa: " << endl << point.kappa << endl;
  out << "Homogenizing Variable tau: " << endl << point.tau << endl;
  out << "##################### Point End" << endl;
}

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
// TODO Shall we move this to Solver header and make it internal
enum class SolverState {
  Feasible,
  Infeasible,
  InProgress
};

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
/**
 * Contains Residual and final Point (Point can be made internal as point is copied(moved) to solution)
 */
class Solution {};
}

#endif  // SOLUTION_HPP
