#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <iostream>

#include <Eigen/Dense>

#include <Core/Residuals.hpp>
#include <Core/Point.hpp>
#include <Problem.hpp>

namespace lp {

/**
 * Contains Residual and final Point (Point can be made internal as point is
 * copied(moved) to solution)
 * FIXME Add more information as part solution which are in residuals if
 * necessary
 */
class Solution {

 public:
  Solution(const Residual& residual, const Point& point)
      : x(getFinalSolution(point.x, point.tau)),
        s(getFinalSolution(point.s, point.tau)),
        y(getFinalSolution(point.y, point.tau)),
        z(getFinalSolution(point.z, point.tau)),
        gap(residual.gap),
        primalObjective(residual.primalObjective),
        dualObjective(residual.dualObjective),
        iterations(residual.iterations) {}

  const Eigen::VectorXd x;
  const Eigen::VectorXd s;
  const Eigen::VectorXd y;
  const Eigen::VectorXd z;

  const double gap;
  const double primalObjective;
  const double dualObjective;
  const int iterations;

 private:
  Eigen::VectorXd getFinalSolution(const Eigen::VectorXd& vector,
                                   const double tau) {
    return vector.unaryExpr([&](const double coeff) { return coeff / tau; });
  }
};

std::ostream& operator<<(std::ostream& out, const Solution& solution) {
  using namespace std;

  out << endl << "##################### Final Solution  " << endl;
  out << "Value of x: " << endl << solution.x << endl;
  out << "Value of s: " << endl << solution.s << endl;
  out << "Value of y: " << endl << solution.y << endl;
  out << "Value of z: " << endl << solution.z << endl;
  out << "-----------------------------------" << endl;
  out << "Value of dual gap: " << solution.gap << endl;
  out << "Value of Primal Objective: " << solution.primalObjective << endl;
  out << "Value of Dual Objective: " << solution.dualObjective << endl;
  out << "Number of Iterations: " << solution.iterations << endl;

  return out;
}
}

#endif  // SOLUTION_HPP
