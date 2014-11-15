#ifndef POINT_HPP
#define POINT_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Problem.hpp>

namespace lp {

// TODO Where should I place this object? In what namespace internal or lp, does
// solution header make sense? If not here what is interface for the user to get
// solution

// TODO First make all variables private and Make Point internal and friend of
// Residual and Solution object

class Point {
 public:
  // Defines sizes, This is default
  Point(const Problem& problem)
      : x(problem.c.rows()),
        s(problem.G.rows()),
        y(problem.A.rows()),
        z(problem.G.rows()) {}

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
}
#endif  // POINT_HPP