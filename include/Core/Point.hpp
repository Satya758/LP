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
  // TODO To use as dummy object is there any better way
  // In operator overload
  Point() {}
  // Used to create only z and tau object in Scalings getAffineDirection and
  // getCombinedDirection methods
  // As these vectors are dense objects, this approach would save memory by not
  // creating zero vectors, TODO Had to check this hypothesis
  Point(int zRows) : z(zRows) {}
  // Defines sizes, This is default
  Point(const Problem& problem)
      : x(problem.c.rows()), s(problem.G.rows()), z(problem.G.rows()) {}

  // I dont have to declare this as friend as all member variables are public
  // anyway
  // friend std::ostream& operator<< (std::ostream &out, const Point& point);

  // Primal Variables
  Eigen::VectorXd x;
  Eigen::VectorXd s;
  // Dual Variables
  Eigen::VectorXd z;
  // Dual slack residual, used in ADMM
  Eigen::VectorXd r;
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
  out << "Dual Variable z:" << endl << point.z << endl;
  out << "Homogenizing Variable kappa: " << endl << point.kappa << endl;
  out << "Homogenizing Variable tau: " << endl << point.tau << endl;
  out << "##################### Point End" << endl;

  return out;
}

// TODO Till below conundrum is resolved do it old school
Point updatePoint(const Point& lhs, const Point& rhs, const double alpha) {
  Point result;

  result.x = lhs.x + alpha * rhs.x;
  result.s = lhs.s + alpha * rhs.s;
  result.z = lhs.z + alpha * rhs.z;
  result.tau = lhs.tau + alpha * rhs.tau;
  result.kappa = lhs.kappa + alpha * rhs.kappa;

  return result;
}

// FIXME These operators are not used as I think they generate temporary objects
// before final result is computed
// All I need is point + alpha * point
// Becuase the way I have implemented it would generate two temporary objects
// before final result is obtained
// TODO I have to do lazy evaluations if its not too much, after going through
// that pain advantage is code would look more nicer
// noalias from Eigen would not help I think but I have to relook again
Point operator+(const Point& lhs, const Point& rhs) {
  Point result;

  result.x = lhs.x + rhs.x;
  result.s = lhs.s + rhs.s;
  result.z = lhs.z + rhs.z;
  result.tau = lhs.tau + rhs.tau;
  result.kappa = lhs.kappa + rhs.kappa;

  return result;
}

// FIXME Explanation given above
Point operator*(const double lhs, const Point& rhs) {
  Point result;

  result.x = lhs * rhs.x;
  result.s = lhs * rhs.s;
  result.z = lhs * rhs.z;
  result.tau = lhs * rhs.tau;
  result.kappa = lhs * rhs.kappa;

  return result;
}
}
#endif  // POINT_HPP