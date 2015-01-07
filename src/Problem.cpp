
#include <Problem.hpp>

namespace lp {

std::ostream& operator<<(std::ostream& out, const Problem& problem) {
  using namespace std;

  out << endl << "##################### Problem Start" << endl;
  out << "Value of c: " << endl << problem.c.rows() << endl;
  out << "Value of G: " << endl << "Rows: " << problem.G.rows()
      << " Cols: " << problem.G.cols() << endl;
  out << "Value of h: " << endl << problem.h.rows() << endl;
  out << "##################### Problem End" << endl;

  return out;
}

NewtonDirection operator*(const double lhs, const NewtonDirection& rhs) {
  NewtonDirection result;

  result.x = lhs * rhs.x;
  result.z = lhs * rhs.z;

  return result;
}

std::ostream& operator<<(std::ostream& out, const NewtonDirection& direction) {
  using namespace std;

  out << endl << "##################### NewtonDirection Start" << endl;
  out << "Value of x: " << endl << direction.x << endl;
  out << "Value of z: " << endl << direction.z << endl;
  out << "##################### NewtonDirection End" << endl;

  return out;
}
}