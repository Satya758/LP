#ifndef NT_SCALINGS_HPP
#define NT_SCALINGS_HPP

#include <Eigen/Dense>

#include <Core/Point.hpp>

namespace lp {
/**
 * Neterov-Todd (NT) scaling is used
 * TODO Solver can accepts scaling as template which can implement any scaling
 * problem demands in future, but for now let us hardcode scaling object in
 * Solver
 *
 * TODO Find home(header) and namespace for this class
 * TODO FIXME Check if using DiagonalWrapper reduces copy of vector into
 * DiagonalMatrix
 * TODO Though semantically NNO Scaling is diagonal in here they are represented
 * as Vectors
 * LinearSolvers should be aware of this and may be use DiagonalWrapper in
 * expressions
 */
class NTScalings {

 public:
  // Identity scaling for Initial Points
  // TODO Currently Identity is not calculated
  NTScalings() {}
  // Compute scalings for given point
  NTScalings(const Point& point)
      : NNO((point.s.cwiseQuotient(point.z)).cwiseSqrt()),
        NNOInverse(NNO.cwiseInverse()),
        NNOLambda((point.s.cwiseProduct(point.z)).cwiseSqrt()) {}
  // NNO -- NonNegativeOrthant, these are DiagonalMatrix but stored as Vectors
  const Eigen::VectorXd NNO;
  const Eigen::VectorXd NNOInverse;
  const Eigen::VectorXd NNOLambda;
};

std::ostream& operator<<(std::ostream& out, const NTScalings& scalings) {
  using namespace std;

  out << endl << "##################### NTScalings Start" << endl;
  out << "Value of NNO: " << endl << scalings.NNO << endl;
  out << "Value of NNOInverse: " << endl << scalings.NNOInverse << endl;
  out << "Value of NNOLambda: " << endl << scalings.NNOLambda << endl;
  out << "##################### NTScalings End" << endl;
}
}

#endif  // NT_SCALINGS_HPP