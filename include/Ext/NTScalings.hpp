#ifndef NT_SCALINGS_HPP
#define NT_SCALINGS_HPP

#include <cmath>

#include <boost/mpl/int.hpp>
#include <boost/log/trivial.hpp>

#include <Eigen/Dense>

#include <Core/Point.hpp>
#include <Core/Residuals.hpp>
#include <Problem.hpp>

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
  NTScalings() : dg(1), dgInverse(1), lambdaG(1) {}
  // Compute scalings for given point
  NTScalings(const Point& point)
      : NNO((point.s.cwiseQuotient(point.z)).cwiseSqrt()),
        NNOInverse(NNO.cwiseInverse()),
        // Lambda belongs to Scalings, as This object knows how to best compute
        // lambda for different cones
        // semantically as well Scaling variable belongs to Scalings object
        // TODO I dont think we are using lambdaSquare anymore
        NNOLambdaSquare(point.s.cwiseProduct(point.z)),
        NNOLambda(NNOLambdaSquare.cwiseSqrt()),
        dg(std::sqrt(point.kappa / point.tau)),
        dgInverse(1 / dg),
        lambdaG(std::sqrt(point.kappa * point.tau)) {}
  // NNO -- NonNegativeOrthant, these are DiagonalMatrix but stored as Vectors
  const Eigen::VectorXd NNO;
  const Eigen::VectorXd NNOInverse;
  // Lambda and lambdaSquare are vectors
  const Eigen::VectorXd NNOLambdaSquare;
  const Eigen::VectorXd NNOLambda;
  // TODO I dont understand this part!!!
  const double dg;
  const double dgInverse;
  const double lambdaG;

  /**
   * Computes Rhs used in computation of affine direction
   * dxRhs = rx
   * dyRhs = ry
   * dzRhs = rz - W{T}(lambda o\ ds)
   *  Where ds = lambda o lambda
   * dtauRhs = rtau - dkappa/tau
   *  Where dkappa = kappa * tau
   *
   *
   * Rhs is used in following system
   *
   * [ 0  A'  G' ][x]   [dxRhs]
   * [-A  0   0  ][y] = [dyRhs]
   * [-G  0   W'W][z]   [dzRhs - W{T}(lambda o ds)]
   *
   * dtauRhs is used calculation of deltaTau
   *
   * Only o operator and diamond operator are done here (This ensures seperation
   *of concerns)
   * So only dzSubRhs = (lambda o\ ds) and dtauSubRhs = dk => kappa * tau are
   *caluclated here
   *
   * Returning point is actually not complete Rhs but subRhs with only z and tau
   *filled
   *
   * Using structure Point even if not all values are filled/used
   */
  Point getAffineSubRhs(const Problem& problem,
                        const Point& currentPoint) const {
    Point subRhs(problem.G.rows());

    subRhs.z = NNOLambda;
    subRhs.tau = currentPoint.kappa * currentPoint.tau;

    return subRhs;
  }

  /**
   * Returns z = lambda o\ ds
   * Returns tau = dk
   *
   */
  Point getCombinedSubRhs(const Problem& problem, const Point& currentPoint,
                          const Point& affinePoint, const double mu,
                          const double sigma) const {
    Point subRhs(problem.G.rows());

    // -sigma * mu * e + (W{-T} * sa) o (W * za) + lambda o lambda
    subRhs.z = NNOLambdaSquare + affinePoint.s.cwiseProduct(affinePoint.z) -
               (sigma * mu * Eigen::VectorXd::Ones(problem.h.rows()));
    // lambda o\ ds = diag(lambda){-1} * ds
    subRhs.z = subRhs.z.cwiseQuotient(NNOLambda);

    // -sigma * mu + ka*ta + k*t
    subRhs.tau = affinePoint.kappa * affinePoint.tau +
                 currentPoint.kappa * currentPoint.tau - (sigma * mu);

    return subRhs;
  }
};

std::ostream& operator<<(std::ostream& out, const NTScalings& scalings) {
  using namespace std;

  out << endl << "##################### NTScalings Start" << endl;
  out << "Value of NNO: " << endl << scalings.NNO << endl;
  out << "Value of NNOInverse: " << endl << scalings.NNOInverse << endl;
  out << "Value of NNOLambda: " << endl << scalings.NNOLambda << endl;
  out << "Value of NNOLambda Square: " << endl << scalings.NNOLambdaSquare
      << endl;
  out << "Value of dg: " << scalings.dg << endl;
  out << "Value of dgInverse: " << scalings.dgInverse << endl;
  out << "Value of lambdaG: " << scalings.lambdaG << endl;
  out << "##################### NTScalings End" << endl;

  return out;
}
}

#endif  // NT_SCALINGS_HPP