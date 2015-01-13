#ifndef ADMM_PROJECTIONS_HPP
#define ADMM_PROJECTIONS_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Core/Point.hpp>
#include <Problem.hpp>

namespace lp {

template <typename SubspaceProjection>
class ADMMProjections {
 public:
  ADMMProjections(const Problem& problem)
      : problem(problem), subspaceProjection(problem) {}

  Point doSubspaceProjection(const Point& currentPoint) {
    Point subspaceProjectedPoint(problem);

    subspaceProjectedPoint.x = currentPoint.x + currentPoint.r;
    subspaceProjectedPoint.z = currentPoint.z + currentPoint.s;
    subspaceProjectedPoint.tau = currentPoint.tau + currentPoint.kappa;

    subspaceProjection.doSubspaceProjection(subspaceProjectedPoint);

    relaxsubspaceProjectedPoint(currentPoint, subspaceProjectedPoint);

    return subspaceProjectedPoint;
  }

  Point doConeProjection(const Point& subspaceProjectedPoint,
                         const Point& currentPoint) {
    Point coneProjectedPoint(problem);

    coneProjectedPoint.x = subspaceProjectedPoint.x - currentPoint.r;
    coneProjectedPoint.z = subspaceProjectedPoint.z - currentPoint.s;
    coneProjectedPoint.tau = subspaceProjectedPoint.tau - currentPoint.kappa;

    // Projection
    coneProjectedPoint.z =
        coneProjectedPoint.z.unaryExpr([](const double & element)->double {
          if (element <= 0) {
            return 0;
          } else {
            return element;
          }
        });

    if (coneProjectedPoint.tau < 0) {
      coneProjectedPoint.tau = 0;
    }

    return coneProjectedPoint;
  }

  void doDualUpdates(const Point& subspaceProjectedPoint,
                     const Point& currentPoint, Point& coneProjectedPoint) {
    // TODO Why x varaible is not changed here
    coneProjectedPoint.r = currentPoint.r;
    coneProjectedPoint.s =
        currentPoint.s + coneProjectedPoint.z - subspaceProjectedPoint.z;
    coneProjectedPoint.kappa = currentPoint.kappa + coneProjectedPoint.tau -
                               subspaceProjectedPoint.tau;
  }

 private:
  const Problem& problem;
  SubspaceProjection subspaceProjection;

  // TODO Check why x is not relaxed here
  void relaxsubspaceProjectedPoint(const Point& currentPoint,
                                   Point& subspaceProjectedPoint) {
    subspaceProjectedPoint.z =
        subspaceProjectedPoint.z * problem.admmRelaxationParameter +
        (1 - problem.admmRelaxationParameter) * currentPoint.z;
    subspaceProjectedPoint.tau =
        subspaceProjectedPoint.tau * problem.admmRelaxationParameter +
        (1 - problem.admmRelaxationParameter) * currentPoint.tau;
  }
};
}
#endif  // ADMM_PROJECTIONS_HPP