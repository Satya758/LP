#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <algorithm>

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>

#include <Core/Point.hpp>
#include <Core/Residuals.hpp>
#include <Core/Timer.hpp>
#include <Problem.hpp>
#include <Solution.hpp>
// FIXME If there is chance of using mathematical symbols in text editor use it
// Copy symbols from http://math.typeit.org/
namespace lp {

// Thought of using in Scalings object to determine computation of Rhs for
// Affine/Combined
// Unable to use due to signature difference between Affine and Combined
// direction
enum class Direction {
  Affine,
  Combined
};

/**
 *Solves Convex problems
 *Linear solver is provided as template argument
 */
template <typename LinearSolver, typename Scalings>
class IPMSolver {
 public:
  IPMSolver(const Problem& problem) : problem(problem), linearSolver(problem) {}

  /**
   *
   */
  Solution solve() {

    Timer& timer = Timer::getIPMInstance();
    // FIXME Remove these lines if not required
    //     const double residualX0{std::max(1.0, problem.c.norm())};
    //     const double residualZ0{std::max(1.0, problem.h.norm())};
    const double residualX0{1.0 + problem.c.norm()};
    const double residualZ0{1.0 + problem.h.norm()};

    timer.start("Initial Point");
    Point currentPoint = getInitialPoint();
    timer.end("Initial Point");

    Residual tolerantResidual(problem);

    // TODO Either start from 1 or use just lessthan instead of lessthan or eual
    // to operator, or maybe I can use these indications to find start/end of
    // loop
    for (int i = 0; i <= problem.ipmMaxIterations; ++i) {

      timer.start("Residual Computation");
      Residual residual(problem, currentPoint, i, residualX0, residualZ0);
      timer.end("Residual Computation");

      BOOST_LOG_TRIVIAL(info) << residual;

      timer.start("SolverState Computation");
      SolverState solverState = getSolverState(residual, tolerantResidual);
      timer.end("SolverState Computation");

      if (solverState == SolverState::Feasible) {
        return Solution(residual, currentPoint, solverState);
      } else if (solverState == SolverState::PrimalInfeasible ||
                 solverState == SolverState::DualInfeasible) {
        //         return Solution(residual, currentPoint, solverState);
        // FIXME return after failure
        BOOST_LOG_TRIVIAL(info) << "Problem is Infeasible";
      } else if (i == problem.ipmMaxIterations) {
        return Solution(residual, currentPoint, SolverState::MaximumIterations);
      }

      timer.start("Scalings Compute");
      // Compute scalings for given point
      Scalings scalings(currentPoint);
      timer.end("Scalings Compute");

      timer.start("Affine Direction");
      NewtonDirection subSolution = getSubSolution(scalings);

      Point affineDirection = this->template getDirection<Direction::Affine>(
          scalings, residual, subSolution, currentPoint);
      timer.end("Affine Direction");

      BOOST_LOG_TRIVIAL(trace) << "Iteration: " << i;
      BOOST_LOG_TRIVIAL(trace) << "Scalings " << scalings;
      BOOST_LOG_TRIVIAL(trace) << "Affine Direction" << affineDirection;

      timer.start("StepSize Compute");
      double alpha = this->template computeAlpha<Direction::Affine>(
          affineDirection, scalings);
      double sigma = std::pow(1 - alpha, exp);

      // mu = lambda' * lambda / (1 + m)
      // TODO m will be different when other cones are introduced and what
      // should be when A is used?
      double mu = (scalings.NNOLambda.squaredNorm() +
                   currentPoint.kappa * currentPoint.tau) /
                  (1 + problem.h.rows());

      timer.end("StepSize Compute");
      BOOST_LOG_TRIVIAL(trace) << "Mu: " << mu;
      BOOST_LOG_TRIVIAL(trace) << "Sigma: " << sigma;
      BOOST_LOG_TRIVIAL(trace) << "1 - Sigma: " << 1 - sigma;

      timer.start("Combined Direction");
      Point combinedDirection =
          this->template getDirection<Direction::Combined>(
              scalings, residual, subSolution, currentPoint, affineDirection,
              mu, sigma);

      BOOST_LOG_TRIVIAL(trace) << "Combined Direction" << combinedDirection;
      // Step size updated
      alpha = this->template computeAlpha<Direction::Combined>(
          combinedDirection, scalings);

      // Unscale s and z for step updates
      combinedDirection.s = scalings.NNO.cwiseProduct(combinedDirection.s);
      combinedDirection.z =
          scalings.NNOInverse.cwiseProduct(combinedDirection.z);
      // TODO Find out how unscaling works
      combinedDirection.tau = combinedDirection.tau * scalings.dgInverse;
      combinedDirection.kappa = combinedDirection.kappa * scalings.dg;
      timer.end("Combined Direction");

      timer.start("Direction Update");
      currentPoint = updatePoint(currentPoint, combinedDirection, alpha);

      BOOST_LOG_TRIVIAL(trace) << "Updated current point" << currentPoint;
      timer.end("Direction Update");
    }
  }

 private:
  const Problem& problem;
  LinearSolver linearSolver;

  constexpr static double exp = 3;
  constexpr static double step = 0.99;

  // As name implies its dummy object as Initial point does not need any
  // scaling point
  // TODO Check if creating and passing dummy object can be problamatic for
  // perforamnce
  // TODO Check what is the cost of creating dummy object, if it is high had
  // to find another way
  // TODO May be not much of a problem as it is used only once
  // TODO As Interface dicatates I have to send Identity Matrix, as I know the
  // implementation inside I have taken shortcut but its wrong
  Point getInitialPoint() {

    Point point(problem);
    Scalings scalings;

    linearSolver.template factor<lp::SolveFor::Initial>(scalings);

    {
      // Get Primal initial points
      Eigen::VectorXd rhsX(problem.c.rows());
      rhsX.setZero();

      NewtonDirection direction =
          linearSolver.template solve<lp::SolveFor::Initial>(rhsX, problem.h,
                                                             scalings);
      // TODO FIXME Can we invoke move constructor instead of copying data?
      point.x = direction.x;
      point.s = -1 * direction.z;

      adjustInitialPoint(point.s);
    }
    {
      // Get Dual initial points
      Eigen::VectorXd rhsX(-1 * problem.c);
      Eigen::VectorXd rhsZ(problem.G.rows());
      rhsZ.setZero();

      NewtonDirection direction =
          linearSolver.template solve<lp::SolveFor::Initial>(rhsX, rhsZ,
                                                             scalings);
      // FIXME Move constructor?
      point.z = direction.z;

      adjustInitialPoint(point.z);
    }
    {
      point.kappa = 1;
      point.tau = 1;
    }
    return point;
  }

  /**
   * vector = vector | alpha < 0
   * vector = vector + (1 + alpha) | otherwise
   *
   */
  void adjustInitialPoint(Eigen::VectorXd& vector) {
    double alpha = getMaxStep(vector);

    if (alpha >= -1e-8 * std::max(vector.norm(), 1.0)) {
      vector = (vector.array() + 1 + alpha).matrix();
    }
  }

  /**
   * Finds solution to sub KKT system
   *
   * [0   G'][x]  [-c]
   * [G -W'W][z]  [h]
   *
   * TODO result is multiplied with dgi...why
   */
  NewtonDirection getSubSolution(const Scalings& scalings) {
    // Only one factorization in loop is done here
    linearSolver.template factor<lp::SolveFor::StepDirection>(scalings);

    Eigen::VectorXd rhsX(-1 * problem.c);

    NewtonDirection direction =
        linearSolver.template solve<lp::SolveFor::StepDirection>(
            rhsX, problem.h, scalings);
    // TODO Why
    return scalings.dgInverse * direction;
  }

  /**
   *
   * o\ is inverse operation for o (diamond is used in cvxopt document)
   * return value of z is scaled, actual value is W * z
   * FIXME Why Scaled z
   *
   * [0    G'][x]    [dx]
   * [G  -W'W][z]    [W{T}(lambda o\ ds) - dz]
   *
   * Final Rhs is computed here with help from Scalings object
   *
   */
  template <lp::Direction direction>
  Point getDirection(const Scalings scalings, const Residual& residual,
                     const NewtonDirection& subSolution,
                     const Point& currentPoint,
                     // Used only for Correction direction
                     const Point& affinePoint = Point(), const double mu = 0,
                     const double sigma = 0) {

    const double residualMul = 1.0 - sigma;

    // Returns sub values which have o and o\ as operators
    // in z we find lambda o\ ds (ds depends on Direction)
    // in tau we find dk (dk depends on Direction)
    Point subRhs;
    if (direction == Direction::Affine) {
      subRhs = scalings.getAffineSubRhs(problem, currentPoint);
    } else if (direction == Direction::Combined) {
      subRhs = scalings.getCombinedSubRhs(problem, currentPoint, affinePoint,
                                          mu, sigma);
    }

    // rhsZ = rz - W{T} * (lambda o\ ds)
    // Notice Scaling matrix is not transposed here, in LP scenario we don't
    // have to
    // transpose Scaling matrix as it is diagonal, TODO but it should be done
    // differently for other cones
    Eigen::VectorXd rhsZ = (scalings.NNO.asDiagonal() * subRhs.z) -
                           residual.residualZ * residualMul;

    BOOST_LOG_TRIVIAL(trace) << "rhsX" << residual.residualX * residualMul;
    BOOST_LOG_TRIVIAL(trace) << "rhsZ" << rhsZ;

    NewtonDirection subSolution2 =
        linearSolver.template solve<lp::SolveFor::StepDirection>(
            residual.residualX * residualMul, rhsZ, scalings);

    BOOST_LOG_TRIVIAL(trace) << "Solution from LS: " << subSolution2;

    Point point(problem);

    // deltaTau = (dTau - dk/tau + c'x2 + b'y2 + W{-T}h'z2)/ kappa/tau +
    // ||Wz||^{2}
    // x2,y2,z2 are newtonDirection calculated above
    // W{-T}h'z2 instead of h'z2 is because of the way solution is calculated
    // x1,y1,z1 are from subSolution
    // TODO Check why do we have to multiply with dgInverse
    point.tau =
        scalings.dgInverse *
        (residual.residualTau * residualMul - subRhs.tau / currentPoint.tau +
         problem.c.dot(subSolution2.x) +
         (scalings.NNOInverse.asDiagonal() * problem.h).dot(subSolution2.z)) /
        (1 + subSolution.z.squaredNorm());

    // x = x2 + tau * x1
    point.x = subSolution2.x + point.tau * subSolution.x;
    // z = z2 + tau * z1
    point.z = subSolution2.z + point.tau * subSolution.z;

    // XXX
    // s = -W{T} * ( lambda o\ ds + W*z)
    // => W{-T}*s = -( lambda o\ ds + W*z)
    // W{-T}*s is used in combined search direction
    // TODO Check this... Like scaled Z here we have scaled S
    // point.s = -1 * scalings.NNO.asDiagonal() * (subRhs.z + point.z);
    point.s = -subRhs.z - point.z;
    // kappa = -(dk + currentKappa * tau)/currentTau;
    //     point.kappa =
    //         -(subRhs.tau + currentPoint.kappa * point.tau) /
    // currentPoint.tau;
    // TODO Check lambdaG
    point.kappa = -(subRhs.tau / scalings.lambdaG) - point.tau;

    return point;
  }

  template <Direction direction>
  double computeAlpha(const Point& searchDirection,
                      const Scalings& scalings) const {
    // Expand the s and z with lambda
    // TODO Check what is theoritcal advantage
    Eigen::VectorXd s = searchDirection.s.cwiseQuotient(scalings.NNOLambda);
    Eigen::VectorXd z = searchDirection.z.cwiseQuotient(scalings.NNOLambda);

    double ts = getMaxStep(s);
    double tz = getMaxStep(z);

    double maxCoeff = std::max({0.0,
                                ts,
                                tz,
                                -searchDirection.kappa / scalings.lambdaG,
                                -searchDirection.tau / scalings.lambdaG});

    if (maxCoeff == 0) {
      return 1.0;
    } else {
      if (direction == Direction::Affine) {
        return std::min(1.0, 1.0 / maxCoeff);
      } else if (direction == Direction::Combined) {
        return std::min(1.0, step / maxCoeff);
      }
    }
  }
};

}  // lp

#endif  // SOLVER_HPP
