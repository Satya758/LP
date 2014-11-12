#ifndef SUITESPARSE_CHOLESKYLLT_HPP
#define SUITESPARSE_CHOLESKYLLT_HPP

#include <boost/log/trivial.hpp>
#include <boost/mpl/int.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include <Problem.hpp>

namespace lp {

// TODO Define new namespace
// TODO Handling different use cases is becoming complex, has to find better way
// for clean code
// TODO Review with peers to make code clean

template <const lp::SolveFor solveFor>
class EnumToType {

};
/**
 *Extending Eigen Cholmod wrapper as Eigen does not expose factor after
 *computing step, which is required in our use case with different inequality
 *and equality constaints
 */
class CholmodSupernodalLLTWithFactor
    : public Eigen::CholmodSupernodalLLT<
          Eigen::SparseMatrix<double, Eigen::Lower>> {

  typedef Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double, Eigen::Lower>>
      Base;

 public:
  CholmodSupernodalLLTWithFactor() : Base() {}

  ~CholmodSupernodalLLTWithFactor() { free(); }

  // TODO Do we need to return anythin in here as we are not using it
  // Free the variables that are created
  void compute(const Eigen::SparseMatrix<double>& matrix) {
    Base::compute(matrix);
    // factor is copied as it is lost for next compute
    // Only used when there are equality constraints in problem
    // See template secialization to achive this
    m_cholmodFactor_old =
        cholmod_copy_factor(Base::m_cholmodFactor, &Base::cholmod());
  }
  // Returns L (factored)
  // copy old factor using cholmod_copy_factor reason is
  // as compute is used again to factor new lhs, old factor
  // is lost
  Eigen::MappedSparseMatrix<double> getFactorAsEigen() {

    if (!m_cholmodFactor_old) {
      BOOST_LOG_TRIVIAL(error)
          << "Somthing wrong in creating/caching old factor";
    }

    Eigen::MappedSparseMatrix<double, Eigen::ColMajor, int> cholmodFactor(
        m_cholmodFactor_old->n, m_cholmodFactor_old->n,
        m_cholmodFactor_old->nzmax, static_cast<int*>(m_cholmodFactor_old->p),
        static_cast<int*>(m_cholmodFactor_old->i),
        static_cast<double*>(m_cholmodFactor_old->x));

    return cholmodFactor;
  }

  // Solves with custom factor
  // Copied from Eigen::CholmodSupport _solve method
  // TODO Had to find a better way!
  Eigen::VectorXd solveWithOldFactor(Eigen::VectorXd& rhs) {
    eigen_assert(m_factorizationIsOk &&
                 "The decomposition is not in a valid state for solving, you "
                 "must first call either compute() or symbolic()/numeric()");
    const Index size = m_cholmodFactor_old->n;
    EIGEN_UNUSED_VARIABLE(size);
    eigen_assert(size == rhs.rows());

    // note: cd stands for Cholmod Dense
    Eigen::VectorXd& rhs_ref(rhs.const_cast_derived());
    cholmod_dense rhs_cd = Eigen::viewAsCholmod(rhs_ref);
    cholmod_dense* x_cd = cholmod_solve(CHOLMOD_A, m_cholmodFactor_old, &rhs_cd,
                                        &Base::cholmod());
    if (!x_cd) {
      this->m_info = Eigen::NumericalIssue;
    }
    // TODO optimize this copy by swapping when possible (be careful with
    // alignment, etc.)
    Eigen::VectorXd solution = Eigen::VectorXd::Map(
        reinterpret_cast<double*>(x_cd->x), rhs.rows(), rhs.cols());
    cholmod_free_dense(&x_cd, &Base::cholmod());

    return solution;
  }

  // its a destructor, as normal destructor is not called
  void free() {
    if (m_cholmodFactor_old) {
      // TODO Check if we have to send address of m_cholmodFactor_old
      // as in analyze method in Base class
      cholmod_free_factor(&m_cholmodFactor_old, &Base::cholmod());
      m_cholmodFactor_old = NULL;
    }
  }

 private:
  cholmod_factor* m_cholmodFactor_old = NULL;
};

/**
 *SuiteSparse cholmod
 *TODO Move below comment to IPM algorithm class
 *Q. Why call factor, why not do it in constructor of Linear solver provided
 *A. Hmmm... Both approaches have their own advantages and disadvantages, 1.
 *factor in constructor is much clean and straight but user (creator of the
 *linear solver) will not have
 *any context of previous iteration as object is destroyed at the end of each
 *iteration in IPM. But advantage is user doesn't have to delete any objects at
 *each iteration that he
 *has created to find solution, he can rely on destructor to do proper cleaning
 *in face of any disater 2. Seperate factor method gives user flexibiltiy to
 *store context information and use it later to improve performace like reuse
 *analyze, update/downdate. But user have to maintain objects created carefully.
 *As performace is utmost importance second option is chosen
 *
 */
class SuiteSparseCholeskyLLT {

 public:
  SuiteSparseCholeskyLLT(const lp::Problem& problem) : problem(problem) {
    BOOST_LOG_TRIVIAL(info) << "Constructor of Linear Solver";
  }

  // Factor matrix omega = G' * W{-1} * W{-T} * G
  // If scalingMatrix W is Identity factor G' * G
  // If A (equality constaints) is not empty then factor A * omega{-1} * A' to
  // find y
  // First find y, then x and finally z
  // Diagonal Scaling matrix is only good for Linear programming
  void factor(
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrix,
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrixInverse,
      constexpr lp::SolveFor solveFor) {
    BOOST_LOG_TRIVIAL(info) << "Step directions Factorization";

    // TODO Duplication of code how do we avoid it!?
    solver.free();
    // First compute G' * W{-1} and store it
    omega = getOmega(scalingMatrixInverse, EnumToType<solveFor>());
     //omegaTilde =
       //  getOmegaTilde(scalingMatrixInverse, boost::mpl::int_<lp::SolveFor::Initial>());

    if (problem.A.size() != 0) {
      // Factor G' * W{-1} * W{-T} * G
      solver.compute(omega * omega.transpose());
      // Returns L
      inequalityFactor = solver.getFactorAsEigen();

      // L{-1} * A'
      // solve does not exist for sparse matrix, using solveInPlace, but
      // solveInPlace does not allow Eigen::OnTheLeft flag!
      // TODO There is some problem when problem.A.transpose() is used directly
      // in solveInPlace method instead of temporary
      // inequalityLInverseATranspose
      // TODO Invoke lazy evaluation if creating temporary variable is pain
      // Check http://eigen.tuxfamily.org/dox/TopicLazyEvaluation.html
      inequalityLInverseATranspose = problem.A.transpose();
      // TODO FIXME Instead of using traingular solve use cholmod_solve check
      // documentation
      // https://www.cise.ufl.edu/research/sparse/cholmod/CHOLMOD/Doc/UserGuide.pdf
      inequalityFactor.triangularView<Eigen::Lower>().solveInPlace(
          inequalityLInverseATranspose);
      // Compute A * L{-T} * L{-1} * A'
      solver.compute(inequalityLInverseATranspose.transpose() *
                     inequalityLInverseATranspose);
    } else {
      solver.compute(omega * omega.transpose());
    }

    if (solver.info() != Eigen::Success) {
      // TODO Raise exception maybe Matrix is singular
      BOOST_LOG_TRIVIAL(error) << "Factorization failed";
    }
  }

  // rhsX, rhsY, rhsZ are not primal/dual variables but symbolic names of right
  // hand side of three equations
  // When A(equality constaints) is empty rhsY is not required/considered
  // When problem is in standard for without inequality constraints rhs has only
  // two vectors
  // First calculate y from rhs
  // (A * L{-T} * L{-1} * A') * y = A * L{-T} * L{-1} * (rhsX + G' * W{-1} *
  // W{-T} * rhsZ + A' * rhsY) - rhsY
  // then calculate x from rhs
  // (L * L') * x = rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
  // finally calculate z
  // z = W{-1} * W{-T} * (G*x - rhsZ)
  lp::NewtonDirection solve(
      const Eigen::VectorXd& rhsX, const Eigen::VectorXd& rhsY,
      const Eigen::VectorXd& rhsZ,
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrixInverse,
      const lp::SolveFor solveFor) const {

    lp::NewtonDirection direction;

    if (problem.A.size() != 0) {
      {
        // L{-1} * (rhsX + G' * W{-1} * W{-T} * rhsZ + A' * rhsY)
        Eigen::VectorXd subRhs =
            rhsX + omegaTilde * rhsZ + problem.A.transpose() * rhsY;
        inequalityFactor.triangularView<Eigen::Lower>().solveInPlace(subRhs);

        Eigen::VectorXd rhs =
            inequalityLInverseATranspose.transpose() * subRhs - rhsY;
        direction.y = solver.solve(rhs);
      }
      {
        // rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
        Eigen::VectorXd rhs = rhsX + omegaTilde * rhsZ +
                              problem.A.transpose() * (rhsY - direction.y);
        // Removing constantness just for this call
        direction.x = const_cast<SuiteSparseCholeskyLLT*>(this)
                          ->solver.solveWithOldFactor(rhs);
      }
      { direction.z = problem.G * direction.x - rhsZ; }
    } else {
      {
        // rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
        Eigen::VectorXd rhs = rhsX + problem.G.transpose() * rhsZ;

        direction.x = solver.solve(rhs);
      }
      { direction.z = problem.G * direction.x - rhsZ; }
    }

    // TODO Is there way to enable this if condition during compile time
    // enable_if only works for methods not conditional statements
    if (solveFor == lp::SolveFor::StepDirection) {
      // TODO check why below statement fails
      // direction.z = scalingMatrixInverse * scalingMatrixInverse *
      // direction.z;
    }

    return direction;
  }

 private:
  // G'
  Eigen::SparseMatrix<double> getOmega(
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrixInverse,
     EnumToType<lp::SolveFor::Initial>) const {
    // TODO For initial, G Transpose is done twice unnecesarly
    // but I think its OK, as it is only one time, for the convinience of code
    return problem.G.transpose();
  }

  // G' * W{-1}
  Eigen::SparseMatrix<double> getOmega(
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrixInverse,
      EnumToType<lp::SolveFor::StepDirection>) const {
    return problem.G.transpose() * scalingMatrixInverse;
  }

  // G' * W{-1} * W{-T}
  // For initial its just G'
  Eigen::SparseMatrix<double> getOmegaTilde(
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrixInverse,
      boost::mpl::int_<lp::SolveFor::Initial>) const {
    return problem.G.transpose();
  }

  // G' * W{-1} * W{-T}
  // Transpose of diagonal matrix alters nothing
  Eigen::SparseMatrix<double> getOmegaTilde(
      const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& scalingMatrixInverse,
      boost::mpl::int_<lp::SolveFor::StepDirection>) const {
    return omega * scalingMatrixInverse;
  }

  const lp::Problem& problem;
  // Custom Solver
  CholmodSupernodalLLTWithFactor solver;
  // All below variables are used to store intermediate results to avoid
  // repeated calculations or copies (Product computaions are coslty)

  // G' * W{-1} * W{-T} Should be used in multiple solves to solve
  Eigen::SparseMatrix<double> omegaTilde;  // rhsZCoefficient;
  // L
  Eigen::SparseMatrix<double> inequalityFactor;
  // L{-1} * A'
  Eigen::SparseMatrix<double> inequalityLInverseATranspose;
  // G' * W{-1} Lets call it Omega
  Eigen::SparseMatrix<double> omega;
};
}

#endif  // SUITESPARSE_CHOLESKYLLT_HPP