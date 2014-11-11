#ifndef SUITESPARSE_CHOLESKYLLT_HPP
#define SUITESPARSE_CHOLESKYLLT_HPP

#include <boost/log/trivial.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

#include <Problem.hpp>

// TODO Define new namespace

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

    m_cholmodFactor_old =
        cholmod_copy_factor(Base::m_cholmodFactor, &Base::cholmod());
  }
  // Returns L (factored)
  // copy old factor using cholmod_copy_factor reason is
  // as compute is used again to factor new lhs, old factor
  // is lost
  Eigen::MappedSparseMatrix<double> getCholmodFactorAsEigen() {

    if (!m_cholmodFactor_old) {
      BOOST_LOG_TRIVIAL(error) << "Somthing wrong in creating factor";
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
  void factor(Eigen::DiagonalMatrix<double, Eigen::Dynamic> scalingMatrix) {
    BOOST_LOG_TRIVIAL(info) << "Factorization";
  }

  // When calculating Initial points use this method as scaling matrix is
  // Identity matrix
  // Q. Why another method instead of above method?
  // A. Performance...We dont have to check for identity condition and interface
  // is clean (ofcourse you have to understand algorithm to see clean
  // interface :-) ).
  // Use this when scaling matrix is identity matrix
  void factor() {

    // First delete any objects from previous iteration
    // this should be done in another factor as well
    // TODO How do we avoid this duplication?
    solver.free();

    forInitial = true;
    // Factor G' * G
    solver.compute(problem.G.transpose() * problem.G);
    // Check for existence of equality constaints
    if (problem.A.size() != 0) {
      // Returns L
      inequalityFactor = solver.getCholmodFactorAsEigen();

      // L{-1} * A'
      // solve does not exist for sparse matrix, using solveInPlace, but
      // solveInPlace does not allow Eigen::OnTheLeft flag!
      inequalityLInverseATranspose = problem.A.transpose();

      inequalityFactor.triangularView<Eigen::Lower>().solveInPlace(
          inequalityLInverseATranspose);
      // Compute A * L{-T} * L{-1} * A'
      solver.compute(inequalityLInverseATranspose.transpose() *
                     inequalityLInverseATranspose);
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
  void solve(const Eigen::VectorXd& rhsX, const Eigen::VectorXd& rhsY,
             const Eigen::VectorXd& rhsZ) {
    // TODO Refactor after testing
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;
    // TODO Initial would not work, i.e. when do we turn off flag after solve...
    // No as solve for initial can be called for multiple sides
    if (forInitial) {
      if (problem.A.size() != 0) {
        {
          // L{-1} * (rhsX + G' * W{-1} * W{-T} * rhsZ + A' * rhsY)
          Eigen::VectorXd subRhs = rhsX + problem.G.transpose() * rhsZ +
                                   problem.A.transpose() * rhsY;
          inequalityFactor.triangularView<Eigen::Lower>().solveInPlace(subRhs);

          Eigen::VectorXd rhs =
              inequalityLInverseATranspose.transpose() * subRhs - rhsY;
          y = solver.solve(rhs);
        }
        {
          // rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
          Eigen::VectorXd rhs = rhsX + problem.G.transpose() * rhsZ +
                                problem.A.transpose() * (rhsY - y);

          x = solver.solveWithOldFactor(rhs);
        }
        { z = problem.G * x - rhsZ; }
      } else {
        {
          // rhsX + G' * W{-1} * W{-T} * rhsZ + A' * (rhsY - y)
          Eigen::VectorXd rhs = rhsX + problem.G.transpose() * rhsZ;

          x = solver.solve(rhs);
        }
        { z = problem.G * x - rhsZ; }
      }
    }

    BOOST_LOG_TRIVIAL(info) << x;
    BOOST_LOG_TRIVIAL(info) << z;
  }

 private:
  lp::Problem problem;
  // Custom Solver
  CholmodSupernodalLLTWithFactor solver;
  // All below variables are used to store intermediate results to avoid
  // repeated calculations or copies
  // Indication for solve method to use rhsZCoefficient or not
  // TODO I dont think forInitial will work
  bool forInitial;
  // W{-1}
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> scalingMatrixInverse;
  // G' * W{-1} * W{-T}
  Eigen::SparseMatrix<double> rhsZCoefficient;
  // L
  Eigen::SparseMatrix<double> inequalityFactor;
  // L{-1} * A' // TODO May be this object is better maintained by
  // CholmodSupernodalLLTWithFactor object
  Eigen::SparseMatrix<double> inequalityLInverseATranspose;
};

#endif  // SUITESPARSE_CHOLESKYLLT_HPP