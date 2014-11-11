
#include <boost/log/trivial.hpp>

#include <Problem.hpp>
#include <Solver.hpp>

#include <SuiteSparseCholeskyLLT.hpp>

// Had to find with equality constraints as well, this is where
// lot of mumbo jumbo of code is present
lp::Problem getInequalityTest() {
  lp::Problem problem(4, 0, 2);

  problem.c(0) = 2;
  problem.c(1) = 1;

  problem.h(0) = 1;
  problem.h(1) = -2;
  problem.h(2) = 0;
  problem.h(3) = 4;

  problem.G.insert(0, 0) = -1;
  problem.G.insert(0, 1) = 1;

  problem.G.insert(1, 0) = -1;
  problem.G.insert(1, 1) = -1;

  problem.G.insert(2, 1) = -1;

  problem.G.insert(3, 0) = 1;
  problem.G.insert(3, 1) = -2;

  problem.G = problem.G * 1;
  BOOST_LOG_TRIVIAL(info) << problem.c;

  return problem;
}

int main(int argc, char **argv) {

  BOOST_LOG_TRIVIAL(info) << "Hello World!";

  // lp::Problem problem(5, 5, 5);
  // TODO Internal is used but should not be
  // Refactor after testing
  lp::internal::Solver<SuiteSparseCholeskyLLT> solver(getInequalityTest());
  return 0;
}
