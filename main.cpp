
// #include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
// #include <boost/log/expressions.hpp>

// #include <Problem.hpp>
// #include <Solution.hpp>
// #include <Core/Solver.hpp>
//
// #include <Ext/NTScalings.hpp>
// #include <Ext/SuiteSparseCholeskyLLT.hpp>

#include <Presolve/LPFormatParser.hpp>

// // Had to find with equality constraints as well, this is where
// // lot of mumbo jumbo of code is present
// lp::Problem getInequalityTest() {
//   lp::Problem problem(4, 0, 2);
//
//   problem.c(0) = 2;
//   problem.c(1) = 1;
//
//   problem.h(0) = 1;
//   problem.h(1) = -2;
//   problem.h(2) = 0;
//   problem.h(3) = 4;
//
//   problem.G.insert(0, 0) = -1;
//   problem.G.insert(0, 1) = 1;
//
//   problem.G.insert(1, 0) = -1;
//   problem.G.insert(1, 1) = -1;
//
//   problem.G.insert(2, 1) = -1;
//
//   problem.G.insert(3, 0) = 1;
//   problem.G.insert(3, 1) = -2;
//
//   return problem;
// }
//
// // Taken from
// // http://www.vitutor.com/alg/linear_programming/problems_solutions.html
// lp::Problem getVitutorTest1() {
//   lp::Problem problem(4, 0, 2);
//
//   problem.c(0) = 30;
//   problem.c(1) = 40;
//
//   problem.h(0) = -3000;
//   problem.h(1) = -4000;
//   problem.h(2) = 0;
//   problem.h(3) = 0;
//
//   problem.G.insert(0, 0) = -20;
//   problem.G.insert(0, 1) = -30;
//
//   problem.G.insert(1, 0) = -40;
//   problem.G.insert(1, 1) = -30;
//
//   problem.G.insert(2, 0) = -1;
//
//   problem.G.insert(3, 1) = -1;
//
//   return problem;
// }
//
// lp::Problem getVitutorTest2() {
//   lp::Problem problem(4, 0, 2);
//
//   problem.c(0) = 600;
//   problem.c(1) = 800;
//
//   problem.h(0) = -400;
//   problem.h(1) = 9;
//   problem.h(2) = 0;
//   problem.h(3) = 0;
//
//   problem.G.insert(0, 0) = -40;
//   problem.G.insert(0, 1) = -50;
//
//   problem.G.insert(1, 0) = 1;
//   problem.G.insert(1, 1) = 1;
//
//   problem.G.insert(2, 0) = -1;
//
//   problem.G.insert(3, 1) = -1;
//
//   return problem;
// }
//
// lp::Problem getVitutorTest3() {
//   lp::Problem problem(4, 0, 2);
//
//   problem.c(0) = 30;
//   problem.c(1) = 50;
//
//   problem.h(0) = 200;
//   problem.h(1) = 100;
//   problem.h(2) = -20;
//   problem.h(3) = -10;
//
//   problem.G.insert(0, 0) = 1;
//   problem.G.insert(0, 1) = 3;
//
//   problem.G.insert(1, 0) = 1;
//   problem.G.insert(1, 1) = 1;
//
//   problem.G.insert(2, 0) = -1;
//
//   problem.G.insert(3, 1) = -1;
//
//   return problem;
// }
//
// void testSolver() {
//   boost::log::core::get()->set_filter(boost::log::trivial::severity >=
//                                       boost::log::trivial::info);
//
//   BOOST_LOG_TRIVIAL(info) << "Started...";
//
//   {
//     lp::Problem problem = getVitutorTest1();
//     lp::Solver<lp::SuiteSparseCholeskyLLT<lp::NTScalings>, lp::NTScalings>
//         solver(problem);
//
//     lp::Solution solution = solver.solve();
//
//     BOOST_LOG_TRIVIAL(info) << solution;
//   }
//
//   {
//     lp::Problem problem = getVitutorTest2();
//     lp::Solver<lp::SuiteSparseCholeskyLLT<lp::NTScalings>, lp::NTScalings>
//         solver(problem);
//
//     lp::Solution solution = solver.solve();
//
//     BOOST_LOG_TRIVIAL(info) << solution;
//   }
//
//   {
//     lp::Problem problem = getVitutorTest3();
//     lp::Solver<lp::SuiteSparseCholeskyLLT<lp::NTScalings>, lp::NTScalings>
//         solver(problem);
//
//     lp::Solution solution = solver.solve();
//
//     BOOST_LOG_TRIVIAL(info) << solution;
//   }
//
//   {
//     lp::Problem problem = getInequalityTest();
//     lp::Solver<lp::SuiteSparseCholeskyLLT<lp::NTScalings>, lp::NTScalings>
//         solver(problem);
//
//     lp::Solution solution = solver.solve();
//
//     BOOST_LOG_TRIVIAL(info) << solution;
//   }
// }

int main(int argc, char **argv) {

  lp::parser::LPFormatParser parser;

  BOOST_LOG_TRIVIAL(info) << "Started";

  //   parser.parse("/home/satya/Desktop/sample.lp");
  parser.parse("/home/satya/Desktop/test.lp");

  BOOST_LOG_TRIVIAL(info) << "Ended";
  return 0;
}
