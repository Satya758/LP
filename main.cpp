
#include <string>
#include <iostream>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>

#include <Problem.hpp>
#include <Solution.hpp>
#include <Core/IPMSolver.hpp>
#include <Core/ADMMSolver.hpp>

#include <Ext/NTScalings.hpp>
#include <Ext/IPMCholeskyLLT.hpp>

#include <Presolve/LPFormatParser.hpp>

// SCS
#include <Ext/ADMMCholeskyLDLT.hpp>

extern "C" {
#include "scs.h"
#include "../linsys/amatrix.h"
}

#include <random>
#include <iostream>

#define EXTRAVERBOSE

void callSCS(lp::Problem& problem) {

  Cone* k;
  Data* d;
  Work* w;
  Sol* sol;

  Info info = {0};

  k = (Cone*)scs_calloc(1, sizeof(Cone));
  d = (Data*)scs_calloc(1, sizeof(Data));
  sol = (Sol*)scs_calloc(1, sizeof(Sol));

  AMatrix* Araw = (AMatrix*)malloc(sizeof(AMatrix));

  problem.G.makeCompressed();

  d->m = problem.G.rows();
  d->n = problem.G.cols();

  d->A = Araw;

  Araw->x = problem.G.valuePtr();
  Araw->i = problem.G.innerIndexPtr();
  Araw->p = problem.G.outerIndexPtr();

  d->b = problem.h.data();
  d->c = problem.c.data();

  k->l = problem.G.rows();

  d->max_iters = 1000000;
  d->eps = 1e-1;
  d->alpha = 1.5;
  d->rho_x = 1e-3;
  d->scale = 1;
  d->verbose = true;
  d->normalize = 1;

  scs(d, k, sol, &info);

  using namespace std;

  cout << info.status << endl;

  cout << sol->x[0] << "----" << sol->x[1] << endl;
  cout << sol->y[0] << "----" << sol->y[1] << endl;
  cout << sol->s[0] << "----" << sol->s[1] << endl;

  cout << "Primal R: " << info.resPri << endl;
  cout << "Dual R: " << info.resDual << endl;
  cout << "Duality Gap: " << info.relGap << endl;
  cout << "Number of iterations: " << info.iter << endl;
  cout << "Time taken: " << info.solveTime << endl;
}

class CommandOptions {
 public:
  std::string fileName;
  std::string choleskySolver = "CholmodLLT";
  int logOptions;
};

CommandOptions getOptions(int argc, char** argv) {
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()("help", "LP Format file name to solve and log options")(
      "file,f", po::value<std::string>(), "Enter file path name")(
      "CholeskySolver,c", po::value<std::string>(),
      "Enter solver name, default is Cholmod LLT");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
  }

  CommandOptions options;

  if (vm.count("file")) {
    options.fileName = vm["file"].as<std::string>();
  } else {
    std::cout << "File name is not given " << std::endl;
  }

  if (vm.count("CholeskySolver")) {
    options.choleskySolver = vm["CholeskySolver"].as<std::string>();
  }

  return options;
}

int main(int argc, char** argv) {
  // testSolver();
  boost::log::core::get()->set_filter(boost::log::trivial::severity >=
                                      boost::log::trivial::info);

  CommandOptions options = getOptions(argc, argv);

  lp::parser::LPFormatParser parser;

  BOOST_LOG_TRIVIAL(info) << "Started";

  if (options.fileName.empty()) {
    std::cout << "File name is not given " << std::endl;
    return 1;
  }

  BOOST_LOG_TRIVIAL(info) << "Parser and Presolve Started";
  lp::Problem problem = parser.parse(options.fileName);
  BOOST_LOG_TRIVIAL(info) << "Parser and Presolve Ended";

  //     typedef lp::IPMCholeskyLLT<
  //         lp::NTScalings,
          Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double,
    Eigen::Lower>>
          CholmodSolver;
  //     typedef lp::IPMCholeskyLLT<lp::NTScalings,
  //                                Eigen::PastixLDLT<Eigen::SparseMatrix<double>,
  //                                                  Eigen::Lower>>
  // PastixSolver;
  //
  //     BOOST_LOG_TRIVIAL(info) << "Started to solve";
  //
  //     if (options.choleskySolver == "CholmodLLT") {
  //       lp::IPMSolver<CholmodSolver, lp::NTScalings> solver(problem);
  //       lp::Solution solution = solver.solve();
  //       BOOST_LOG_TRIVIAL(info) << solution;
  //     } else {
  //       lp::IPMSolver<PastixSolver, lp::NTScalings> solver(problem);
  //       lp::Solution solution = solver.solve();
  //       BOOST_LOG_TRIVIAL(info) << solution;
  //     }
  //
  //     BOOST_LOG_TRIVIAL(info) << "Ended";

  //   typedef lp::ADMMCholeskyLDLT<Eigen::PastixLDLT<
  //       Eigen::SparseMatrix<double>, Eigen::Lower>> ADMMPastixSolver;
  //
  //   lp::ADMMSolver<ADMMPastixSolver> solver(problem);
  //   lp::Solution solution = solver.solve();
  //   BOOST_LOG_TRIVIAL(info) << solution;

  callSCS(problem);

  return 0;
}
