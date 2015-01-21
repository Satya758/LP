
#include <string>
#include <iostream>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>

#include <Core/Timer.hpp>
#include <Problem.hpp>
#include <Solution.hpp>
#include <Core/IPMSolver.hpp>

#include <Core/ADMMSolver.hpp>

#include <Ext/NTScalings.hpp>
#include <Ext/IPMCholeskyLLT.hpp>

#include <Presolve/LPFormatParser.hpp>

#include <random>
#include <iostream>

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

  lp::Timer& timer = lp::Timer::getIPMInstance();
  timer.start("Over All", true);

  BOOST_LOG_TRIVIAL(info) << "Parser and Presolve Started";
  lp::Problem problem = parser.parse(options.fileName);
  BOOST_LOG_TRIVIAL(info) << "Parser and Presolve Ended";

  //   problem.normalize();

  // std::ptrdiff as index is not working
  typedef lp::IPMCholeskyLLT<
      lp::NTScalings, lp::CholmodSupernodalLLT<Eigen::SparseMatrix<double>,
                                               Eigen::Lower>> CholmodSolver;

  BOOST_LOG_TRIVIAL(info) << "Started to solve";

  lp::IPMSolver<CholmodSolver, lp::NTScalings> solver(problem);
  lp::Solution solution = solver.solve();
  BOOST_LOG_TRIVIAL(info) << solution;

  BOOST_LOG_TRIVIAL(info) << "Ended";
  timer.end("Over All");

  BOOST_LOG_TRIVIAL(info) << timer;

  lp::Timer& timer2 = lp::Timer::getADMMInstance();
  timer2.start("Over All", true);

  lp::ADMMSolver solver2(problem);
  lp::Solution solution2 = solver2.solve();

  BOOST_LOG_TRIVIAL(info) << solution2;

  timer2.end("Over All");

  BOOST_LOG_TRIVIAL(info) << timer2;
  return 0;
}
