
#include <string>
#include <iostream>

// #include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/program_options.hpp>

#include <Problem.hpp>
#include <Solution.hpp>
#include <Core/Solver.hpp>

#include <Ext/NTScalings.hpp> 
#include <Ext/SuiteSparseCholeskyLLT.hpp>

#include <Presolve/LPFormatParser.hpp>

class CommandOptions {
 public:
  std::string fileName;
  int logOptions;
};

CommandOptions getOptions(int argc, char **argv) {
  namespace po = boost::program_options;

  po::options_description desc("Allowed options");
  desc.add_options()("help", "LP Format file name to solve and log options")(
      "file,f", po::value<std::string>(), "Enter file path name");

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

  return options;
}

int main(int argc, char **argv) {
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
  //   lp::Problem problem = parser.parse("/home/satya/Desktop/QiTest.lp");
  lp::Problem problem = parser.parse(options.fileName);

  lp::Solver<lp::SuiteSparseCholeskyLLT<lp::NTScalings>, lp::NTScalings> solver(
      problem);

  BOOST_LOG_TRIVIAL(info) << "Started to solve";
  lp::Solution solution = solver.solve();

  BOOST_LOG_TRIVIAL(info) << solution;

  BOOST_LOG_TRIVIAL(info) << "Ended";

  return 0;
}
