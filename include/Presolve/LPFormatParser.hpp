#ifndef LP_FORMAT_PARSER_HPP
#define LP_FORMAT_PARSER_HPP

#define BOOST_SPIRIT_DEBUG
//#define BOOST_SPIRIT_DEBUG_OUT

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>
#include <type_traits>
#include <cmath>

#include <boost/log/trivial.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SPQRSupport>

#include <Problem.hpp>

namespace lp {
namespace parser {

namespace qi = boost::spirit::qi;
namespace fu = boost::fusion;
namespace ascii = boost::spirit::ascii;

enum class Operators {
  lt,
  le,
  gt,
  ge,
  eq,
  invalid,
  free  // Added to handle free varaibles in Bounds
};

std::ostream& operator<<(std::ostream& out, const Operators& op) {
  out << static_cast<std::underlying_type<Operators>::type>(op);
  return out;
}

class OPSymbols : public qi::symbols<char, Operators> {
 public:
  OPSymbols() {
    add("<", Operators::le)("<=", Operators::le)(">", Operators::gt)(
        ">=", Operators::ge)("=", Operators::eq)("free", Operators::free);
  }
};

class BoundSymbols : public qi::symbols<char, double> {
 public:
  BoundSymbols() {
    add("-inf", -std::numeric_limits<double>::infinity())(
        "+inf", std::numeric_limits<double>::infinity());
  }
};
// SCN means SignCoefficientName
class SCS {
 public:
  // Initialize default values
  char sign = '+';
  double value = 1.0;
  std::string name;
};

class Constraints {
 public:
  std::vector<SCS> lhs;
  Operators op;
  double rhs;
};

class Bounds {
 public:
  double lower;
  // Lower Operator
  Operators lOp = Operators::invalid;
  std::string name;
  // Upper Operator
  Operators uOp = Operators::invalid;
  double upper;
};

// Output Structure
class ParsedObject {
 public:
  std::vector<SCS> objective;
  std::vector<Constraints> constraints;
  std::vector<Bounds> bounds;
};

std::ostream& operator<<(std::ostream& out, const ParsedObject& po) {
  using namespace std;

  out << endl << "##################### PO Start" << endl;
  out << "Objective : " << po.objective.size() << endl;
  out << "Constraints: " << po.constraints.size() << endl;
  out << "Bounds: " << po.bounds.size() << endl;
  out << "##################### PO End" << endl;

  return out;
}
}
}

// Fusion does like to be in global namespace
BOOST_FUSION_ADAPT_STRUCT(lp::parser::SCS,
                          (char, sign)(double, value)(std::string, name))

BOOST_FUSION_ADAPT_STRUCT(lp::parser::Constraints,
                          (std::vector<lp::parser::SCS>,
                           lhs)(lp::parser::Operators, op)(double, rhs))

BOOST_FUSION_ADAPT_STRUCT(lp::parser::Bounds,
                          (double, lower)(lp::parser::Operators,
                                          lOp)(std::string,
                                               name)(lp::parser::Operators,
                                                     uOp)(double, upper))

BOOST_FUSION_ADAPT_STRUCT(
    lp::parser::ParsedObject,
    (std::vector<lp::parser::SCS>,
     objective)(std::vector<lp::parser::Constraints>,
                constraints)(std::vector<lp::parser::Bounds>, bounds))

namespace lp {
namespace parser {

namespace qi = boost::spirit::qi;
namespace fu = boost::fusion;
namespace ascii = boost::spirit::ascii;

// Comment and Space Skipper
template <typename Iterator>
class Skipper : public qi::grammar<Iterator> {
 public:
  Skipper() : Skipper::base_type(commentSpaceSkip) {
    using namespace qi;

    commentSpaceSkip =
        ascii::space | (lit('\\') >> *(char_ - eol) >> eol | blank);
  }

  qi::rule<Iterator> commentSpaceSkip;
};

// Grammer
template <typename Iterator, typename Skip = Skipper<Iterator>>
class LPFormatGrammar : public qi::grammar<Iterator, ParsedObject(), Skip> {
 public:
  LPFormatGrammar() : LPFormatGrammar::base_type(parsedObj) {
    using namespace qi;

    objectiveText = lit("Minimize") >> *(+alnum >> lit(":"));
    // TODO Default value is + here
    sign = char_("+") | char_("-");  // | attr("+");

    objectiveName = +(alnum[_val += _1] - "Subject To");

    // TODO Default value does not work attr(1)?
    // Had to add Default value 1 to double, when converting from my own data
    // structures
    objective = -sign >> -double_ >> objectiveName;

    constraintsText = lit("Subject To");
    constraintName = (+alnum >> lit(":"));
    // Here we dont need use difference operator as in objectiveName because
    // constraints end with double after inequality/equality operator
    constraintVarName = +(alnum[_val += _1] - ("Bounds"));

    // We are not storing constraint name
    constraintsLHS = -constraintName >> -sign >> -double_ >> constraintVarName;
    // symbols should be at class scope earlier I have kept scope at method and
    // parser started to crash at runtime
    constraints = +constraintsLHS >> opSymbols >> double_;

    boundText = lit("Bounds");
    boundVarName = +(alnum[_val += _1] - "End");
    luBounds = double_ | bSymbols;
    bounds = -luBounds >> -opSymbols >> boundVarName >> -opSymbols >> -luBounds;

    parsedObj = objectiveText >> +objective >> constraintsText >>
                +constraints >> -boundText >> *bounds;

    //     BOOST_SPIRIT_DEBUG_NODE(parsedObj);
  }

  qi::rule<Iterator, ParsedObject(), Skip> parsedObj;
  qi::rule<Iterator, Skip> objectiveText;
  qi::rule<Iterator, SCS(), Skip> objective;
  qi::rule<Iterator, std::string(), Skip> objectiveName;
  qi::rule<Iterator, char(), Skip> sign;
  qi::rule<Iterator, Skip> constraintsText;
  qi::rule<Iterator, Skip> constraintName;
  qi::rule<Iterator, std::string(), Skip> constraintVarName;
  qi::rule<Iterator, SCS(), Skip> constraintsLHS;
  qi::rule<Iterator, Constraints(), Skip> constraints;
  qi::rule<Iterator, Skip> boundText;
  qi::rule<Iterator, std::string(), Skip> boundVarName;
  // lower upper Bounds
  qi::rule<Iterator, double(), Skip> luBounds;
  qi::rule<Iterator, Bounds(), Skip> bounds;
  // Need to be member variable otherwise parser crashes if it is defined in
  // local scope
  OPSymbols opSymbols;
  BoundSymbols bSymbols;
};

class LPFormatParser {
 public:
  Problem parse(const std::string& fileName) const {
    using namespace qi;

    std::ifstream fileStream(fileName, std::ios_base::in);

    if (!fileStream) {
      BOOST_LOG_TRIVIAL(info) << "Failed to read file";
      // FIXME Had to return problem, but problem is not created, may be raise
      // exception

      // return;
    }

    std::string fileStorage;

    fileStream.unsetf(std::ios::skipws);

    std::copy(std::istream_iterator<char>(fileStream),
              std::istream_iterator<char>(), std::back_inserter(fileStorage));

    std::string::const_iterator firstIterator = fileStorage.begin();
    std::string::const_iterator lastIterator = fileStorage.end();

    LPFormatGrammar<std::string::const_iterator> lpGrammar;
    Skipper<std::string::const_iterator> skipper;

    ParsedObject obj;

    bool isParsingSuccess =
        phrase_parse(firstIterator, lastIterator, lpGrammar, skipper, obj);

    if (!isParsingSuccess) {
      BOOST_LOG_TRIVIAL(info) << "Failed";
    }

    BOOST_LOG_TRIVIAL(info) << obj;

    Problem problem = getProblem(obj);

    BOOST_LOG_TRIVIAL(info) << "After presolve: " << problem;

    return problem;
  }

 private:
  typedef Eigen::Triplet<double> CT;
  /**
   * Add default values when data is copied to problem
   *
   * Idea is
   * 1. First push objective along with default values and proper signs
   * 2. Store names of all variables in unordered_map to have O(1) access along
   *    with index, e.g. {Name, colIndex>}, colIndex is just incremented int for
   *    each name
   * 3. Push the constraints into Triplet along with signs and default values,
   *    based on operator convert into standard form required by Problem, fill
   *    RHS h vector as well constraint by constrant, If it is eq operator
   *    duplicate the entries with le and ge operators.
   *    e.g. g1 + g2 = h1 results in g1 + g2 >= h1; g1 + g2 <= h2 (ofcourse h1
   *    == h2)
   * 4. Finally check the bounds and add them to consraints based on operator
   * 5. Position of constraint coeff (colIndex) can obtained from unordered_map
   * 6. First data is copied to temporary objects (Vectors, Vector<Triplets>)
   *    so that we know the size, and create Problem object with known sizes
   *
   */
  Problem getProblem(const ParsedObject& parsedObj) const {

    std::vector<double> c;
    //     c.reserve(parsedObj.objective.size());
    // Approximate size, worst case actual size will be double if there are
    // equality constraints
    // Reseve some size
    std::vector<double> h;
    //     h.reserve(2 * parsedObj.constraints.size() + 2 *
    // parsedObj.bounds.size());
    std::unordered_map<std::string, int> nameMap;
    //     nameMap.reserve(parsedObj.objective.size());

    int index = 0;
    for (const SCS& scs : parsedObj.objective) {

      switch (scs.sign) {
        case '+':
          c.push_back(scs.value);
          break;
        case '-':
          c.push_back(-scs.value);
          break;
      }

      nameMap.insert({scs.name, index});

      ++index;
    }

    BOOST_LOG_TRIVIAL(info) << "Got the names " << nameMap.size();
    std::vector<CT> constraintTriplets =
        getConstraintTriplets(parsedObj, nameMap, c, h);

    BOOST_LOG_TRIVIAL(info) << "Got the constraintTriplets "
                            << constraintTriplets.size();

    addBounds(parsedObj, nameMap, constraintTriplets, c, h);
    BOOST_LOG_TRIVIAL(info) << "Got the bounds " << constraintTriplets.size();
    return getProblem(constraintTriplets, h, c);
  }

  // Convert to problem
  Problem getProblem(const std::vector<CT>& constraintTriplets,
                     const std::vector<double>& h,
                     const std::vector<double>& c) const {

    Problem problem(h.size(), 0, c.size());
    {
      // Make objective
      int index = 0;
      for (const double value : c) {
        problem.c(index) = value;
        ++index;
      }
    }

    {  // Make inequality rhs
      int index = 0;
      for (const double value : h) {
        problem.h(index) = value;
        ++index;
      }
    }

    {
      // Make inequality consraints lhs
      problem.G.setFromTriplets(constraintTriplets.begin(),
                                constraintTriplets.end());
    }

    BOOST_LOG_TRIVIAL(info) << "Before QR: " << problem;

    removeRedundantRows(problem);

    return problem;
  }

  /**
   *
   * Full row rank is required for linear programming, as all proofs are
   * considered with full row rank of a constrained matrix
   *
   * Determination of full row rank
   * 1. Using QR (SPQR) from suitesparse as it has more chance of getting
   *    parallelized (in documentaion), Eigen has also QR implemented I do not
   *    know the performace comparision between Eigen QR and suitesparse QR
   *    (TODO Get the benchmarks for both functions)
   * 2. QR is used to find the column rank of matrix (I did not find any way to
   *    do it to get for row rank ).
   * 3. First matrix is transposed to get rows to columns (we should copy
   *    constraint
   *    matrix to avoid changing original ), the QR is called to get rank and
   *    permutation matrix.
   * 4. Permuation is multiplied from left hand side to bring up all independent
   *    rows to top, and then Rank is used to resize (shrink) the row size of
   *    matrix, which will also remove all dependent rows.
   * 5. Same is done for RHS
   * 6. After this original sizes given to Problem constructor are not any more
   *valid
   * FIXME Currently using 4.2.1 as rank is not return in cholmod_cc (common)
   * structure for some reason, but m_rank is Eigen class has rank value, we
   *have
   * to inherit that class and modify getRank method to get actual rank to
   *delete
   * dependent columns
   *
   * As equality for we don't delete rows, we check here for dependent columns
   */
  void removeRedundantRows(Problem& problem) const {
    Eigen::SPQR<Eigen::SparseMatrix<double>> qrSolver;
    // TODO Use constants from common header
    qrSolver.setPivotThreshold(1.0e-10);

    qrSolver.compute(problem.G);
    BOOST_LOG_TRIVIAL(info) << "Rank: " << qrSolver.rank();
    // TODO When I use same matrix on both side resultant matrix is full of
    // zeros, this is eigen way of doing, have to find reason and better way to
    // do it
    BOOST_LOG_TRIVIAL(info) << "Before Copy: ";
//     Eigen::SparseMatrix<double> copyOFG(problem.G);
    problem.G = problem.G * qrSolver.colsPermutation();
    BOOST_LOG_TRIVIAL(info) << "After Perm " << problem.G.cols();
    problem.G.conservativeResize(problem.G.rows(), qrSolver.rank());
    BOOST_LOG_TRIVIAL(info) << "After Perm 2";
    problem.c = qrSolver.colsPermutation() * problem.c;
    // FIXME We have to keep track of removed columns, how do we compute final
    // value in post solve?
    BOOST_LOG_TRIVIAL(info) << "After Perm 3";
    problem.c.conservativeResize(qrSolver.rank());
    BOOST_LOG_TRIVIAL(info) << "After Perm 4";
  }

  // FIXME Raw loops refactor the code!!
  std::vector<CT> getConstraintTriplets(
      const ParsedObject& parsedObj,
      std::unordered_map<std::string, int>& nameMap, std::vector<double>& c,
      std::vector<double>& h) const {

    std::vector<CT> constraintTriplets;
    // Number of elements in constraints are equal to row * columns, in this
    // case rows are in constraints and columns are in objective, to
    // consider additional added rows (because of eqaulity and bounds)
    // we just double the size of row * columns
    // as more is better than less
    //     constraintTriplets.reserve(2 * parsedObj.objective.size() *
    //                                parsedObj.constraints.size());

    int rowIndex = 0;
    for (const Constraints& ct : parsedObj.constraints) {
      int signModifier;
      bool isEqual = false;

      // There are no lt and gt?
      switch (ct.op) {
        case Operators::le:
          signModifier = 1;
          break;
        case Operators::ge:
          signModifier = -1;
          break;
        case Operators::eq:
          isEqual = true;
          break;
      }

      if (!isEqual) {
        for (const SCS& scs : ct.lhs) {
          int colIndex = getColIndex(scs.name, nameMap, c);

          switch (scs.sign) {
            case '+':
              constraintTriplets.push_back(
                  CT(rowIndex, colIndex, signModifier * scs.value));
              break;
            case '-':
              constraintTriplets.push_back(
                  CT(rowIndex, colIndex, -(signModifier * scs.value)));
              break;
          }
        }
        h.push_back(signModifier * ct.rhs);

      } else if (isEqual) {
        int doubleInserter = 0;
        // For greater than or equal to Sign
        int geSign = 1;
        while (doubleInserter < 2) {
          if (doubleInserter == 1) {
            geSign = -1;
            ++rowIndex;  // Increment row for duplicate entry in next line
          }

          for (const SCS& scs : ct.lhs) {
            int colIndex = getColIndex(scs.name, nameMap, c);

            switch (scs.sign) {
              case '+':
                constraintTriplets.push_back(
                    CT(rowIndex, colIndex, geSign * scs.value));
                break;
              case '-':
                constraintTriplets.push_back(
                    CT(rowIndex, colIndex, -(geSign * scs.value)));
                break;
            }
          }
          h.push_back(geSign * ct.rhs);
          ++doubleInserter;
        }
      }

      // BOOST_LOG_TRIVIAL(info) << "Added row: " << rowIndex;
      ++rowIndex;
    }

    return constraintTriplets;
  }

  // Returns column index to insert into
  // We need objective as we have to add zero to objective for varaible that are
  // not found in objective definition
  int getColIndex(const std::string& name,
                  std::unordered_map<std::string, int>& nameMap,
                  std::vector<double>& c) const {
    int colIndex;

    std::unordered_map<std::string, int>::const_iterator elem =
        nameMap.find(name);

    if (elem == nameMap.end()) {
      c.push_back(0);
      // We will add it to the end
      colIndex = c.size() - 1;
      // Insert for next reference
      nameMap.insert({name, colIndex});
    } else {
      colIndex = elem->second;
    }

    return colIndex;
  }

  // // FIXME Raw loops refactor the code!!
  void addBounds(const ParsedObject& parsedObj,
                 std::unordered_map<std::string, int>& nameMap,
                 std::vector<CT>& constraintTriplets, std::vector<double>& c,
                 std::vector<double>& h) const {
    // rowIndex is incremented for every addition, so there is so much of
    // FIXME duplication

    int rowIndex = h.size() - 1;
    for (const Bounds& bounds : parsedObj.bounds) {
      int colIndex = getColIndex(bounds.name, nameMap, c);

      if (bounds.uOp != Operators::invalid && !std::isinf(bounds.upper)) {
        // There is no lt and gt
        switch (bounds.uOp) {
          // Push as is with coefficient equal to 1, no requirement of sign
          // change
          case Operators::le:
            ++rowIndex;
            constraintTriplets.push_back(CT(rowIndex, colIndex, 1));
            h.push_back(bounds.upper);
            break;
          // Sign change
          case Operators::ge:
            ++rowIndex;
            constraintTriplets.push_back(CT(rowIndex, colIndex, -1));
            h.push_back(-bounds.upper);
            break;
          // Insert two records with different signs
          case Operators::eq:
            ++rowIndex;
            constraintTriplets.push_back(CT(rowIndex, colIndex, 1));
            h.push_back(bounds.upper);
            ++rowIndex;
            constraintTriplets.push_back(CT(rowIndex, colIndex, -1));
            h.push_back(-bounds.upper);
            break;
        }
      }
      if (bounds.lOp != Operators::invalid && !std::isinf(bounds.lower)) {
        switch (bounds.lOp) {
          // Push as is with coefficient equal to 1, no requirement of sign
          // change
          case Operators::ge:
            ++rowIndex;
            constraintTriplets.push_back(CT(rowIndex, colIndex, 1));
            h.push_back(bounds.lower);
            break;
          // Sign change
          case Operators::le:
            ++rowIndex;
            constraintTriplets.push_back(CT(rowIndex, colIndex, -1));
            h.push_back(-bounds.lower);
            break;
          // This scenario will not be encountered in our parser
          case Operators::eq:
            // DO nothing
            break;
        }
      }
    }
  }
};
}
}

#endif  // LP_FORMAT_PARSER_HPP