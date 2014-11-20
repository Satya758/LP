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

#include <boost/log/trivial.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

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
  eq
};

std::ostream& operator<<(std::ostream& out, const Operators& op) {
  out << static_cast<std::underlying_type<Operators>::type>(op);
  return out;
}

class OPSymbols : public qi::symbols<char, Operators> {
 public:
  OPSymbols() {
    add("<", Operators::le)("<=", Operators::le)(">", Operators::gt)(
        ">=", Operators::ge)("=", Operators::eq);
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
  char sign;
  double value;
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
  Operators lOp;
  std::string name;
  // Upper Operator
  Operators uOp;
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
                          (double, lower) (lp::parser::Operators, lOp) (std::string, name) (lp::parser::Operators, uOp) (double, upper))

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

    //BOOST_SPIRIT_DEBUG_NODE(parsedObj);
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
  void parse(const std::string& fileName) const {
    using namespace qi;

    std::ifstream fileStream(fileName, std::ios_base::in);

    if (!fileStream) {
      BOOST_LOG_TRIVIAL(info) << "Failed to read file";

      return;
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
  }
};
}
}

#endif  // LP_FORMAT_PARSER_HPP