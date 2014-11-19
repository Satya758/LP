#ifndef LP_FORMAT_PARSER_HPP
#define LP_FORMAT_PARSER_HPP

#define BOOST_SPIRIT_DEBUG
//#define BOOST_SPIRIT_DEBUG_OUT

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>

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

class Objective {
public:
  char sign;
  double value;
  std::string name;
};
// Output Structure
class ParsedObject {
 public:
  std::vector<Objective> objective;
};

std::ostream& operator<<(std::ostream& out, const ParsedObject& po) {
  using namespace std;

  out << endl << "##################### PO Start" << endl;
  out << "Objective : " << po.objective.size() << endl;
  out << "##################### PO End" << endl;

  return out;
}
}
}

// Fusion does like to be in global namespace
BOOST_FUSION_ADAPT_STRUCT(lp::parser::Objective,
			  (char, sign)
                          (double, value)
                           (std::string,name))

BOOST_FUSION_ADAPT_STRUCT(lp::parser::ParsedObject,
                          (std::vector<lp::parser::Objective>,
                           objective))

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

    commentSpaceSkip = ascii::space | (lit('\\') >> *(char_ - eol) >> eol | blank );
  }

  qi::rule<Iterator> commentSpaceSkip;
};

void print(ParsedObject i) {
  //BOOST_LOG_TRIVIAL(info) << i;
}
// Grammer
template <typename Iterator, typename Skip = Skipper<Iterator>>
class LPFormatGrammar
    : public qi::grammar<Iterator, ParsedObject(), Skip> {
 public:
  LPFormatGrammar() : LPFormatGrammar::base_type(parsedObj) {
    using namespace qi;

    objectiveText = lit("Minimize") >> *(+alpha >> lit(":"));
    sign = char_("+") | char_("-");

    name = +alnum[_val += _1];
    objective = -sign >> double_  >> name;

    parsedObj = objectiveText >> +objective - lit("Subject To");

    BOOST_SPIRIT_DEBUG_NODE(parsedObj);
  }

  qi::rule<Iterator, ParsedObject(), Skip> parsedObj;
  qi::rule<Iterator, Skip> objectiveText;
  qi::rule<Iterator, Objective(), Skip> objective;
  qi::rule<Iterator, std::string(), Skip> name;
  qi::rule<Iterator, char(), Skip> sign;
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