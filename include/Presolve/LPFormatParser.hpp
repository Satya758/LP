#ifndef LP_FORMAT_PARSER_HPP
#define LP_FORMAT_PARSER_HPP

#include <vector>
#include <unordered_map>

#include <Problem.hpp>
namespace lp {
namespace parser {

/**
 * Forward declaration for output structure of raw parser, ofcourse methods does
 * not have to be private members (what about accessing private varibles
 * check!!!)
 */
class ParsedObject;

/**
 * Coverts LP Format into Problem definition, like removing equality
 * constraints, converting greater than to less than inequalities, bounds to
 * constraints, extending objective with missing varaibles
 *
 */
class LPFormatParser {
 public:
  Problem parse(const std::string& fileName) const;

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
  Problem getProblem(ParsedObject& parsedObj) const;

  // Convert to problem
  Problem getProblem(const std::vector<CT>& constraintTriplets,
                     const std::vector<double>& h,
                     const std::vector<double>& c) const;

  /**
   *
   * Full column rank is required for linear programming, as all proofs are
   * considered with full column rank of a constrained matrix
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
  void removeRedundantColumns(Problem& problem) const;

  std::vector<CT> getConstraintTriplets(
      const ParsedObject& parsedObj,
      std::unordered_map<std::string, int>& nameMap, std::vector<double>& c,
      std::vector<double>& h) const;

  // Returns column index to insert into
  // We need objective as we have to add zero to objective for varaible that are
  // not found in objective definition
  int getColIndex(const std::string& name,
                  std::unordered_map<std::string, int>& nameMap,
                  std::vector<double>& c) const;

  /**
   * In LP format if varaible is not specified in the bounds
   * or not bounded in constraints, its considered to be in positive orthant
   * i.e. by default x >= 0
   *
   * By this time all names are in nameMap, so missing names in bounds will be
   * added to bounds, here nameMap is source and bounds is target
   */
  void addBounds(const ParsedObject& parsedObj,
                 std::unordered_map<std::string, int>& nameMap,
                 std::vector<CT>& constraintTriplets, std::vector<double>& c,
                 std::vector<double>& h) const;

  /**
   * Missing bounds in LP format should be positive orthant
   * Unfortunately for us, bounds can also be represted in constraints
   * so we have to check row signleton constriants before adding missing
   * pieces to bounds
   */
  void addMissingBounds(const std::unordered_map<std::string, int>& nameMap,
                        ParsedObject& parsedObject) const;

  /**
   * Returns row singletons used to find missing bounds
   * As there will be less or no singleton constraints, we can store it in
   * vector for search later
   */
  std::vector<std::string> getRowSingletonConstraints(
      const ParsedObject& parsedObject) const;

  // FIXME Remove, just for testing
  void printInCVXSparseForm(const Problem& problem1) const;
};
}
}

#endif  // LP_FORMAT_PARSER_HPP