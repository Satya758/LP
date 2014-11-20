#ifndef PRESOLVE_HPP
#define PRESOLVE_HPP

#include <unordered_map>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Problem.hpp>

namespace lp {
// TODO reserve size, use guess work
// TODO For now I dont need another intermediate result I will to presolve in
// parser class itself
// But it cant be easily extendible to other formats, but for now its OK, we
// will worry about that later
// Having Intermediate structure increases copy and also might not be compatible
// with other formats out there
}
#endif  // PRESOLVE_HPP