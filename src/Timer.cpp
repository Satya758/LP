
#include <Core/Timer.hpp>

namespace lp {
namespace internal {
std::ostream& operator<<(std::ostream& out, const Clock& clock) {
  using namespace std;

  out << setw(20) << "Duration: " << setw(20) << clock.totalDuration.count()
      << setw(20) << " Percentage % : " << setw(20) << clock.percentage;

  return out;
}
}

std::ostream& operator<<(std::ostream& out, Timer& timer) {
  using namespace std;

  timer.computePercentage();

  out << endl << setw(40) << timer.overAllFragment
      << timer.fragmentTime.at(timer.overAllFragment);

  for (auto fragment : timer.fragmentTime) {
    if (fragment.first != timer.overAllFragment) {
      out << endl << setw(40) << fragment.first
          << timer.fragmentTime.at(fragment.first);
    }
  }

  return out;
}
}