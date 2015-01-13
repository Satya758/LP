
#include <vector>
#include <algorithm>

#include <Core/Timer.hpp>

namespace lp {
int Timer::depth = 0;
int Timer::sequence = 0;
namespace internal {

std::ostream& operator<<(std::ostream& out, const Clock& clock) {
  using namespace std;

  out << setw(5) << clock.sequence << setw(5) << clock.depth << setw(10)
      << clock.iterations << setw(20) << clock.totalDuration.count() << setw(20)
      << clock.percentage;

  return out;
}
}

std::ostream& operator<<(std::ostream& out, Timer& timer) {
  using namespace std;

  timer.computePercentage();

  //   typedef std::pair<std::string, internal::Clock> FragmentPair;
  //   std::vector<FragmentPair> mapElements(timer.fragmentTime.begin(),
  //                                         timer.fragmentTime.end());
  //
  //   std::sort(mapElements.begin(), mapElements.end(),
  //             [](const FragmentPair& lhs, const FragmentPair& rhs) {
  //     return lhs.second.sequence > rhs.second.sequence;
  //   });

  out << endl << setw(50) << "Fragment Name" << setw(5) << "Seq" << setw(5)
      << "Dpt" << setw(10) << "Iters" << setw(20) << "Duration" << setw(20)
      << "Percentage %" << endl;

  out << setw(50) << timer.module + "_" + timer.overAllFragment
      << timer.fragmentTime.at(timer.overAllFragment);

  for (auto fragment : timer.fragmentTime) {
    if (fragment.first != timer.overAllFragment) {
      out << endl << setw(50) << timer.module + "_" + fragment.first
          << timer.fragmentTime.at(fragment.first);
    }
  }

  return out;
}
}