#ifndef TIMER_HPP
#define TIMER_HPP

#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <chrono>

namespace lp {

namespace internal {

class Clock {
 public:
  Clock(const int sequence, const int depth)
      : sequence(sequence), depth(depth) {}

  void setStart() {
    if (started && !ended) {
      throw std::runtime_error("Fragment already started but not ended");
    }
    started = true;
    ++iterations;
    start = std::chrono::system_clock::now();
  }

  void setEnd() {
    if (!started && ended) {
      throw std::runtime_error("Fragment not started but trying to end");
    }
    ended = true;

    end = std::chrono::system_clock::now();

    std::chrono::duration<double> duration = end - start;

    totalDuration += duration;

    started = false;
    ended = true;
  }

  void computePercentage(Clock& overAllClock) {
    double totalTime = overAllClock.totalDuration.count();
    percentage = 100 * totalDuration.count() / totalTime;
  }

  const int sequence;
  const int depth;

 private:
  std::chrono::time_point<std::chrono::system_clock> start;
  std::chrono::time_point<std::chrono::system_clock> end;
  std::chrono::duration<double> totalDuration;
  double percentage;

  bool started = false;
  bool ended = true;

  int iterations = 0;

  friend std::ostream& operator<<(std::ostream& out, const Clock& clock);
};

std::ostream& operator<<(std::ostream& out, const Clock& clock);
}

class Timer {
 public:
  static Timer& getIPMInstance() {
    static Timer ipmTimer("IPM");

    return ipmTimer;
  }

  static Timer& getADMMInstance() {
    static Timer admmTimer("ADMM");

    return admmTimer;
  }

  void start(std::string fragment, bool isOverAll) {
    try {
      fragmentTime.at(fragment).setStart();
      ++depth;
    }
    catch (const std::out_of_range& e) {
      ++depth;
      ++sequence;
      fragmentTime.emplace(fragment, internal::Clock(sequence, depth));
      fragmentTime.at(fragment).setStart();
    }

    if (isOverAll) {
      overAllFragment = fragment;
    }
  }

  void start(std::string fragment) { start(fragment, false); }

  void end(std::string fragment) {
    --depth;
    fragmentTime.at(fragment).setEnd();
  }

 private:
  Timer(const Timer& other) = delete;
  Timer& operator=(const Timer& other) = delete;
  // TODO Not thread safe
  static int depth;
  static int sequence;

  Timer(const std::string module) : module(module) {}

  std::map<std::string, internal::Clock> fragmentTime;
  std::string overAllFragment;

  void computePercentage() {
    internal::Clock& overAllClock = fragmentTime.at(overAllFragment);

    for (auto& fragmentPair : fragmentTime) {
      fragmentPair.second.computePercentage(overAllClock);
    }
  }

  const std::string module;

  friend std::ostream& operator<<(std::ostream& out, Timer& timer);
};
}
#endif  // TIMER_HPP
