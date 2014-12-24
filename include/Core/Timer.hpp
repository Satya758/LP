#ifndef TIMER_HPP
#define TIMER_HPP

#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <chrono>

namespace lp {

enum class Fragments {
  Parsing,
  RemoveRedundantCols,
  Presolve,
  // Cholesky
  CholeskyAnalyze,
  CholeskyFactorize,
  CholeskySolve,
  // PSD Matrix
  PSDMatrix,
  // Residual
  ResidualComputation,
  InitialPoint,
  SolverStateComputation,
  ScalingsCompute,
  StepSizeCompute,
  AffineDirection,
  CombinedDirection,
  DirectionUpdate,
  OverAll
};

namespace internal {

class Clock {
 public:
  Clock() {}

  void setStart() {
    if (started && !ended) {
      throw std::runtime_error("Fragment already started but not ended");
    }
    started = true;

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

 private:
  std::chrono::time_point<std::chrono::system_clock> start;
  std::chrono::time_point<std::chrono::system_clock> end;
  std::chrono::duration<double> totalDuration;
  double percentage;

  bool started = false;
  bool ended = true;

  friend std::ostream& operator<<(std::ostream& out, const Clock& clock);
};

std::ostream& operator<<(std::ostream& out, const Clock& clock) {
  using namespace std;

  out << setw(20) << "Duration: " << setw(20) << clock.totalDuration.count()
      << setw(20) << " Percentage % : " << setw(20) << clock.percentage;

  return out;
}
}

class Timer {
 public:
  static Timer& getInstance() {
    static Timer timer;

    return timer;
  }

  void start(Fragments fragment) { fragmentTime.at(fragment).setStart(); }

  void end(Fragments fragment) { fragmentTime.at(fragment).setEnd(); }

 private:
  Timer(const Timer& other) = delete;
  Timer& operator=(const Timer& other) = delete;

  Timer() {
    fragmentTime.emplace(std::make_pair(Fragments::Parsing, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::RemoveRedundantCols, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::Presolve, internal::Clock()));

    fragmentTime.emplace(
        std::make_pair(Fragments::CholeskyAnalyze, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::CholeskyFactorize, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::CholeskySolve, internal::Clock()));

    fragmentTime.emplace(
        std::make_pair(Fragments::PSDMatrix, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::ResidualComputation, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::InitialPoint, internal::Clock()));

    fragmentTime.emplace(
        std::make_pair(Fragments::SolverStateComputation, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::ScalingsCompute, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::StepSizeCompute, internal::Clock()));

    fragmentTime.emplace(
        std::make_pair(Fragments::AffineDirection, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::CombinedDirection, internal::Clock()));
    fragmentTime.emplace(
        std::make_pair(Fragments::DirectionUpdate, internal::Clock()));

    fragmentTime.emplace(std::make_pair(Fragments::OverAll, internal::Clock()));
  }

  std::map<Fragments, internal::Clock> fragmentTime;

  void computePercentage() {
    internal::Clock& overAllClock = fragmentTime.at(Fragments::OverAll);

    for (auto& fragmentPair : fragmentTime) {
      fragmentPair.second.computePercentage(overAllClock);
    }
  }

  friend std::ostream& operator<<(std::ostream& out, Timer& timer);
};

std::ostream& operator<<(std::ostream& out, Timer& timer) {
  using namespace std;

  timer.computePercentage();

  out << endl << setw(40)
      << "Over All: " << timer.fragmentTime.at(Fragments::OverAll);

  out << endl << setw(40)
      << "Parsing: " << timer.fragmentTime.at(Fragments::Parsing);
  out << endl << setw(40) << "Remove Redundant Cols: "
      << timer.fragmentTime.at(Fragments::RemoveRedundantCols);

  out << endl << setw(40) << "Cholesky Analyze: "
      << timer.fragmentTime.at(Fragments::CholeskyAnalyze);
  out << endl << setw(40) << "Cholesky Factorize: "
      << timer.fragmentTime.at(Fragments::CholeskyFactorize);
  out << endl << setw(40)
      << "Cholesky Solve: " << timer.fragmentTime.at(Fragments::CholeskySolve);
  out << endl << setw(40)
      << "PSDMatrix: " << timer.fragmentTime.at(Fragments::PSDMatrix);

  out << endl << setw(40) << "Residual Computation: "
      << timer.fragmentTime.at(Fragments::ResidualComputation);
  out << endl << setw(40)
      << "Initial Point: " << timer.fragmentTime.at(Fragments::InitialPoint);
  out << endl << setw(40) << "SolverState Computation: "
      << timer.fragmentTime.at(Fragments::SolverStateComputation);

  out << endl << setw(40) << "Scalings Compute: "
      << timer.fragmentTime.at(Fragments::ScalingsCompute);
  out << endl << setw(40) << "StepSize Compute: "
      << timer.fragmentTime.at(Fragments::StepSizeCompute);

  out << endl << setw(40) << "Affine Direction: "
      << timer.fragmentTime.at(Fragments::AffineDirection);
  out << endl << setw(40) << "Combined Direction: "
      << timer.fragmentTime.at(Fragments::CombinedDirection);
  out << endl << setw(40) << "Direction Update: "
      << timer.fragmentTime.at(Fragments::DirectionUpdate);

  return out;
}
}
#endif  // TIMER_HPP