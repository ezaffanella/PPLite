/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2026 Enea Zaffanella <enea.zaffanella@unipr.it>

This file is part of PPLite.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef pplite_clock_hh
#define pplite_clock_hh 1

#include "globals.hh"

#include <chrono>
#include <iostream>
#include <iomanip>

namespace pplite {

// This is a steady clock: it starts automatically when created;
// it can be restarted using restart(); it can be queried by calling
// elapsed_time(), returning the time elapsed since last restart.
// When printed, elapsed time is expressed in seconds, with digits
// up to microseconds (no matter what is the actual clock resolution).
class Clock {
public:
  using Steady_Clock = std::chrono::steady_clock;
  using Duration = Steady_Clock::duration;

  Clock() noexcept : saved_time(Steady_Clock::now()) {}
  // Not copyable, only moveable.
  Clock(const Clock&) = delete;
  Clock& operator=(const Clock&) = delete;
  Clock(Clock&&) = default;
  Clock& operator=(Clock&&) = default;
  ~Clock() = default;

  void restart() { saved_time = Steady_Clock::now(); }
  Duration elapsed_time() const { return Steady_Clock::now() - saved_time; }
  void print_elapsed(std::ostream& os) const { print(os, elapsed_time()); }

  static void print(std::ostream& os, Duration d) {
    using namespace std::chrono;
    // Cast duration to (floating point) seconds.
    auto secs = duration_cast<duration<double>>(d);
    // Print up to microseconds
    auto saved_flags = os.flags();
    os << std::fixed << std::setprecision(6) << secs.count();
    os.flags(saved_flags);
  }

private:
  Steady_Clock::time_point saved_time;
}; // class Clock

// A data structure collecting simple time stats.
struct Time_Stats {
  unsigned long counter; // Number of calls
  Clock::Duration time;  // Cumulative time spent in calls
  Time_Stats() : counter(0), time(0) {}

  // Adds `time_incr' time to stats,
  // increasing by 1 the number of calls.
  void incr(Clock::Duration time_incr) {
    ++counter;
    time += time_incr;
  }
}; // struct Time_Stats

} // namespace pplite

#endif // !defined(pplite_clock_hh)
