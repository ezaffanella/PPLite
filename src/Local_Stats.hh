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

#ifndef pplite_Local_Stats_hh
#define pplite_Local_Stats_hh 1

#include "globals.hh"
#include "clock.hh"
#include <string>

namespace pplite {

// Local_Stats and Local_Clock are used to simplify collection of
// execution time statistics for a single chunk of code.
// Example:
// void my_fun(...) {
//   static PPLITE_TLS Local_Stats my_stats("my custom name");
//   Local_Clock my_clock(my_stats);
//   <code whose efficiency is to be measured>
// }
// The first time my_fun is called, my_stats gets initialized;
// on each call, the ctor/dtor of my_clock will upgrade the statistics;
// at the end of the program, the destructor of my_stats will print
// the overall statistics.
//
// NOTE: when the library is configured to be thread safe, then macro
// PPLITE_TLS expands to "thread_local". Hence, each thread will compute
// and print its own stats (you may also want to customize the stat name,
// e.g., by using the thread id).

struct Local_Stats {
  std::string name;
  Time_Stats time_stats;

  explicit Local_Stats(std::string s) : name(std::move(s)), time_stats() {}
  ~Local_Stats() { dump(std::cerr); }

  // Neither copyable nor moveable.
  Local_Stats(const Local_Stats&) = delete;
  Local_Stats& operator=(const Local_Stats&) = delete;
  Local_Stats(Local_Stats&&) = delete;
  Local_Stats& operator=(Local_Stats&&) = delete;

  void incr(Clock::Duration time_incr) { time_stats.incr(time_incr); }

  void dump(std::ostream& os) const {
    using namespace std;
    os << name << " ";
    const Time_Stats& ts = time_stats;
    os << ts.counter << " ";
    Clock::print(os, ts.time);
    os << std::endl;
  }

}; // struct Local_Stats

struct Local_Clock {
  Clock clock;
  Local_Stats& local_stats;
  bool noisy;

  explicit Local_Clock(Local_Stats& ls, bool noisy_dtor = false)
    : local_stats(ls), noisy(noisy_dtor) {}
  ~Local_Clock() { stop_clock(); }

  // Neither copyable nor moveable.
  Local_Clock(const Local_Clock&) = delete;
  Local_Clock& operator=(const Local_Clock&) = delete;
  Local_Clock(Local_Clock&&) = delete;
  Local_Clock& operator=(Local_Clock&&) = delete;

  void stop_clock() const {
    auto elapsed = clock.elapsed_time();
    local_stats.incr(elapsed);
    if (noisy) {
      std::cerr << local_stats.name << ": ";
      std::cerr << "call number " << local_stats.time_stats.counter << ": ";
      Clock::print(std::cerr, elapsed);
      std::cerr << std::endl;
    }
  }
}; // struct Local_Clock

} // namespace pplite

#endif // !defined(pplite_Local_Stats_hh)
