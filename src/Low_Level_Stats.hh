/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#ifndef pplite_Low_Level_Stats_hh
#define pplite_Low_Level_Stats_hh 1

#include "globals.hh"
#include "clock.hh"

#include <iostream>

namespace pplite {

enum class LLOp {
  CONV_C2G, CONV_G2C, INCR_C2G, INCR_G2C, NUM_OPS
};

constexpr auto llop_size = static_cast<unsigned>(LLOp::NUM_OPS);

// The array of operation names (indexed by LLOp).
// Careful: these must match the ordering of LLOp enumeration values.
constexpr const char* llop_names[llop_size] = {
  "conv_con_to_gen", "conv_gen_to_con",
  "incr_con_to_gen", "incr_gen_to_con",
};

// A data structure collecting stats for each LLOp.
struct LLOp_Stats {
  // Allow to dynamically enable/disable counters.
  static PPLITE_TLS bool are_enabled;
  // These are the global counters.
  static PPLITE_TLS Counter scalar_prod;
  static PPLITE_TLS Counter linear_comb;
  static PPLITE_TLS Counter sat_count;
  static PPLITE_TLS Counter sat_incl;
  static PPLITE_TLS Counter sat_union;

  static void enable_counters() { are_enabled = true; }
  static void disable_counters() { are_enabled = false; }
  static void bump_scalar_prod() { if (are_enabled) ++scalar_prod; }
  static void bump_linear_comb() { if (are_enabled) ++linear_comb; }
  static void bump_sat_count() { if (are_enabled) ++sat_count; }
  static void bump_sat_incl() { if (are_enabled) ++sat_incl; }
  static void bump_sat_union() { if (are_enabled) ++sat_union; }

  static void print_stats(std::ostream& os) {
    os << " SP " << scalar_prod;
    os << " LC " << linear_comb;
    os << " SC " << sat_count;
    os << " SI " << sat_incl;
    os << " SU " << sat_union;
  }
  static void reset_stats() {
    scalar_prod = 0;
    linear_comb = 0;
    sat_count = 0;
    sat_incl = 0;
    sat_union = 0;
  }

  // This collects statistics info for a single operation.
  struct Stats {
    Counter calls;         // Number of calls
    Clock::Duration time;  // Cumulative time spent in calls
    Counter scalar_prod;   // Number of scalar products.
    Counter linear_comb;   // Number of linear combines.
    Counter sat_count;     // Number of sat population counts.
    Counter sat_incl;      // Number of sat inclusion tests.
    Counter sat_union;     // Number of sat rows unions.
    Stats() : calls(0), time(0),
              scalar_prod(0), linear_comb(0),
              sat_count(0), sat_incl(0), sat_union(0) {}
  };

  // The collected statistics.
  Stats stats[llop_size];

  // Whether or not the destructor should output statistics info.
  bool noisy_dtor_;
  // Flag to dynamically enable/disable the tracking of statistics.
  bool track;

  explicit LLOp_Stats(bool noisy_dtor = false)
    : noisy_dtor_(noisy_dtor), track(true) {}

  // Clears all statistics.
  void reset();

  // Adds `incr_time' to stats for `op',
  // increasing by 1 the number of calls.
  void incr(LLOp op, Clock::Duration incr_time,
            Counter incr_sp, Counter incr_lc,
            Counter incr_sc, Counter incr_si, Counter incr_su) {
    // Do nothing if we are not tracking stats.
    if (!track)
      return;
    Stats& s = stats[static_cast<unsigned>(op)];
    ++s.calls;
    // Time
    s.time += incr_time;
    // Low level stats.
    s.scalar_prod += incr_sp;
    s.linear_comb += incr_lc;
    s.sat_count += incr_sc;
    s.sat_incl += incr_si;
    s.sat_union += incr_su;
  }

  void dump_op(std::ostream& os, LLOp op, bool dump_times = false) const;
  void dump(std::ostream& os, bool dump_times = false) const;
  ~LLOp_Stats();
};

extern PPLITE_TLS LLOp_Stats llop_stats;

class LLOp_Clock {
public:
  LLOp_Clock(LLOp op, bool noisy = false)
    : clock(), op_(op), noisy_(noisy) { init(); }
  ~LLOp_Clock() { stop_clock(); }
  // Neither copyable nor moveable.
  LLOp_Clock(const LLOp_Clock&) = delete;
  LLOp_Clock& operator=(const LLOp_Clock&) = delete;
  LLOp_Clock(LLOp_Clock&&) = delete;
  LLOp_Clock& operator=(LLOp_Clock&&) = delete;

  void init();
  void stop_clock() const;
  void elapsed_ops(Counter& incr_sp, Counter& incr_sat) const;

private:
  Clock clock;
  Counter saved_scalar_prod;
  Counter saved_linear_comb;
  Counter saved_sat_count;
  Counter saved_sat_incl;
  Counter saved_sat_union;
  const LLOp op_;
  const bool noisy_;
}; // class LLOp_Clock

} // namespace pplite

#endif // !defined(pplite_Low_Level_Stats_hh)
