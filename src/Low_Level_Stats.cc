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

#include "pplite-config.h"
#include "Low_Level_Stats.hh"

#include "utils.hh"
#include <iostream>
#include <utility>

#if PPLITE_CONVERSION_TIME_STATS
PPLITE_TLS pplite::LLOp_Stats pplite::llop_stats(true);
#else
PPLITE_TLS pplite::LLOp_Stats pplite::llop_stats(false);
#endif

PPLITE_TLS bool pplite::LLOp_Stats::are_enabled = true;
PPLITE_TLS pplite::Counter pplite::LLOp_Stats::scalar_prod = 0;
PPLITE_TLS pplite::Counter pplite::LLOp_Stats::linear_comb = 0;
PPLITE_TLS pplite::Counter pplite::LLOp_Stats::sat_count = 0;
PPLITE_TLS pplite::Counter pplite::LLOp_Stats::sat_incl = 0;
PPLITE_TLS pplite::Counter pplite::LLOp_Stats::sat_union = 0;

namespace pplite {

void
LLOp_Stats::dump_op(std::ostream& os, LLOp op, bool dump_times) const {
  auto op_val = static_cast<unsigned>(op);
  // Do not dump if operation was never called.
  if (stats[op_val].calls == 0)
    return;
  using namespace std;
  const char fill_char = os.fill();
  os << "\n" << llop_names[op_val] << " ";
  const Stats& st = stats[op_val];
  os << st.calls;
  if (dump_times) {
    os << " ";
    Clock::print(os, st.time);
  }
  os.fill(fill_char);
#if PPLITE_LOW_LEVEL_COUNTERS
  os << " SP " << st.scalar_prod;
  os << " LC " << st.linear_comb;
  os << " SC " << st.sat_count;
  os << " SI " << st.sat_incl;
  os << " SU " << st.sat_union;
#endif
  os << flush;
}

void
LLOp_Stats::dump(std::ostream& os, bool dump_times) const {
  os << "=== Dumping time stats ===\n";
  for (auto i : range(llop_size))
    dump_op(os, static_cast<LLOp>(i), dump_times);
  os << std::endl;
}

void
LLOp_Stats::reset() {
  for (auto& stat : stats) {
    stat.calls = 0;
    stat.time = Clock::Duration::zero();
    stat.scalar_prod = 0;
    stat.linear_comb = 0;
    stat.sat_count = 0;
    stat.sat_incl = 0;
    stat.sat_union = 0;
  }
}

LLOp_Stats::~LLOp_Stats() {
  if (noisy_dtor_)
    dump(std::cerr, true);
}

void LLOp_Clock::init() {
#if PPLITE_LOW_LEVEL_COUNTERS
  saved_scalar_prod = LLOp_Stats::scalar_prod;
  saved_linear_comb = LLOp_Stats::linear_comb;
  saved_sat_count = LLOp_Stats::sat_count;
  saved_sat_incl = LLOp_Stats::sat_incl;
  saved_sat_union = LLOp_Stats::sat_union;
#endif
}

void LLOp_Clock::stop_clock() const {
  auto elapsed = clock.elapsed_time();
  if (noisy_) {
    std::cerr << llop_names[static_cast<unsigned>(op_)] << " ";
    Clock::print(std::cerr, elapsed);
    std::cerr << std::endl;
  }
  Counter incr_sp = LLOp_Stats::scalar_prod - saved_scalar_prod;
  Counter incr_lc = LLOp_Stats::linear_comb - saved_linear_comb;
  Counter incr_sc = LLOp_Stats::sat_count - saved_sat_count;
  Counter incr_si = LLOp_Stats::sat_incl - saved_sat_incl;
  Counter incr_su = LLOp_Stats::sat_union - saved_sat_union;
  llop_stats.incr(op_, elapsed,
                  incr_sp, incr_lc, incr_sc, incr_si, incr_su);
}

} // namespace pplite
