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

#include "pplite-config.h"
#include "Poly_Stats.hh"

#include <iostream>
#include <iomanip>
#include <cstring>

pplite::AbsOp_Stats pplite::absop_stats;

namespace pplite {

void
AbsOp_Stats::dump_op(std::ostream& os, AbsOp op, bool dump_times) const {
  auto op_val = static_cast<unsigned>(op);
  using namespace std;
  os << "\n" << absop_names[op_val] << " ";
  const Time_Stats& ts = time_stats[op_val];
  os << ts.counter;
  if (dump_times) {
    os << " ";
    Clock::print(os, ts.time);
  }
  os << flush;
}

void
AbsOp_Stats::dump(std::ostream& os, bool dump_times) const {
  os << "=== Dumping time stats for PPLite's operations ===\n";
  for (auto i : range(absop_size))
    dump_op(os, static_cast<AbsOp>(i), dump_times);
  os << std::endl;
  dump_overall_time(os);
  os << std::endl;
}

void
AbsOp_Stats::dump_overall_time(std::ostream& os) const {
  auto overall_time = Clock::Duration::zero();
  for (const auto& ts : time_stats)
    overall_time += ts.time;
  os << "Total AbsOp time: ";
  Clock::print(os, overall_time);
  os << std::endl;
}

void
AbsOp_Stats::reset() {
  for (auto& ts : time_stats) {
    ts.counter = 0;
    ts.time = Clock::Duration::zero();
  }
}

AbsOp_Stats::~AbsOp_Stats() {
  if (noisy_dtor)
    dump(std::cerr, true);
}

} // namespace pplite
