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

#include "pplite_test.hh"

bool
test01() {
  const auto dim = 8;
  Integer big = Integer("12345678901234");
  big *= big;
  Poly ph(dim);
  for (auto i = 0; i < dim; ++i) {
    ph.add_con(Var(i) >= -big);
    ph.add_con(Var(i) <= big);
  }
  ph.minimize();
  nout << "total_memory_in_bytes(ph) = "
       << total_memory_in_bytes(ph) << std::endl;
  nout << "external_memory_in_bytes(ph) = "
       << external_memory_in_bytes(ph) << std::endl;
  nout << "sizeof(ph) = " << sizeof(ph) << std::endl;

  nout << "external_memory_in_bytes(ph.impl().cs) = "
       << external_memory_in_bytes(ph.impl().cs) << std::endl;
  nout << "external_memory_in_bytes(ph.impl().gs) = "
       << external_memory_in_bytes(ph.impl().gs) << std::endl;
  nout << "external_memory_in_bytes(ph.impl().sat_c) = "
       << external_memory_in_bytes(ph.impl().sat_c) << std::endl;
  nout << "external_memory_in_bytes(ph.impl().sat_g) = "
       << external_memory_in_bytes(ph.impl().sat_g) << std::endl;
  nout << "external_memory_in_bytes(ph.impl().cs_pending) = "
       << external_memory_in_bytes(ph.impl().cs_pending) << std::endl;
  nout << "external_memory_in_bytes(ph.impl().gs_pending) = "
       << external_memory_in_bytes(ph.impl().gs_pending) << std::endl;

  const auto& c = ph.impl().cs.sk_rows[0];
  nout << "total_memory_in_bytes(Integer) = "
       << total_memory_in_bytes(c.inhomo_term()) << std::endl;
  nout << "external_memory_in_bytes(Integer) = "
       << external_memory_in_bytes(c.inhomo_term()) << std::endl;

  return true;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
