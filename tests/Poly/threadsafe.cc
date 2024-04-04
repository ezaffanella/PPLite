/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
   Copyright (C) 2018-2023 Enea Zaffanella <enea.zaffanella@unipr.it>

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
#include <thread>

Poly
build_polyhedron() {
  const dim_type dim = 10;
  Poly ph(dim);
  for (dim_type i = 0; i < dim; ++i) {
    ph.add_con(Var(i) >= 0);
    ph.add_con(Var(i) <= 1);
  }
  return ph;
}

// The work to be done in each thread.
void
compute_image(Poly* ph, int tid) {
  std::stringstream ss;
  ss << "start of thread " << tid << "\n";
  nout << ss.str();
  ph->affine_image(Var(0), Linear_Expr(), 5);
  ph->minimize();
  ss = std::stringstream();
  ss << "end of thread " << tid << "\n";
  nout << ss.str();
}

bool
test01() {
  const int num_threads = 8;
  std::vector<std::thread> tds { num_threads };
  std::vector<Poly> phs { num_threads, build_polyhedron() };

  for (int i = 0; i < num_threads; ++i)
    tds[i] = std::thread(compute_image, &phs[i], i);

  for (int i = 0; i < num_threads; ++i)
    tds[i].join();

  return phs.end() == std::adjacent_find(phs.begin(), phs.end(),
                                         std::not_equal_to<Poly>());
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
