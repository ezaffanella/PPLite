/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
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

#include "pplite_test.hh"

Poly
point_poly(bool pos, const dim_type k, const dim_type dim) {
  Poly pp(dim);
  pp.add_con(Var(k) == (pos ? 1 : -1));
  for (dim_type i = 0; i < dim; ++i) {
    if (i == k)
      continue;
    pp.add_con(Var(i) == 0);
  }
  return pp;
}

bool
test01() {
  const dim_type min_dim = 12;
#ifdef NDEBUG
  const dim_type max_dim = 18;
#else
  const dim_type max_dim = 15;
#endif

  for (auto dim = min_dim; dim <= max_dim; ++dim) {
    Clock clock;
    Poly ph(dim, Spec_Elem::EMPTY);
    for (dim_type i = 0; i < dim; ++i) {
      ph.poly_hull_assign(point_poly(true, i, dim));
      ph.poly_hull_assign(point_poly(false, i, dim));
    }
    ph.minimize();
    clock.print_elapsed(nout);
    nout << " " << ph.num_min_cons()
         << " " << ph.num_min_gens()
         << endl;
  }
  return true;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
