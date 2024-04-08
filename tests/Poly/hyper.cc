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

bool
test01() {
  const dim_type first_dim = 10;
#ifdef NDEBUG
  const dim_type last_dim = 16;
#else
  const dim_type last_dim = 12;
#endif

  for (dim_type dim = first_dim; dim <= last_dim; ++dim) {
    Cons cs;
    cs.reserve(2 * dim);
    for (dim_type i = 0; i < dim; ++i) {
      cs.push_back(Var(i) >= 0);
      cs.push_back(Var(i) <= 1);
    }
    nout << "Converting hypercube of dim " << dim << " | ";
    Poly ph(dim);
    ph.add_cons(cs);
    {
      Clock clock;
      ph.minimize();
      clock.print_elapsed(nout);
    }
    nout << " " << ph.num_min_cons()
         << " " << ph.num_min_gens()
         << endl;
  }
  return true;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
