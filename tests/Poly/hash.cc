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
hyper_cons(dim_type dim) {
  Poly ph(dim);
  for (auto i = 0; i < dim; ++i) {
    ph.add_con(Var(i) <= i);
    ph.add_con(Var(i) >= -i);
  }
  return ph;
}

bool
test01() {
  const dim_type min_dim = 4;
  const dim_type max_dim = 12;

  for (auto dim = min_dim; dim <= max_dim; ++dim) {
    Poly ph = hyper_cons(dim);
    Clock clock;
    ph.minimize();
    nout << "minimization for dim " << dim << " computed in ";
    clock.print_elapsed(nout);
    nout << " secs\n";
    clock.restart();
    auto ph_hash = ph.hash();
    nout << "hash_code " << ph_hash
         << " for dim " << dim << " computed in ";
    clock.print_elapsed(nout);
    nout << " secs\n";
  }
  return true;
}

bool
test02() {
  Var A(0); Var B(1);

  Poly x(2);
  x.add_con(B == -2);
  x.add_con(A <= 0);

  Poly y(2);
  y.add_con(A <= 0);
  y.add_con(B <= 0);

  auto x_hash = x.hash();

  // Const preserving operation
  x.is_disjoint_from(y);

  bool ok = (x_hash == x.hash());

  return ok;
}

bool
test03() {
  Var A(0); Var B(1);

  Poly x1(2);
  x1.add_con(A <= 1);
  x1.add_con(B <= 0);

  Poly x2(2);
  x2.add_con(B <= 0);
  x2.add_con(A <= 1);

  bool ok = (x1.hash() == x2.hash());

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
