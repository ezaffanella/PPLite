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

// The work to be done in a each thread.
void
compute_image(Poly* ph) {
  ph->affine_image(Var(0), Linear_Expr(), 5);
}

bool
test01() {
  Poly ph1 = build_polyhedron();
  Poly ph2(ph1);
  Poly ph3(ph1);
  Poly ph4(ph1);

  std::thread t1(compute_image, &ph1);
  std::thread t2(compute_image, &ph2);
  std::thread t3(compute_image, &ph3);
  std::thread t4(compute_image, &ph4);
  t1.join();
  t2.join();
  t3.join();
  t4.join();
  return (ph1 == ph2) && (ph2 == ph3) && (ph3 == ph4);
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
