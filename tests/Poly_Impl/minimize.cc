/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2023-2024 Enea Zaffanella <enea.zaffanella@unipr.it>

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

bool test01() {
  auto A = Var(0);
  auto B = Var(1);
  auto C = Var(2);
  auto D = Var(3);

  Cons cs = {
    D >= 0,
    2*A - C >= 0,
    2*A - B >= 0,
    2*A - D >= 0,
    2*B - D >= 0,
    2*C - D >= 0,
  };

  Poly ph(4, Topol::NNC);
  ph.add_cons(cs);
  ph.minimize();

  nout << "\n=== initialization ===\n";
  ph.ascii_dump(nout);

  ph.add_con(D > 3);
  nout << "\n\n=== added (pending) constraint ===\n";
  ph.ascii_dump(nout);

  ph.minimize();
  nout << "\n\n=== minimized ===\n";
  ph.ascii_dump(nout);

  return ph.check_inv();
}

bool test02() {
  auto A = Var(0);
  auto B = Var(1);
  auto C = Var(2);
  auto D = Var(3);

  Cons cs = {
    2*A - C >= 0,
    2*A - B - D >= 0,
    - B - D + 1 >= 0,
    B >= 0,
    - B + 2*C >= 0,
    2*A - 6*B + 11*C - 6*D >= 0,

    (- B - D + 1) + (B) + (- B + 2*C) + (2*A - 6*B + 11*C - 6*D) > 0
  };

  Gens gs = {
    ray(A + 2*B + 2*C + 2*D),
    ray(2*A + B + 2*C + 2*D)
  };

  Poly ph(4, Topol::NNC);
  ph.add_cons(cs);
  ph.add_gens(gs);
  ph.minimize();
  assert(ph.check_inv());

  nout << "\n\n=== added (pending) constraint ===\n";
  ph.add_gen(ray(A + B + 2*C + 2*D));
  ph.ascii_dump(nout);
  assert(ph.check_inv());

  ph.minimize();
  nout << "\n\n=== minimized ===\n";
  ph.ascii_dump(nout);

  return ph.check_inv();
}

bool test03() {
  auto A = Var(0);
  auto B = Var(1);
  auto C = Var(2);
  auto D = Var(3);

  Gens gs = {
    point(),
    point(11*A + 23*B + 23*C, 23),
    closure_point(66*A + 23*D, 23),
    ray(A),
    ray(11*A + 23*B),
  };

  Poly ph(4, Spec_Elem::EMPTY, Topol::NNC);
  ph.add_gens(gs);
  ph.minimize();
  nout << "\n\n=== minimized ===\n";
  ph.ascii_dump(nout);

  ph.add_gen(point(11*A + 23*B + 23*D, 23));
  ph.minimize();
  nout << "\n\n=== added g ===\n";
  ph.ascii_dump(nout);

  return ph.check_inv(true);
}


bool test04() {
  auto A = Var(0);
  auto B = Var(1);
  auto C = Var(2);

  Cons cs = {
    A + B >= 1,
    A <= 2,
    B >= 0,
    B <= 1,
    A + 2*B > 1,
    C >= 0
  };

  Poly ph(3, Topol::NNC);
  ph.add_cons(cs);
  ph.minimize();

  nout << "\n\n=== minimized ===\n";
  ph.ascii_dump(nout);

  ph.add_gen(point());
  ph.minimize();
  nout << "\n\n=== added g ===\n";
  ph.ascii_dump(nout);

  return ph.check_inv(true);
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
