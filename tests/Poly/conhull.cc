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
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A == 0);
  ph1.add_con(B == 1);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(A == 0);
  ph2.add_con(B == 4);

  print_cons(ph2, "*** ph2 ***");

  ph1.con_hull_assign(ph2);

  Poly kr(2);
  kr.add_con(A == 0);
  kr.add_con(B >= 1);
  kr.add_con(B <= 4);

  print_cons(ph1, "*** ph1.con_hull_assign(ph2) ***");
  print_cons(kr, "*** kr ***");

  return ph1 == kr;
}

bool
test02() {
  Var A(0);
  Var B(1);

  Poly ph1(2);
  ph1.add_con(A == 0);
  ph1.add_con(B == 1);

  print_cons(ph1, "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(A == -2);
  ph2.add_con(B == 4);

  print_cons(ph2, "*** ph2 ***");

  ph1.con_hull_assign(ph2);

  Poly kr(2);
  kr.add_con(A >= -2);
  kr.add_con(A <= 0);
  kr.add_con(B >= 1);
  kr.add_con(B <= 4);

  print_cons(ph1, "*** ph1.con_hull_assign(ph2) ***");
  print_cons(kr, "*** kr ***");

  return ph1 == kr;
}

bool
test03() {
  Var A(0);
  Var B(1);

  std::vector<Poly> ph;

  Poly ph0(2);
  ph0.add_con(A >= 0);
  ph0.add_con(B <= 2);
  ph0.add_con(B >= A + 1);
  ph.push_back(std::move(ph0));

  print_cons(ph[0], "*** ph0 ***");

  Poly ph1(2);
  ph1.add_con(A >= 0);
  ph1.add_con(A <= 1);
  ph1.add_con(B <= 2);
  ph1.add_con(B >= 1);
  ph.push_back(std::move(ph1));

  print_cons(ph[1], "*** ph1 ***");

  Poly ph2(2);
  ph2.add_con(A >= 0);
  ph2.add_con(B >= 0);
  ph2.add_con(B <= -A + 1);
  ph.push_back(std::move(ph2));

  print_cons(ph[2], "*** ph2 ***");

  Poly nary_chull(2, Spec_Elem::EMPTY);
  con_hull(nary_chull, ph.begin(), ph.end());
  nary_chull.minimize();

  print_cons(nary_chull, "*** nary_con_hull_assign ***");

  Poly kr(2);
  kr.add_con(A >= 0);
  kr.add_con(A <= 1);
  kr.add_con(B >= 0);
  kr.add_con(B <= 2);
  kr.minimize();

  print_cons(kr, "*** kr ***");

  bool ok = (nary_chull == kr);
  if (!ok)
    return false;

  Poly seq_chull(2, Spec_Elem::EMPTY);
  seq_chull.con_hull_assign(ph[0]);
  seq_chull.con_hull_assign(ph[1]);
  seq_chull.con_hull_assign(ph[2]);
  seq_chull.minimize();

  print_cons(seq_chull, "*** sequential1 con_hull_assign ***");

  ok = (seq_chull == nary_chull);
  if (!ok)
    return false;

  Poly seq_chull2(2, Spec_Elem::EMPTY);
  seq_chull2.con_hull_assign(ph[0]);
  seq_chull2.con_hull_assign(ph[2]);
  seq_chull2.con_hull_assign(ph[1]);
  seq_chull2.minimize();

  print_cons(seq_chull2, "*** sequential2 con_hull_assign ***");

  Poly kr2(2);
  kr2.add_con(A >= 0);
  kr2.add_con(B >= 0);
  kr2.add_con(B <= 2);
  kr2.add_con(B <= -A + 3);
  kr2.add_con(B >= A - 1);
  kr2.minimize();

  print_cons(kr2, "*** kr2 con_hull_assign ***");

  ok = (seq_chull2 == kr2);
  if (!ok)
    return false;

  ok = seq_chull2.contains(seq_chull) && !seq_chull.contains(seq_chull2);

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
