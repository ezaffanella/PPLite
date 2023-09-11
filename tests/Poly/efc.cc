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

bool
test01() {
  Var A(0);
  Var B(1);
  Poly ph1(Topol::NNC, 2);
  ph1.add_con(A >= 0);
  ph1.add_con(B >= 0);
  ph1.add_con(A + B > 0);
  ph1.minimize();

  Poly ph2(ph1);
  ph2.add_con(A + B <= 1);
  ph2.add_gen(ray(A));
  ph2.add_gen(ray(B));
  ph2.minimize();

  return ph1 == ph2 && ph1.check_inv() && ph2.check_inv();
}

bool
test02() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph(Topol::NNC, 4, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(A));
  ph.add_gen(point(B));
  ph.add_gen(point(C));
  ph.minimize();
  ph.add_gen(ray(A+B+C));
  ph.minimize();
  ph.ascii_dump(nout);

  ph.add_gen(ray(A));
  ph.minimize();
  ph.ascii_dump(nout);

  ph.add_gen(ray(B));
  ph.minimize();
  ph.ascii_dump(nout);

  return ph.check_inv();
}

bool
test03() {
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);

  Poly ph(Topol::NNC, 4, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(A));
  ph.add_gen(line(B));
  ph.add_gen(line(C));
  ph.add_gen(line(D));
  ph.minimize();
  ph.ascii_dump(nout);
  ph.add_con(B > 0);
  ph.minimize();
  ph.ascii_dump(nout);

  return ph.check_inv();
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph(Topol::NNC, 2, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(point(A));
  ph.add_gen(point(B));
  ph.add_gen(point(A + B));
  ph.minimize();
  nout << "Num square cons: " << ph.num_min_cons() << std::endl;
  print_cons(ph, "Square cons");

  ph.add_gen(ray(B));
  ph.minimize();
  nout << "Num half-stripe cons: " << ph.num_min_cons() << std::endl;
  print_cons(ph, "Half-stripe cons");

  ph.add_gen(ray(A));
  ph.minimize();
  nout << "Num quarter cons: " << ph.num_min_cons() << std::endl;
  print_cons(ph, "Quarter cons");

  return ph.check_inv();
}

bool
test05() {
  Var A(0);
  Var B(1);
  Var C(2);

  Cons cs;
  cs.push_back(A >= 0);
  cs.push_back(A <= 2);
  cs.push_back(B >= 0);
  cs.push_back(B <= 2);
  cs.push_back(C >= 0);
  cs.push_back(A + B > 0);
  Poly ph(Topol::NNC, 3);
  ph.add_cons(cs);
  ph.minimize();

  return ph.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
END_MAIN

