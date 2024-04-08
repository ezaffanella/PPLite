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

#define TEST_PREDICATE_TRUE(pred)               \
  if (!pred) {                                  \
    nout << "!" #pred << endl;                  \
    ok = false;                                 \
  }

#define TEST_PREDICATE_FALSE(pred)              \
  if (pred) {                                   \
    nout << #pred << endl;                      \
    ok = false;                                 \
  }

bool
test01() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(0*A == 0);

  bool ok = true;

  TEST_PREDICATE_FALSE(ph.constrains(A));
  TEST_PREDICATE_FALSE(ph.constrains(B));

  ph.add_con(0*A == 1);

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  ph.add_gen(point());
  ph.add_gen(line(A+B));

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  ph.add_gen(line(A-B));

  TEST_PREDICATE_FALSE(ph.constrains(A));
  TEST_PREDICATE_FALSE(ph.constrains(B));

  ph.add_con(A >= 1);

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_FALSE(ph.constrains(B));

  ph.add_con(B >= 2);

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  ph.add_con(A <= B);

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3, Spec_Elem::EMPTY);

  bool ok = true;

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));
  TEST_PREDICATE_TRUE(ph.constrains(C));

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  ph.add_con(A >= B);
  ph.add_gen(point());
  ph.copy_gens();
  ph.add_con(A - B <= -1);

  bool ok = true;

  TEST_PREDICATE_TRUE(ph.constrains(C));
  TEST_PREDICATE_TRUE(ph.constrains(B));
  TEST_PREDICATE_TRUE(ph.constrains(A));

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(B >= 0);
  ph.copy_gens();
  ph.add_con(B >= 7);

  bool ok = true;

  TEST_PREDICATE_FALSE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.copy_cons();
  ph.add_gen(ray(A));
  ph.add_gen(ray(-A));

  bool ok = true;

  TEST_PREDICATE_FALSE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  return ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);

  bool ok = true;

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_FALSE(ph.constrains(B));

  ph.add_con(B >= A);

  TEST_PREDICATE_TRUE(ph.constrains(A));
  TEST_PREDICATE_TRUE(ph.constrains(B));

  return ok;
}

bool
test07() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  ph.add_con(A - B >= 1);
  ph.add_con(A <= B);

  bool ok = true;

  TEST_PREDICATE_TRUE(ph.constrains(C));
  TEST_PREDICATE_TRUE(ph.constrains(B));
  TEST_PREDICATE_TRUE(ph.constrains(A));

  return ok;
}

bool
test08() {
  Var A(0);
  Poly ph(3);
  ph.add_con(A >= 10);
  ph.add_con(A <= 0);
  auto uncon = ph.get_unconstrained();
  return uncon.empty();
}

bool
test09() {
  Poly ph(3);
  auto uncon = ph.get_unconstrained();
  Index_Set kres;
  kres.set_until(3);
  return (uncon == kres);
}

bool
test10() {
  Poly ph(20);
  ph.add_con(Var(5) == Var(10));
  auto uncon = ph.get_unconstrained();
  Index_Set kres;
  kres.set_until(20);
  kres.reset(5);
  kres.reset(10);
  return (uncon == kres);
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test07);
  DO_TEST(test08);
  DO_TEST(test09);
  DO_TEST(test10);
END_MAIN
