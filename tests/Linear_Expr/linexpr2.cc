/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

//start testing sign_normalize
bool
test01() {
  Linear_Expr e(5);
  e.set(0,-5);
  e.sign_normalize();
  return check_print(e, "5*A");
}

bool
test02() {
  Linear_Expr e(5);
  e.set(0,5);
  e.sign_normalize();
  return check_print(e, "5*A");
}

bool
test03() {
  Linear_Expr e(5);
  e.set(0,-5);
  e.set(1,3);
  e.set(2,-6);
  e.set(3,4);
  e.sign_normalize();
  return check_print(e, "5*A - 3*B + 6*C - 4*D");
}

bool
test04() {
  Linear_Expr e(5);
  e.set(0,5);
  e.set(1,-3);
  e.set(2,6);
  e.set(3,-4);
  e.sign_normalize();
  return check_print(e, "5*A - 3*B + 6*C - 4*D");
}
//end testing sign_normalize

//start testing get
bool
test05() {
  Linear_Expr e(5);
  e.set(0,5);
  return e.get(0) == 5;
}

bool
test06() {
  Linear_Expr e(5);
  e.set(0,5);
  return e.get(Var(0)) == 5;
}
//end testing get

//start testing all_zeroes
bool
test07() {
  Linear_Expr e(5);
  e.set(3,5);
  return !e.all_zeroes(2,4);
}

bool
test08() {
  Linear_Expr e(5);
  return e.all_zeroes(0,5);
}
//end testing all_zeroes

//testin num_zeroes
bool
test09() {
  Linear_Expr e(5);
  e.set(3,5);
  return e.num_zeroes(0,5) == 4;
}

//start testing gcd
bool
test10() {
  Linear_Expr e(5);
  return e.gcd(0,5) == 0;
}

bool
test11() {
  Linear_Expr e(5);
  e.set(1,21);
  e.set(2,70);
  e.set(4,-42);
  return e.gcd(0,5) == 7;
}
//end testing gcd

//testing mul_assign
bool
test12() {
  Linear_Expr e(5);
  e.set(1,21);
  e.set(2,7);
  e.set(4,-4);
  e.mul_assign(2, 0, e.space_dim());
  return check_print(e, "42*B + 14*C - 8*E");
}

//testing last_nonzero
bool
test13() {
  Linear_Expr e(5);
  e.set(1,21);
  e.set(2,7);
  return e.last_nonzero() == 2;
}

bool
test14() {
  Linear_Expr e(5);
  return e.last_nonzero() == 5;
}

//testing first_nonzero
bool
test15() {
  Linear_Expr e(5);
  e.set(1,21);
  e.set(2,7);
  return e.first_nonzero() == 1;
}

bool
test16() {
  Linear_Expr e(5);
  return e.first_nonzero() == 5;
}

//testing negate
bool
test17() {
  Linear_Expr e(5);
  e.set(1,21);
  e.set(2,7);
  e.negate(0, e.space_dim());
  return check_print(e, "-21*B - 7*C");
}

bool
test18() {
  Linear_Expr e1 = Var(0) + Var(2);
  Linear_Expr e2 = Var(0) + Var(4);
  e1.shift_space_dims(Var(1), 2);

  return e1.is_equal_to(e2);
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
  DO_TEST(test11);
  DO_TEST(test12);
  DO_TEST(test13);
  DO_TEST(test14);
  DO_TEST(test15);
  DO_TEST(test16);
  DO_TEST(test17);
  DO_TEST(test18);
END_MAIN
