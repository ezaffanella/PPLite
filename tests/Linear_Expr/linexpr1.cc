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

bool
test01() {
  Linear_Expr e;
  return e.space_dim() == 0;
}

//start testing print function
bool
test02() {
  Linear_Expr e;
  return check_print(e, "0");
}

bool
test03() {
  Linear_Expr e(5);
  return check_print(e, "0");
}

bool
test04() {
  Linear_Expr e(5);
  e.set(0,5);
  return check_print(e, "5*A");
}

bool
test05() {
  Linear_Expr e(5);
  e.set(0,-5);
  return check_print(e, "-5*A");
}

bool
test06() {
  Linear_Expr e(5);
  e.set(0,5);
  e.set(2,-3);
  e.set(4,-7);
  return check_print(e, "5*A - 3*C - 7*E");
}

bool
test07() {
  Linear_Expr e(30);
  e.set(0,5);
  e.set(2,-3);
  e.set(4,-7);
  e.set(25,7);
  e.set(26,8);
  e.set(27,9);
  return check_print(e, "5*A - 3*C - 7*E + 7*Z + 8*A1 + 9*B1");
}
// end print function testing

//testing set_space_dim
bool
test08() {
  Linear_Expr e;
  e.set_space_dim(5);
  return e.space_dim() == 5;
}

//testing swap_space_dims
bool
test09() {
  Linear_Expr e(5);
  e.set(0,5);
  e.set(2,-3);
  e.set(4,-7);
  e.swap_space_dims(2, 4);
  return check_print(e, "5*A - 7*C - 3*E");
}

//start testing shift_space_dims
bool
test10() {
  Linear_Expr e(5);
  Var v(4);
  e.shift_space_dims(v,5);
  return e.space_dim() == 10;
}

bool
test11() {
  Linear_Expr e(5);
  Var v(2);
  e.set(0,5);
  e.set(2,-3);
  e.set(4,-7);
  e.shift_space_dims(v,5);
  return check_print(e, "5*A - 3*H - 7*J");
}
//end testing shift_space_dims

//start testing is_zero
bool
test12() {
  Linear_Expr e(5);
  return e.is_zero();
}

bool
test13() {
  Linear_Expr e(5);
  e.set(0,5);
  return !e.is_zero();
}
//end testing is_zero

//testing is_equal_to
bool
test14() {
  Linear_Expr e(6);
  Linear_Expr i(6);
  e.set(0,5);
  e.set(1,4);
  e.set(5,6);
  i.set(0,5);
  i.set(1,4);
  i.set(5,6);
  return e.is_equal_to(i);
}

bool
test15() {
  Linear_Expr e(6);
  Linear_Expr i(8);
  e.set(0,5);
  e.set(1,4);
  e.set(5,6);
  i.set(0,5);
  i.set(1,4);
  i.set(5,6);
  return e.is_equal_to(i);
}

bool
test16() {
  Linear_Expr e(6);
  Linear_Expr i(8);
  e.set(0,5);
  e.set(1,4);
  e.set(5,6);
  i.set(0,5);
  i.set(1,4);
  i.set(5,6);
  i.set(7,6);
  return !i.is_equal_to(e);
}

bool
test17() {
  Linear_Expr e(8);
  Linear_Expr i(6);
  e.set(0,5);
  e.set(1,4);
  e.set(5,6);
  i.set(0,5);
  i.set(1,4);
  i.set(5,6);
  e.set(7,6);
  return !i.is_equal_to(e);
}
//end testing is_equal_to

//testing is_equal_to with c1 and c2
bool
test18() {
  Linear_Expr e(5), i(5);
  e.set(1,21);
  e.set(2,7);
  i.set(1,42);
  i.set(2,14);
  return Linear_Expr::is_equal_to(e, i, 2, 1, 0, e.space_dim());
}

bool
test19() {
  Linear_Expr e(5), i(7);
  e.set(1,21);
  e.set(2,7);
  i.set(1,42);
  i.set(2,14);
  return e.is_equal_to(i, 2, 1);
}

bool
test20() {
  Linear_Expr e(7), i(5);
  e.set(1,21);
  e.set(2,7);
  i.set(1,42);
  i.set(2,14);
  return e.is_equal_to(i, 2, 1);
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
  DO_TEST(test19);
  DO_TEST(test20);
END_MAIN
