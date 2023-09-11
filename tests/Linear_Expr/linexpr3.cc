/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

//testing constructor with and expression and a dim
bool
test01() {
  Linear_Expr e(5);
  e.set(0,-5);
  e.set(4,12);
  Linear_Expr i(e, 10);
  return Linear_Expr::is_equal_to(e, i, 0, e.space_dim())
    && i.space_dim() == 10;
}

//start testing permute_space_dims_cycle
bool
test02() {
  Linear_Expr e(5);
  e.set(0,-15);
  e.set(1,7);
  e.set(2,22);
  e.set(3,16);
  e.set(4,12);
  Dims cycle;
  cycle.push_back(0);
  cycle.push_back(2);
  cycle.push_back(4);
  e.permute_space_dims_cycle(cycle, e.space_dim());
  return check_print(e, "12*A + 7*B - 15*C + 16*D + 22*E");
}

bool
test03() {
  Linear_Expr e(5);
  e.set(0,-15);
  e.set(1,7);
  e.set(2,22);
  e.set(3,16);
  e.set(4,12);
  Dims cycle;
  cycle.push_back(0);
  cycle.push_back(2);
  e.permute_space_dims_cycle(cycle, e.space_dim());
  return check_print(e, "22*A + 7*B - 15*C + 16*D + 12*E");
}

bool
test04() {
  Linear_Expr e(5);
  e.set(0,-15);
  e.set(1,7);
  e.set(2,22);
  e.set(3,16);
  e.set(Var(4),12);
  Dims cycle;
  cycle.push_back(0);
  e.permute_space_dims_cycle(cycle, e.space_dim());
  return check_print(e, "-15*A + 7*B + 22*C + 16*D + 12*E");
}
//end testing permute_space_dims_cycle

//start testing m_swap
bool
test05() {
  Linear_Expr e(5);
  Linear_Expr i(5);
  e.set(0,5);
  e.set(1,4);
  e.m_swap(i);
  return check_print(i, "5*A + 4*B") && e.is_zero();
}

bool
test06() {
  Linear_Expr e(6);
  Linear_Expr i(4);
  e.set(0,5);
  e.set(1,4);
  e.set(5,6);
  e.m_swap(i);
  return check_print(i, "5*A + 4*B + 6*F") && e.is_zero();
}
//end testing m_swap

//testing +=
bool
test07() {
  Linear_Expr e(6);
  Var v(6);
  e.set(0,5);
  e.set(1,4);

  e += v;
  return check_print(e, "5*A + 4*B + G");
}

//start testing -=
bool
test08() {
  Linear_Expr e(6);
  Var v(6);
  e.set(0,5);
  e.set(1,4);

  e -= v;
  return check_print(e, "5*A + 4*B - G");
}

bool
test09() {
  Linear_Expr e(6), i(7);
  e.set(0,5);
  e.set(1,4);

  i.set(0,5);
  i.set(1,4);
  e -= i;
  return check_print(e, "0");
}
//end testing -=

//start testing add_mul_assign
bool
test10() {
  Linear_Expr e(6);
  e.set(0,5);
  e.set(1,4);
  add_mul_assign(e, 5, Var(6));
  return check_print(e, "5*A + 4*B + 5*G");
}

bool
test11() {
  Linear_Expr e(6), i(7);
  e.set(0,5);
  e.set(1,4);
  e.set(2,7);
  i.set(0,1);
  i.set(1,2);
  add_mul_assign(e, 5, i);
  return check_print(e, "10*A + 14*B + 7*C");
}
//end testing add_mul_assign

//start testing sub_mul_assign
bool
test12() {
  Linear_Expr e(6);
  e.set(0,5);
  e.set(1,4);
  sub_mul_assign(e, 5, Var(6));
  return check_print(e, "5*A + 4*B - 5*G");
}

bool
test13() {
  Linear_Expr e(6), i(7);
  e.set(0,5);
  e.set(1,4);
  e.set(2,7);
  i.set(0,1);
  i.set(1,2);
  sub_mul_assign(e, 5, i);
  return check_print(e, "-6*B + 7*C");
}
//end testing sub_mul_assign

//testing linear_combine
bool
test14() {
  Linear_Expr e(3), i(3);
  e.set(0,5);
  e.set(1,4);
  e.set(2,7);
  i.set(0,3);
  i.set(1,4);
  i.set(2,3);
  Integer i1(0), i2(0);
  e.linear_combine(i,1,i1,i2);
  return check_print(e, "2*A + 4*C");
}

bool
test15() {
  Linear_Expr e(3), i(3);
  e.set(0,5);
  e.set(1,4);
  e.set(2,7);
  i.set(0,3);
  i.set(1,4);
  i.set(2,3);
  Integer i1(0), i2(0);
  e.linear_combine(i,1,i1,i2);
  return check_print(e, "2*A + 4*C");
}

//testing - with 2 variable
bool
test16() {
  Var B(1), C(2);
  Linear_Expr e = B - C;
  return check_print(e, "B - C");
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
END_MAIN
