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

using pplite::IO_Operators::operator<<;

bool
test01() {
  Var x1(0);
  Var x2(1);

  Poly ph(2);
  ph.add_con(-2*x1-x2 >= -5);
  ph.add_con(4*x1-4*x2 >= -5);
  ph.add_con(x1 >= 0);
  ph.add_con(x2 >= 0);

  print_cons(ph, "*** ph ***");

  auto ae = Affine_Expr(x1 - 2*x2);
  Rational value;
  bool included;
  Gen g(point());
  bool ok = ph.max(ae, value, &included, &g)
    && value == Rational(5, 2)
    && included
    && g.is_point()
    && g.coeff(x1) == 5 && g.coeff(x2) == 0
    && g.divisor() == 2;

  nout << (included ? "maximum" : "supremum") << " = " << value;
  nout << " @ " << g << endl;

  if (!ok)
    return false;

  ok = ph.min(ae, value, &included, &g)
    && value == Rational(-15, 4)
    && included
    && g.is_point()
    && g.coeff(x1) == 5 && g.coeff(x2) == 10
    && g.divisor() == 4;

  nout << (included ? "minimum" : "infimum") << " = " << value;
  nout << " @ " << g << endl;

  return ok;
}

bool
test02() {
  Var x1(0);
  Var x2(1);
  Var x3(2);

  Poly ph(3);
  ph.add_con(-x1-x2-x3 >= -100);
  ph.add_con(-10*x1-4*x2-5*x3 >= -600);
  ph.add_con(-x1-x2-3*x3 >= -150);
  ph.add_con(x1 >= 0);
  ph.add_con(x2 >= 0);
  ph.add_con(x3 >= 0);

  print_cons(ph, "*** ph ***");

  auto ae = Affine_Expr(-10*x1 - 6*x2 - 4*x3, 4);

  Rational value;
  bool included;
  Gen g(point());
  bool ok = ph.max(ae, value, &included, &g)
    && value == Rational(4)
    && included
    && g.is_point()
    && g.coeff(x1) == 0
    && g.coeff(x2) == 0
    && g.coeff(x3) == 0
    && g.divisor() == 1;

  nout << (included ? "maximum" : "supremum") << " = " << value;
  nout << " @ " << g << endl;

  if (!ok)
    return false;

  ok = ph.min(ae, value, &included, &g)
    && value == Rational(-2188, 3)
    && included
    && g.is_point()
    && g.coeff(x1) == 100
    && g.coeff(x2) == 200
    && g.coeff(x3) == 0
    && g.divisor() == 3;

  nout << (included ? "minimum" : "infimum") << " = " << value;
  nout << " @ " << g << endl;

  return ok;
}

bool
test03() {
  Poly ph(0);

  print_cons(ph, "*** ph ***");

  Rational value;
  bool included;
  Gen g(point());
  Affine_Expr ae;
  bool ok = ph.max(ae, value, &included, &g)
    && value == Rational::zero()
    && included
    && g.is_point()
    && g.divisor() == 1;

  nout << (included ? "maximum" : "supremum") << " = " << value;
  nout << " @ " << g << endl;

  if (!ok)
    return false;

  ok = ph.min(ae, value, &included, &g)
    && value == Rational::zero()
    && included
    && g.is_point()
    && g.divisor() == 1;

  nout << (included ? "minimum" : "infimum") << " = " << value;
  nout << " @ " << g << endl;

  return ok;
}

bool
test04() {
  Poly ph(2, Spec_Elem::EMPTY);

  print_cons(ph, "*** ph ***");

  Rational value;
  bool included = false;
  Gen g(point());
  Affine_Expr ae;
  bool ok = !ph.max(ae, value, &included, &g)
    && !ph.min(ae, value, &included, &g);

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY);
  ph.add_gen(point());
  ph.add_gen(ray(2*A + B));

  Itv itv = ph.get_bounds(A);
  nout << itv << std::endl;

  return not itv.is_empty()
    && itv.has_lb() && itv.lb == Rational::zero()
    && itv.inf_ub();
}

bool
test06() {
  Var A(0);

  Poly ph(1);
  ph.add_con(3*A >= 1);
  ph.add_con(5*A <= 2);

  Rational value;
  bool has_min, has_max, included;
  Gen g = point();
  has_min = ph.min(3*A + 5, value, &included, &g);
  nout << "has_min = " << has_min << ", "
       << "min = " << value << ", "
       << "included = " << included << ", "
       << "g = " << g << std::endl;

  if (not (has_min && value == Rational(6) &&
           included && (g == point(A, 3))))
    return false;

  has_max = ph.max(3*A + 5, value, &included, &g);
  nout << "has_max = " << has_max << ", "
       << "max = " << value << ", "
       << "included = " << included << ", "
       << "g = " << g << std::endl;

  return has_max && value == Rational(31, 5) &&
    included && (g == point(2*A, 5));
}

bool
test07() {
  Var A(0);

  Poly ph(1);
  ph.add_con(A >= 0);
  ph.add_con(A <= 1);

  Rational value;
  bool has_min, has_max, included;
  Gen g = point();
  has_min = ph.min(-A + 5, value, &included, &g);
  nout << "has_min = " << has_min << ", "
       << "min = " << value << ", "
       << "included = " << included << ", "
       << "g = " << g << std::endl;

  if (not (has_min && value == Rational(4) &&
           included && (g == point(A))))
    return false;

  has_max = ph.max(-A + 5, value, &included, &g);
  nout << "has_max = " << has_max << ", "
       << "max = " << value << ", "
       << "included = " << included << ", "
       << "g = " << g << std::endl;

  return has_max && value == Rational(5) &&
    included && (g == point());
}

bool
test08() {
  Var A(0);

  Poly ph(1);
  ph.add_con(A >= 0);
  ph.add_con(A <= 1);

  Itv itv = ph.get_bounds(-A + 5);
  nout << itv << std::endl;
  return not itv.is_empty()
    && itv.has_lb() && itv.lb == Rational(4)
    && itv.has_ub() && itv.ub == Rational(5);
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(A >= 0);
  ph.add_con(A <= 1);
  ph.add_con(B == A + 6);

  Itv itv = ph.get_bounds(A + 4);
  nout << itv << std::endl;
  return not itv.is_empty()
    && itv.has_lb() && itv.lb == Rational(4)
    && itv.has_ub() && itv.ub == Rational(5);
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
END_MAIN
