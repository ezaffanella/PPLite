/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2026 Enea Zaffanella <enea.zaffanella@unipr.it>

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
  // Testing sum of hypercubes
  for (auto dim = 4; dim != 8; ++dim) {

    // first hypercube has sides [0, 5]
    Poly x(dim);
    for (auto i : range(0, dim)) {
      x.add_con(Var(i) >= 0);
      x.add_con(Var(i) <= 5);
    }
    print_cons(x, "=== x cons ===");

    // second hypercube has sides [0, 1]
    Poly y(dim);
    for (auto i : range(0, dim)) {
      y.add_con(Var(i) >= 0);
      y.add_con(Var(i) <= 1);
    }
    print_cons(y, "=== y cons ===");

    // known results has sides [0, 6]
    Poly kr(dim);
    for (auto i : range(0, dim)) {
      kr.add_con(Var(i) >= 0);
      kr.add_con(Var(i) <= 6);
    }
    print_cons(kr, "=== known result cons ===");

    auto sum = x.minkowski_sum(y);
    print_cons(sum, "=== sum cons ===");

    if (not sum.equals(kr))
      return false;
  }
  return true;
}

bool
test02() {
  // Testing sum of (almost) hypercubes
  for (auto dim = 3; dim != 6; ++dim) {

    // first hypercube has sides [0, 5]
    Poly x(dim);
    for (auto i : range(0, dim)) {
      x.add_con(Var(i) >= 0);
      x.add_con(Var(i) <= 5);
    }
    // cut away two vertices
    Linear_Expr expr;
    for (auto i : range(0, dim))
      expr += Var(i);
    x.add_con(expr >= 0*dim + 1);
    x.add_con(expr <= 5*dim - 1);
    print_cons(x, "=== x cons ===");

    // second hypercube has sides [-3, 3]
    Poly y(dim);
    for (auto i : range(0, dim)) {
      y.add_con(Var(i) >= -3);
      y.add_con(Var(i) <= 3);
    }
    // cut away two vertices
    y.add_con(expr >= -3*dim + 1);
    y.add_con(expr <= 3*dim - 1);
    print_cons(y, "=== y cons ===");

    // known result has sides [-3, 8]
    Poly kr(dim);
    for (auto i : range(0, dim)) {
      kr.add_con(Var(i) >= -3);
      kr.add_con(Var(i) <= 8);
    }
    // cut away two vertices
    kr.add_con(expr >= -3*dim + 1 + 1);
    kr.add_con(expr <= 8*dim - 1 - 1);
    print_cons(kr, "=== kr cons ===");

    auto sum = x.minkowski_sum(y);
    print_cons(sum, "=== sum cons ===");

    if (not sum.equals(kr))
      return false;
  }
  return true;
}

bool
test03() {
  // Testing octagon + octagon
  auto dim = 2;
  Var A(0);
  Var B(1);

  Poly x(dim);
  x.add_con(A >= 0);
  x.add_con(A <= 3);
  x.add_con(B >= 0);
  x.add_con(B <= 3);
  x.add_con(A + B >= 1);
  x.add_con(A + B <= 5);
  x.add_con(A - B >= -2);
  x.add_con(A - B <= 2);
  print_cons(x, "=== x cons ===");

  Poly y = x;
  y.affine_image(A, A, 10);
  y.affine_image(B, B, 10);
  print_cons(y, "=== y cons ===");

  Poly kr(dim);
  kr.add_con(A >= 10);
  kr.add_con(A <= 16);
  kr.add_con(B >= 10);
  kr.add_con(B <= 16);
  kr.add_con(A + B >= 22);
  kr.add_con(A + B <= 30);
  kr.add_con(A - B >= -4);
  kr.add_con(A - B <= 4);
  print_cons(kr, "=== known result cons ===");

  auto sum = x.minkowski_sum(y);
  print_cons(sum, "=== sum cons ===");
  print_gens(sum, "=== sum gens ===");

  return sum.equals(kr);
}

bool
test04() {
  // Testing half-stripe + halp stripe
  auto dim = 2;
  Var A(0);
  Var B(1);

  Poly x(dim);
  x.add_con(A >= 0);
  x.add_con(A <= 3);
  x.add_con(B >= 0);
  print_cons(x, "=== x cons ===");

  Poly y(dim);
  y.add_con(B >= 0);
  y.add_con(B <= 3);
  y.add_con(A >= 0);
  print_cons(y, "=== y cons ===");

  Poly kr(dim);
  kr.add_con(A >= 0);
  kr.add_con(B >= 0);
  print_cons(kr, "=== known result cons ===");

  auto sum = x.minkowski_sum(y);
  print_cons(sum, "=== sum cons ===");
  print_gens(sum, "=== sum gens ===");

  return sum.equals(kr);
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN

