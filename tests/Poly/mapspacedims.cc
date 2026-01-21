/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2001-2010 Roberto <bagnara@cs.unipr.it>
   Copyright (C) 2010-2018 BUGSENG srl (http://bugseng.com)
   Copyright (C) 2018-2026 Enea Zaffanella <enea.zaffanella@unipr.it>

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

namespace {

struct PFunc : public Dims {
  using Dims::Dims;
  PFunc(dim_type dim) : Dims(dim, not_a_dim()) {}
};

void
print_function(const Dims& pfunc,
               const std::string& intro,
               std::ostream& os = nout) {
  using namespace pplite::IO_Operators;
  if (!intro.empty())
    os << intro << "\n";
  for (dim_type i = 0, i_end = num_rows(pfunc); i < i_end; ++i) {
    if (pfunc[i] == not_a_dim())
      os << Var(i) << " --> not mapped\n";
    else
      os << Var(i) << " --> " << Var(pfunc[i]) << "\n";
  }
}

} // namespace

bool
test01() {
  PFunc function(3);

  Poly ph1(3);

  print_function(function, "*** function ***");
  print_cons(ph1, "*** ph1 ***");

  ph1.map_space_dims(function);

  Poly known_result;

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.map_space_dims(function) ***");

  return ok;
}

bool
test02() {
  PFunc function(3);

  Poly ph1(3, Spec_Elem::EMPTY);

  print_function(function, "*** function ***");
  print_cons(ph1, "*** ph1 ***");

  ph1.map_space_dims(function);

  Poly known_result(0, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.map_space_dims(function) ***");

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);
  Var C(2);

  PFunc function(3);
  function[0] = 2;
  function[2] = 0;
  function[1] = 1;

  Gens gs;
  gs.push_back(point(2*C));
  gs.push_back(line(A + B));
  gs.push_back(ray(A + C));

  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gens(gs);

  print_function(function, "*** function ***");
  print_gens(ph1, "*** ph1 ***");

  ph1.map_space_dims(function);

  Gens known_gs;
  known_gs.push_back(point(2*A));
  known_gs.push_back(line(C + B));
  known_gs.push_back(ray(C + A));
  Poly known_result(3, Spec_Elem::EMPTY);
  known_result.add_gens(known_gs);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after ph1.map_space_dims(function) ***");

  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);
  Var C(2);

  PFunc function(3);
  function[0] = 1;
  function[2] = 0;

  Gens gs;
  gs.push_back(point());
  gs.push_back(ray(A + B));
  gs.push_back(ray(A - C));

  Poly ph1(3, Spec_Elem::EMPTY);
  ph1.add_gens(gs);

  print_function(function, "*** function ***");
  print_gens(ph1, "*** ph1 ***");

  ph1.map_space_dims(function);

  Gens known_gs;
  known_gs.push_back(point());
  known_gs.push_back(ray(B));
  known_gs.push_back(ray(B - A));
  Poly known_result(2, Spec_Elem::EMPTY);
  known_result.add_gens(known_gs);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after ph1.map_space_dims(function) ***");

  return ok;
}

bool
test05() {
  Var A(0);
  Var B(1);

  PFunc function(2);
  function[0] = 0;
  function[1] = 1;

  Gens gs;
  gs.push_back(point());
  gs.push_back(point(A));
  gs.push_back(point(B));
  gs.push_back(point(A + B));

  Poly ph1(2, Spec_Elem::EMPTY);
  ph1.add_gens(gs);

  Poly known_result(ph1);

  print_function(function, "*** function ***");
  print_gens(ph1, "*** ph1 ***");

  ph1.map_space_dims(function);

  bool ok = (ph1 == known_result);

  print_gens(ph1, "*** after ph1.map_space_dims(function) ***");

  return ok;
}

bool
test06() {
  PFunc function(3);
  function[0] = 1;
  function[1] = 0;

  Poly ph1(3, Spec_Elem::EMPTY);

  print_function(function, "*** function ***");
  print_cons(ph1, "*** ph1 ***");

  ph1.map_space_dims(function);

  Poly known_result(2, Spec_Elem::EMPTY);

  bool ok = (ph1 == known_result);

  print_cons(ph1, "*** after ph1.map_space_dims(function) ***");

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);
  Var z(2);

  PFunc rotate_right(3);
  rotate_right[0] = 1;
  rotate_right[1] = 2;
  rotate_right[2] = 0;

  PFunc rotate_left(3);
  rotate_left[0] = 2;
  rotate_left[1] = 0;
  rotate_left[2] = 1;

  Poly ph(3);
  ph.add_con(-4*x - 2*y + z >= -8);
  ph.add_con(-4*x + 2*y + z >= 4);
  ph.add_con(-2*x - y + 2*z >= -1);
  ph.add_con(-2*x + y + 2*z >= 5);
  ph.add_con(-x - y - 2*z >= -13);
  ph.add_con(-x - z >= -5);
  ph.add_con(-x >= -1);
  ph.add_con(-x + y - 2*z >= -7);
  ph.add_con(-y >= -4);
  ph.add_con(y >= 2);
  ph.add_con(x >= 0);

  print_cons(ph, "*** ph ***");
  print_function(rotate_right, "*** rotate_right ***");
  print_function(rotate_left, "*** rotate_left ***");

  Poly rs[4];
  rs[0] = ph;

  print_cons(rs[0], "*** rs[0] ***");

  for (int i = 1; i <= 3; ++i) {
    rs[i] = rs[i-1];
    rs[i].map_space_dims(rotate_right);
    auto msg_i = std::string("*** rs[") + std::to_string(i) + "] ***";
    print_cons(rs[i], msg_i);

  }

  Poly ls[4];
  ls[3] = ph;

  print_cons(ls[3], "*** ls[3] ***");

  for (int i = 2; i >= 0; --i) {
    ls[i] = ls[i+1];
    // Force generators to be up-to-date, for a change.
    (void) ls[i].copy_gens();
    ls[i].map_space_dims(rotate_left);

    auto msg_i = std::string("*** ls[") + std::to_string(i) + "] ***";
    print_cons(ls[i], msg_i);

  }

  for (int i = 0; i <= 3; ++i)
    if (rs[i] != ls[i]) {
      nout << "rs[" << i << "] != ls[" << i << "]" << endl;
      return false;
    }

  return true;
}

bool
test08() {
  Var A(0);
  Var B(1);
  Var C(2);

  Poly ph(3);
  ph.add_con(A >= 2);
  ph.add_con(B >= 1);
  ph.add_con(C >= 0);

  PFunc rotate_right(3);
  rotate_right[0] = 1;
  rotate_right[1] = 2;
  rotate_right[2] = 0;

  print_cons(ph, "*** ph ***");
  print_function(rotate_right, "*** rotate_right ***");

  ph.map_space_dims(rotate_right);

  Poly known_result(3);
  known_result.add_con(A >= 0);
  known_result.add_con(B >= 2);
  known_result.add_con(C >= 1);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.map_space_dims(rotate_right) ***");

  return ok;
}

bool
test09() {
  Var A(0);
  Var B(1);

  Poly ph(2);
  ph.add_con(-A + B == 0);

  PFunc rotate_right(2);
  rotate_right[0] = 1;
  rotate_right[1] = 0;

  print_cons(ph, "*** ph ***");
  print_function(rotate_right, "*** rotate_right ***");

  ph.map_space_dims(rotate_right);

  Poly known_result(2);
  known_result.add_con(A == B);

  bool ok = (ph == known_result);

  print_cons(ph, "*** after ph.map_space_dims(rotate_right) ***");

  return ok;
}

bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph(5);
  ph.add_con(A >= 1);
  ph.add_con(A <= 0);

  PFunc f(5);
  f[0] = 4;
  f[1] = 3;
  f[2] = 2;
  f[3] = 1;
  f[4] = 0;

  ph.map_space_dims(f);
  ph.minimize();
  ph.ascii_dump(nout);

  return ph.check_inv();
}

bool
test11() {
  Var A(0);
  Var B(1);

  Poly ph(5, Spec_Elem::EMPTY);

  PFunc f(5);
  f[0] = 2;
  f[1] = 1;
  f[2] = 0;
  f[3] = 3;
  f[4] = 4;

  ph.map_space_dims(f);
  ph.ascii_dump(nout);

  return ph.check_inv();
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
END_MAIN
