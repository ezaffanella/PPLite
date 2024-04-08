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
  Poly ph(0, Spec_Elem::UNIVERSE, Topol::NNC);
  BBox bbox = ph.get_bounding_box();

  BBox expected(0, Spec_Elem::UNIVERSE);

  bool ok = (bbox == expected);

  return bbox.check_inv() && expected.check_inv() && ok;
}

bool
test02() {
  Poly ph(0, Spec_Elem::UNIVERSE);
  BBox bbox = ph.get_bounding_box();

  BBox expected(0, Spec_Elem::UNIVERSE);

  bool ok = (bbox == expected);

  return bbox.check_inv() && expected.check_inv() && ok;
}

bool
test03() {
  Poly ph(0, Spec_Elem::EMPTY, Topol::NNC);
  BBox bbox = ph.get_bounding_box();

  BBox expected(0, Spec_Elem::EMPTY);

  bool ok = (bbox == expected);

  return bbox.check_inv() && expected.check_inv() && ok;
}

bool
test04() {
  Poly ph(0, Spec_Elem::EMPTY, Topol::NNC);
  ph.add_gen(point());
  BBox bbox = ph.get_bounding_box();

  BBox expected(0, Spec_Elem::UNIVERSE);

  bool ok = (bbox == expected);

  return bbox.check_inv() && expected.check_inv() && ok;
}

bool
test05() {
  BBox bbox(0, Spec_Elem::EMPTY);
  bbox.add_gen(point());

  BBox expected(0, Spec_Elem::UNIVERSE);

  bool ok = (bbox == expected);

  return bbox.check_inv() && expected.check_inv() && ok;
}

bool
test06() {
  Var A(0);
  Var B(1);

  Poly ph(2, Spec_Elem::EMPTY, Topol::NNC);
  BBox bbox = ph.get_bounding_box();

  Gens gs = { point(), point(A), closure_point(B), closure_point(A + B) };

  bbox.add_gens(gs);

  ph.add_gens(gs);
  ph.minimize();
  BBox expected = ph.get_bounding_box();

  bool ok = (bbox == expected);

  return bbox.check_inv() && expected.check_inv() && ok;
}


BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
END_MAIN
