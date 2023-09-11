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
  Poly ph;
  return ph.space_dim() == 0 && ph.affine_dim() == 0;
}

bool
test02() {
  Poly ph(0, Spec_Elem::EMPTY);
  return ph.space_dim() == 0 && ph.affine_dim() == 0;
}

bool
test03() {
  Poly ph(3);
  return ph.space_dim() == 3 && ph.affine_dim() == 3;
}

bool
test04() {
  Poly ph(3, Spec_Elem::EMPTY);
  return ph.space_dim() == 3 && ph.affine_dim() == 0;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN

