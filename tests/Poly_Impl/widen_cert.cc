/* PPLite: a lightweight library for convex polyhedra derived from PPL.
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

bool
test01() {
  Cons cs1;
  cs1.push_back(Var(0) == 0);
  cs1.push_back(Var(1) == 0);
  cs1.push_back(Var(2) == 0);
  cs1.push_back(Var(3) >= 0);
  cs1.push_back(Var(3) <= 4);
  Poly ph1(Topol::NNC, 4);
  ph1.add_cons(cs1);
  ph1.minimize();

  bhrz03_widen::Cert cert = bhrz03_widen::Cert(ph1.impl());

  return cert.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
