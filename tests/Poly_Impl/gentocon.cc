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
  Var A(0);
  Var B(1);
  Var C(2);
  Var D(3);
  Var E(4);
  Var F(5);

  Poly ph(6, Topol::NNC);
  ph.add_con(A == 1);
  ph.add_con(B == 0);
  ph.add_con(E -10*F == 0);
  ph.add_con(E < 0);
  ph.minimize();

  ph.add_gen(line(D));
  ph.add_gen(line(C));
  ph.add_gen(closure_point(A + B));
  ph.add_gen(ray(-10*E - F));
  Index_Set ns;
  ns.set(0);
  ns.set(1);
  ph.impl().gs_pending.ns_rows.push_back(ns);

  ph.minimize();

  return (ph.check_inv());
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
