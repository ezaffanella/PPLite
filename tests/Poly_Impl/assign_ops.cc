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

// Recognizing efc from ns input
bool
test10() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Spec_Elem::EMPTY, Topol::NNC);
  ph1.add_gen(point(A));
  ph1.add_gen(point(2*A));
  ph1.add_gen(line(B));
  // ph1_efc = { A>=1, A<=2 }

  Poly ph2(2, Spec_Elem::EMPTY, Topol::NNC);
  ph2.add_gen(point());
  ph2.add_gen(point(4*A));
  ph2.add_gen(ray(B));
  // In gen2con efc maximal repr is ensured.
  // ph2_efc = { A>=0, A<=4, B>=0 }

  NS_Rows& ph2_ns = ph2.impl().cs.ns_rows;

  ph2.intersection_assign(ph1);
  ph2.minimize();
  // ns-loop avoids to introduce a efc different repr:
  // ph1_efc is not introduced and ph2_efc is not made redundant:
  // we still have the maximal repr.
  // ph2_efc = { A>=1, A<=2 , B>=0 }

  Poly known_result(2, Spec_Elem::EMPTY, Topol::NNC);
  known_result.add_gen(point(A));
  known_result.add_gen(point(2*A));
  known_result.add_gen(ray(B));
  known_result.minimize();
  // kr_efc = { A>=1, A<=2, B>=0 }

  NS_Rows& kr_ns = known_result.impl().cs.ns_rows;

  bool ok = (ph2 == known_result &&
             ph2_ns.size() == kr_ns.size() &&
             ph2_ns.size() == 1 &&
             ph2_ns[0].size() == kr_ns[0].size());
  return ok;
}

BEGIN_MAIN
  DO_TEST(test10);
END_MAIN
