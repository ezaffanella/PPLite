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
  Var A(0);
  Var B(1);

  MyDisj ph1(3);
  print_gens(ph1, "*** ph1 ***");

  MyPSet ps { ph1 };

  Vars_Set to_fold;
  to_fold.insert(A);

  ps.fold_space_dims(to_fold, B);
  ph1.fold_space_dims(to_fold, B);

  MyDisj known_result(2);

  bool ok = (ph1 == known_result)
    && ps.size() == 1 && (ps.seq().front() == known_result);

  print_gens(ph1, "*** after folding {A} into B ***");
  print_gens(ps.seq().front(), "*** after folding {A} into B ***");

  return ok;
}

bool
test02() {
  Var A(0);
  Var B(1);

  MyDisj ph1(2);
  ph1.add_con(B >= 0);
  ph1.add_con(B <= 2);

  MyDisj ph2(2);
  ph2.add_con(B >= 10);
  ph2.add_con(B <= 12);

  MyPSet ps { ph1 };
  ps.add_disjunct(ph2);
  ps.minimize();

  Vars_Set to_fold;
  to_fold.insert(B);
  ps.fold_space_dims(to_fold, A);

  return ps.is_universe();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
