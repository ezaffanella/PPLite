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
  // 5-dim hypercube centered in the origin
  const auto sdim = 5;
  Cons cs;
  for (auto i = 0; i < sdim; ++i) {
    cs.push_back(Var(i) >= -10);
    cs.push_back(Var(i) <= 10);
  }
  // add a skeleton strict constraint
  cs.push_back(Var(0) > -10);

  // This is the known result for skel constraints
  Poly kr(sdim, Topol::NNC);
  kr.add_cons(cs);
  kr.minimize();

  // now add non-skeleton constraints
  for (auto i = 0; i < sdim; ++i)
    for (auto j = i + 1; j < sdim; ++j)
      cs.push_back(Var(i) + Var(j) < 10 + 10);

  Poly ph(sdim, Topol::NNC);
  ph.add_cons(cs);
  ph.minimize();

  auto skel_cons = ph.skeleton_cons();
  Poly ph_sk(sdim, Topol::NNC);
  ph_sk.add_cons(skel_cons.begin(), skel_cons.end());
  ph_sk.minimize();

  print_cons(ph, "*** ph cons ***");
  print_cons(ph_sk, "*** ph skeleton cons ***");
  print_cons(kr, "*** known result cons ***");

  return (ph_sk == kr);
}

bool
test02() {
  // 3-dim open square centered in the origin
  const auto sdim = 3;
  Cons cs;
  for (auto i = 0; i < sdim; ++i) {
    cs.push_back(Var(i) > -10);
    cs.push_back(Var(i) < 10);
  }

  // This is the known result for skel constraints
  Poly kr(sdim, Topol::NNC);
  kr.add_cons(cs);
  // add a skel generator
  auto vertex = point(10*Var(0) + 10*Var(1) + 10*Var(2));
  kr.add_gen(vertex);
  kr.minimize();

  Poly ph = kr;
  // now add non-skeleton generators
  for (auto i = 0; i < sdim; ++i)
    for (auto j = i+1; j < sdim; ++j)
      ph.add_gen(point(-10*Var(i) - 10*Var(j)));
  ph.minimize();

  auto skel_gens = ph.skeleton_gens();
  Gens gs(skel_gens.begin(), skel_gens.end());
  Poly ph_sk(sdim, Topol::NNC, Spec_Elem::EMPTY);
  ph_sk.add_gens(gs);
  ph_sk.minimize();

  print_gens(ph, "*** ph gens ***");
  print_gens(ph_sk, "*** ph skeleton gens ***");
  print_gens(kr, "*** known result gens ***");

  return (ph_sk == kr);
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
