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
#include <vector>
#include <algorithm>

namespace {

int SUCC;
int TOT;

struct BPoly {
  Poly ph;
  mutable BBox bbox;
  mutable bool valid;
  BPoly(const Poly& _ph)
    : ph(_ph), bbox(ph.space_dim()), valid(false) {}

  void maybe_compute_bbox() const {
    if (valid) return;
    bbox = ph.get_bounding_box();
    valid = true;
  }

  void intersection_assign(const BPoly& y) {
    maybe_compute_bbox();
    y.maybe_compute_bbox();
    if (bbox.is_disjoint_from(y.bbox)) {
      ph.set_empty();
      bbox.set_empty();
      ++SUCC;
      return;
    }
    ph.intersection_assign(y.ph);
    valid = false;
  }
};

bool operator==(const Poly& ph, const BPoly& bph) {
  return ph == bph.ph;
}

using P_Vec = std::vector<Poly>;
using BP_Vec = std::vector<BPoly>;

template <typename Set>
void intersection_assign(Set& xps, const Set& yps) {
  for (auto& x : xps)
    for (const auto& y : yps) {
      ++TOT;
      x.intersection_assign(y);
    }
}

Poly
build_hyper(dim_type dim, dim_type lb, dim_type ub) {
  Cons cs;
  for (dim_type i = dim; i-- > 0; ) {
    cs.push_back(Var(i) <= ub);
    cs.push_back(Var(i) >= lb);
  }
  Poly ph(dim, Topol::NNC);
  ph.add_cons(cs);
  ph.minimize();
  return ph;
}

P_Vec
build_powerset(dim_type dim, dim_type card,
               dim_type size, dim_type step) {
  P_Vec ps;
  ps.reserve(card);
  dim_type lb = 0;
  dim_type ub = size;
  for (dim_type i = 0; i < card; ++i) {
    ps.push_back(build_hyper(dim, lb, ub));
    lb += step;
    ub += step;
  }
  return ps;
}

bool
test01() {
  SUCC = 0;
  dim_type n1 = 50, size1 = 10, step1 = 10;
  dim_type n2 = 100, size2 = 5, step2 = 5;
  dim_type dim = 2;

  nout << "HyperP_Vec1: card: " << n1 << "  dim: " << dim
       << "  size: " << size1 << "  step: " << step1 << endl;
  nout << "HyperP_Vec2: card: " << n2 << "  dim: " << dim
       << "  size: " << size2 << "  step: " << step2 << endl;
  P_Vec ps1 = build_powerset(dim, n1, size1, step1);
  P_Vec ps2 = build_powerset(dim, n2, size2, step2);
  BP_Vec bps1(ps1.begin(), ps1.end());
  BP_Vec bps2(ps2.begin(), ps2.end());

  {
    nout << "Boxed intersection:\t";
    Clock clock;
    intersection_assign(bps1, bps2);
    clock.print_elapsed(nout);
    nout << endl;
  }

  TOT = 0;
  {
    nout << "Standard intersection:\t";
    Clock clock;
    intersection_assign(ps1, ps2);
    clock.print_elapsed(nout);
    nout << endl;
  }

  nout << "*** Test avoided " << SUCC << "/" << TOT << " times." << endl;

  bool ok = std::equal(ps1.begin(), ps1.end(), bps1.begin());

  return ok;
}

bool
test02() {
  SUCC = 0;
  TOT = 0;
  dim_type n1 = 20, size1 = 10, step1 = 0;
  dim_type n2 = 20, size2 = 10, step2 = 0;
  dim_type dim = 5;

  // Measuring the overhead.
  // No empty intersections.
  nout << "HyperP_Vec1: card: " << n1 << "  dim: " << dim
       << "  size: " << size1 << "  step: " << step1 << endl;
  nout << "HyperP_Vec2: card: " << n2 << "  dim: " << dim
       << "  size: " << size2 << "  step: " << step2 << endl;
  P_Vec ps1 = build_powerset(dim, n1, size1, step1);
  P_Vec ps2 = build_powerset(dim, n2, size2, step2);
  BP_Vec bps1(ps1.begin(), ps1.end());
  BP_Vec bps2(ps2.begin(), ps2.end());

  {
    nout << "Boxed intersection:\t";
    Clock clock;
    intersection_assign(bps1, bps2);
    clock.print_elapsed(nout);
    nout << endl;
  }

  TOT = 0;
  {
    nout << "Standard intersection:\t";
    Clock clock;
    intersection_assign(ps1, ps2);
    clock.print_elapsed(nout);
    nout << endl;
  }

  bool ok = std::equal(ps1.begin(), ps1.end(), bps1.begin());

  nout << "*** Test avoided " << SUCC << "/" << TOT << " times." << endl;

  return ok;
}

bool
test03() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gen(point());
  ph1.add_gen(closure_point(2*A));
  ph1.add_gen(closure_point(2*B));
  ph1.minimize();

  Poly ph2(2, Topol::NNC, Spec_Elem::EMPTY);;
  ph2.add_gen(point(A + 2*B));
  ph2.add_gen(point(2*A + B));
  ph2.add_gen(point(A + B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();
  bool ok = !box1.is_disjoint_from(box2);
  if (!ok) return false;

  ok = boxed_is_disjoint_from(ph1, ph2, box1, box2);
  return ok;
}

bool
test04() {
  Var A(0);
  Var B(1);

  Poly ph1(2, Topol::NNC, Spec_Elem::EMPTY);
  ph1.add_gen(point(4*A + 2*B));
  ph1.add_gen(point(2*A + 4*B));
  ph1.add_gen(ray(A + B));
  ph1.minimize();

  Poly ph2(2, Topol::NNC, Spec_Elem::EMPTY);;
  ph2.add_gen(point());
  ph2.add_gen(point(A));
  ph2.add_gen(point(B));
  ph2.minimize();

  BBox box1 = ph1.get_bounding_box();
  BBox box2 = ph2.get_bounding_box();

  bool ok = boxed_is_disjoint_from(ph1, ph2, box1, box2);
  return ok;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
