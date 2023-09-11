#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Cons cs;
  MyPSet c_ps(1, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj ph(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  c_ps.add_con(x == 1);

  MyPSet nnc_ps = c_ps;
  nnc_ps.set_topology(Topol::NNC);

  for (const auto& c_phi : c_ps)
    print_cons(c_phi, "*** c_phi ***");
  for (const auto& nnc_phi : nnc_ps)
    print_cons(nnc_phi, "*** nnc_phi ***");

  return c_ps.check_inv() && nnc_ps.check_inv();
}

bool
test02() {
  Var x(0);
  Cons cs;
  MyPSet nnc_ps(1, Topol::NNC, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x > 0);
  cs.push_back(x <= 1);
  MyDisj ph = MyDisj(1, Topol::NNC);
  ph.add_cons(cs);
  nnc_ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x < 1);
  ph = MyDisj(1, Topol::NNC);
  ph.add_cons(cs);
  nnc_ps.add_disjunct(ph);

  MyPSet c_ps = nnc_ps;
  c_ps.topological_closure_assign();
  c_ps.set_topology(Topol::NNC);

  for (const auto& c_phi : c_ps)
    print_cons(c_phi, "*** c_phi ***");
  for (const auto& nnc_phi : nnc_ps)
    print_cons(nnc_phi, "*** nnc_phi ***");

  return nnc_ps.check_inv() && c_ps.check_inv();
}

bool
test03() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;
  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj ph1(1);
  ph1.add_cons(cs);
  c_ps.add_disjunct(ph1);
  c_ps.seq().erase(c_ps.begin());

  bool ok = c_ps.is_empty();

  MyDisj ph2(1);
  ph2.add_cons(cs);
  c_ps.add_disjunct(ph2);

  cs.push_back(x >= 0);
  cs.push_back(x <= 3);
  MyDisj ph3(1);
  ph3.add_cons(cs);
  c_ps.add_disjunct(ph3);
  c_ps.seq().erase(c_ps.begin(), c_ps.end());

  bool ok1 = c_ps.is_empty();

  return ok && ok1;
}

bool
test04() {
  Var x(0);
  Var y(1);
  MyPSet c_ps(2, Spec_Elem::EMPTY);
  Cons cs1;
  cs1.push_back(x >= 0);
  cs1.push_back(x <= 2);
  MyDisj ph1(2);
  ph1.add_cons(cs1);
  c_ps.add_disjunct(ph1);

  Cons cs2;
  cs2.push_back(y >= 3);
  cs2.push_back(y <= 5);
  MyDisj ph2(2);
  ph2.add_cons(cs2);
  c_ps.add_disjunct(ph2);

  auto i = c_ps.begin();
  auto phi1 = *i;

  bool ok1 = phi1.check_inv() && phi1 == ph1;

  print_cons(ph1, "*** ph1 ***");
  print_cons(phi1, "*** phi1 ***");

  ++i;
  auto phi2 = *i;
  bool ok2 = phi2.check_inv() && phi2 == ph2;

  print_cons(ph2, "*** ph2 ***");
  print_cons(phi2, "*** phi2 ***");

  phi1.poly_hull_assign(phi2);
  print_cons(phi1, "*** phi1 ***");

  return ok1 && ok2;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
