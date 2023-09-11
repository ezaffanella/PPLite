#include "pplite_test.hh"

namespace {

// Powerset of C polyhedra: intersection_assign().
bool
test01() {
  auto x = Var {0};
  MyPSet c_ps { 1, Spec_Elem::EMPTY, Topol::CLOSED };
  Cons cs;

  cs.emplace_back(x >= 0);
  cs.emplace_back(x <= 2);
  MyDisj ph { 1, Spec_Elem::EMPTY, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  MyPSet c_ps1 { 1, Spec_Elem::EMPTY, Topol::CLOSED };

  cs.clear();
  cs.emplace_back(x >= 1);
  cs.emplace_back(x <= 3);
  ph = MyDisj { 1, Spec_Elem::EMPTY, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps1.add_disjunct(std::move(ph));

  c_ps.meet_assign(c_ps1);

  cs.clear();
  cs.emplace_back(x >= 1);
  cs.emplace_back(x <= 2);
  MyPSet c_ps_expected { 1, Spec_Elem::EMPTY, Topol::CLOSED };
  ph = MyDisj { 1, Spec_Elem::EMPTY, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps_expected.add_disjunct(std::move(ph));

  bool ok = c_ps.definitely_entails(c_ps_expected);
  bool ok1 = c_ps_expected.definitely_entails(c_ps);

  return ok && ok1 && c_ps.check_inv() && c_ps1.check_inv();
}

// Powerset of C polyhedra: intersection_assign().
bool
test02() {
  auto x = Var { 0 };
  MyPSet c_ps { 1, Topol::CLOSED, Spec_Elem::EMPTY };
  Cons cs;

  cs.emplace_back(x >= 0);
  cs.emplace_back(x <= 2);
  MyDisj ph { 1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  MyPSet c_ps1 { 1, Topol::CLOSED, Spec_Elem::EMPTY };

  cs.clear();
  cs.emplace_back(x >= 1);
  cs.emplace_back(x <= 3);
  ph = MyDisj { 1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps1.add_disjunct(std::move(ph));

  c_ps.intersection_assign(c_ps1);
  bool ok = !c_ps.is_empty();

  cs.clear();
  cs.emplace_back(x >= 1);
  cs.emplace_back(x <= 2);
  MyPSet c_ps_expected { 1, Topol::CLOSED, Spec_Elem::EMPTY };
  ph = MyDisj { 1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps_expected.add_disjunct(std::move(ph));

  bool ok1 = c_ps.definitely_entails(c_ps_expected);
  bool ok2 = c_ps_expected.definitely_entails(c_ps);

  MyPSet c_ps2 { 1, Topol::CLOSED, Spec_Elem::EMPTY };
  cs.clear();
  cs.emplace_back(x == 4);
  ph = MyDisj { 1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps2.add_disjunct(std::move(ph));

  c_ps2.intersection_assign(c_ps1);
  bool ok3 = c_ps2.is_empty();

  return ok && ok1 && ok2 && ok3 &&
    c_ps.check_inv() && c_ps1.check_inv() && c_ps2.check_inv();
}

bool
test03() {
  auto x = Var { 0 };
  auto y = Var { 1 };
  auto ph = MyDisj { 2, Topol::CLOSED };
  auto cs = Cons { x == 0, y == 0 };

  MyPSet ps1 { 2, Topol::CLOSED, Spec_Elem::EMPTY };

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 0, y == 0 };
  ph.add_cons(std::move(cs));
  ps1.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 1, y == 1 };
  ph.add_cons(std::move(cs));
  ps1.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 2, y == 2 };
  ph.add_cons(std::move(cs));
  ps1.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == y, x >= 3 };
  ph.add_cons(std::move(cs));
  ps1.add_disjunct(std::move(ph));

  MyPSet ps2 { 2, Topol::CLOSED, Spec_Elem::EMPTY };

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 0, y == 0 };
  ph.add_cons(std::move(cs));
  ps2.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 1, y == 1 };
  ph.add_cons(std::move(cs));
  ps2.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 2, y == 2 };
  ph.add_cons(std::move(cs));
  ps2.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == 3, y == 3 };
  ph.add_cons(std::move(cs));
  ps2.add_disjunct(std::move(ph));

  ph = MyDisj { 2, Topol::CLOSED };
  cs = Cons { x == y, x >= 4 };
  ph.add_cons(std::move(cs));
  ps2.add_disjunct(std::move(ph));

  ps1.intersection_assign(ps2);

  return ps1 == ps2;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
