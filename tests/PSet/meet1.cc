#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;

  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj p1 {1};
  p1.add_cons(cs);
  c_ps.add_disjunct(p1);

  MyPSet c_ps1(1, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  MyDisj p2 {1};
  p2.add_cons(cs);
  c_ps1.add_disjunct(p2);

  c_ps.meet_assign(c_ps1);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 2);
  MyPSet c_ps_expected(1, Spec_Elem::EMPTY);
  MyDisj p3 {1};
  p3.add_cons(cs);
  c_ps_expected.add_disjunct(p3);

  bool ok = c_ps.definitely_entails(c_ps_expected);
  bool ok1 = c_ps_expected.definitely_entails(c_ps);

  return ok && ok1 && c_ps.check_inv() && c_ps1.check_inv();
}

bool
test02() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;

  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj p1 {1};
  p1.add_cons(cs);
  c_ps.add_disjunct(p1);

  MyPSet c_ps1(1, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  MyDisj p2 {1};
  p2.add_cons(cs);
  c_ps1.add_disjunct(p2);

  c_ps.meet_assign(c_ps1);
  bool ok = !c_ps.is_empty();

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 2);
  MyPSet c_ps_expected(1, Spec_Elem::EMPTY);
  MyDisj p3 {1};
  p3.add_cons(cs);
  c_ps_expected.add_disjunct(p3);

  bool ok1 = c_ps.definitely_entails(c_ps_expected);
  bool ok2 = c_ps_expected.definitely_entails(c_ps);

  MyPSet c_ps2(1, Spec_Elem::EMPTY);
  cs.clear();
  cs.push_back(x == 4);
  MyDisj p4 {1};
  p4.add_cons(cs);
  c_ps2.add_disjunct(p4);

  c_ps2.meet_assign(c_ps1);
  bool ok3 = c_ps2.is_empty();

  return ok && ok1 && ok2 && ok3
    && c_ps.check_inv() && c_ps1.check_inv() && c_ps2.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
