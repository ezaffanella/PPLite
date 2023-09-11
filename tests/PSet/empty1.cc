#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;

  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj ph(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  MyPSet c_ps1;
  c_ps1 = c_ps;

  bool ok = !c_ps.is_empty();
  return ok;
}

bool
test02() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;

  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj ph(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  MyPSet c_ps1(1, Spec_Elem::EMPTY);
  swap(c_ps, c_ps1);

  bool ok = (c_ps.is_empty() && !c_ps1.is_empty());
  return ok;
}

bool
test03() {
  MyPSet ps(0, Spec_Elem::EMPTY);
  bool b = ps.is_empty();

  ps.add_disjunct(MyDisj(0));
  bool b1 = !ps.is_empty();
  return b && b1;
}

bool
test04() {
  Var x(0);
  Cons cs;
  MyPSet ps(1, Topol::NNC, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x > 0);
  cs.push_back(x <= 1);
  MyDisj ph(1, Topol::NNC);
  ph.add_cons(cs);
  ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x < 1);
  ph = MyDisj(1, Topol::NNC);
  ph.add_cons(cs);
  ps.add_disjunct(ph);

  bool b = !ps.is_empty();
  return b;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
