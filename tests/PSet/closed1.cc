#include "pplite_test.hh"

bool
test01() {
  MyPSet ps(0, Spec_Elem::EMPTY);
  bool b = ps.is_topologically_closed();

  ps.add_disjunct(MyDisj(0));
  bool b1 = ps.is_topologically_closed();
  return b && b1;
}

bool
test02() {
  Var x(0);
  Cons cs;
  MyPSet ps(1, Topol::NNC, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x > 0);
  cs.push_back(x <= 1);
  MyDisj ph(1, Topol::NNC, Spec_Elem::UNIVERSE);
  ph.add_cons(cs);
  ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  ph = MyDisj(1, Topol::NNC, Spec_Elem::UNIVERSE);
  ph.add_cons(cs);
  ps.add_disjunct(ph);

  bool b = ps.is_topologically_closed();
  return b;
}

bool
test03() {
  Var x(0);
  Cons cs;
  MyPSet ps(1, Topol::NNC, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x <= 1);
  MyDisj ph = MyDisj(1, Topol::NNC, Spec_Elem::UNIVERSE);
  ph.add_cons(cs);
  ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x > 0);
  cs.push_back(x < 2);
  ph = MyDisj(1, Topol::NNC, Spec_Elem::UNIVERSE);
  ph.add_cons(cs);
  ps.add_disjunct(ph);

  bool b = !ps.is_topologically_closed();
  return b;
}

bool
test04() {
  MyPSet ps(1, Spec_Elem::EMPTY);
  bool b = ps.is_topologically_closed();

  ps.add_disjunct(MyDisj(1));

  bool b1 = ps.is_topologically_closed();
  return b && b1;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
END_MAIN
