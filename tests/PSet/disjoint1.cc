#include "pplite_test.hh"

bool
test01() {
  MyPSet ps1(0, Spec_Elem::EMPTY);

  MyPSet ps2(0, Spec_Elem::EMPTY);
  bool b = ps1.is_disjoint_from(ps2);
  bool c = ps2.is_disjoint_from(ps1);

  ps1.add_disjunct(MyDisj(0));
  bool b1 = ps1.is_disjoint_from(ps2);
  bool c1 = ps2.is_disjoint_from(ps1);

  ps2.add_disjunct(MyDisj(0));
  bool b2 = !ps1.is_disjoint_from(ps2);
  bool c2 = !ps2.is_disjoint_from(ps1);

  return b && c && b1 && c1 && b2 && c2;
}

bool
test02() {
  Var x(0);
  Cons cs;
  MyPSet ps1(1, Topol::NNC, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x > 0);
  cs.push_back(x <= 1);
  MyDisj ph(1, Topol::NNC);
  ph.add_cons(cs);
  ps1.add_disjunct(ph);

  cs.clear();
  cs.push_back(x == 2);
  ph = MyDisj(1, Topol::NNC);
  ph.add_cons(cs);
  ps1.add_disjunct(ph);

  MyPSet ps2(1, Topol::NNC, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x > 2);
  cs.push_back(x <= 6);
  ph = MyDisj(1, Topol::NNC);
  ph.add_cons(cs);
  ps2.add_disjunct(ph);

  bool b = ps1.is_disjoint_from(ps2);
  bool c = ps2.is_disjoint_from(ps1);

  cs.clear();
  cs.push_back(x >= 2);
  ph = MyDisj(1, Topol::NNC);
  ph.add_cons(cs);
  ps2.add_disjunct(ph);

  bool b1 = !ps1.is_disjoint_from(ps2);
  bool c1 = !ps2.is_disjoint_from(ps1);

  return b && c && b1 && c1;
}

bool
test03() {
  Var x(0);
  Cons cs;
  MyPSet ps1(1, Spec_Elem::EMPTY);

  MyPSet ps2(1, Spec_Elem::EMPTY);
  bool b = ps1.is_disjoint_from(ps2);
  bool c = ps2.is_disjoint_from(ps1);

  ps1.add_disjunct(MyDisj(1));

  bool b1 = ps1.is_disjoint_from(ps2);
  bool c1 = ps2.is_disjoint_from(ps1);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x <= 1);
  MyDisj ph(1);
  ph.add_cons(cs);
  ps2.add_disjunct(ph);

  bool b2 = !ps1.is_disjoint_from(ps2);
  bool c2 = !ps2.is_disjoint_from(ps1);

  return b && c && b1 && c1 && b2 && c2;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
