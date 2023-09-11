#include "pplite_test.hh"

bool
test01() {
  MyPSet ps(0, Spec_Elem::EMPTY);
  bool b = !ps.is_universe();

  ps.add_disjunct(MyDisj(0));
  bool b1 = ps.is_universe();
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
  MyDisj ph(1, Topol::NNC, Spec_Elem::EMPTY);
  ph.add_cons(std::move(cs));
  ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x < 1);
  ph = MyDisj(1, Topol::NNC, Spec_Elem::EMPTY);
  ps.add_disjunct(ph);

  bool b = !ps.is_universe();
  return b;
}

bool
test03() {
  MyPSet ps(1, Spec_Elem::EMPTY);
  bool b = !ps.is_universe();

  ps.add_disjunct(MyDisj(1));

  bool b1 = ps.is_universe();
  return b && b1;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
END_MAIN
