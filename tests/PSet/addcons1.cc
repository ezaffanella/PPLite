#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Con c = (x >= 0);
  MyPSet ps(1, Spec_Elem::EMPTY);
  ps.add_disjunct(MyDisj(1));
  ps.add_con(c);
  Con c1 = (x >= 1);
  ps.add_con(c1);
  bool ok = !ps.is_empty();

  return ok && ps.check_inv();
}

bool
test02() {
  Var x(0);
  Cons cs;
  cs.push_back(x >= 3);
  cs.push_back(x <= 4);
  MyPSet ps(1, Spec_Elem::EMPTY);
  ps.add_disjunct(MyDisj(1, Spec_Elem::UNIVERSE));
  ps.add_cons(cs);
  cs.push_back(x <= 3);
  ps.add_cons(cs);
  bool ok = !ps.is_empty();
  cs.push_back(x <= 2);
  ps.add_cons(cs);

  return ok && ps.is_empty() && ps.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
