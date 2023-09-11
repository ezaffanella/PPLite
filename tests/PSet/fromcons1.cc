#include "pplite_test.hh"

bool
test01() {
  Cons cs;
  MyPSet ps {0, Spec_Elem::EMPTY};
  ps.add_cons(cs);
  return ps.check_inv() && ps.is_empty() && ps.space_dim() == 0;
}

bool
test02() {
  Var x(0);
  Var y(1);
  Var z(2);
  Cons cs;

  cs.push_back(x + y > 0);
  cs.push_back(x - 2*y <= 2);
  cs.push_back(z == 2);
  MyPSet ps {3, Spec_Elem::EMPTY};
  ps.add_cons(cs);
  return ps.check_inv() && ps.space_dim() == 3;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
