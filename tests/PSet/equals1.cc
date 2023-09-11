#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet ps(1);
  Cons cs;
  cs.push_back(x >= 5);
  cs.push_back(x <= 3);
  ps.add_cons(cs);

  MyPSet ps1(1, Spec_Elem::EMPTY);

  bool ok = ps.geom_equals(ps1);
  bool ok1 = ps1.geom_equals(ps);

  return ok && ok1;
}

bool
test02() {
  Var x(0);
  MyPSet ps(1);
  Cons cs;
  cs.push_back(x >= 5);
  cs.push_back(x >= 8);
  ps.add_cons(cs);

  MyPSet ps1(1);
  cs.clear();
  cs.push_back(x >= 8);
  ps1.add_cons(cs);

  bool ok = ps.geom_equals(ps1);
  bool ok1 = ps1.geom_equals(ps);

  return ok && ok1;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
