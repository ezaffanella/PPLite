#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet ps(1, Topol::NNC);
  Cons cs;
  cs.clear();
  cs.push_back(x > 5);
  cs.push_back(x > 8);
  ps.add_cons(cs);

  ps.topological_closure_assign();

  bool ok = ps.check_inv();

  MyPSet known_ps(1, Topol::NNC);
  cs.clear();
  cs.push_back(x >= 5);
  cs.push_back(x >= 8);
  known_ps.add_cons(cs);

  ok = ok && ps.contains(known_ps);

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
