#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;

  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj ph {1};
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  ph = MyDisj {1};
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  bool ok = (c_ps.size() == 2);

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
