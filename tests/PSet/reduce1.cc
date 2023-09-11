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
  cs.push_back(x >= 0);
  cs.push_back(x <= 3);
  ph = MyDisj{1};
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));
  c_ps.omega_reduce();

  bool ok = (c_ps.size() == 1);

  return ok;
}

bool
test02() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);

  MyDisj ph {1};
  ph.add_con(x >= 0);
  ph.add_con(x <= 10);
  c_ps.add_disjunct(std::move(ph));

  ph = MyDisj{1};
  ph.add_con(x >= 2);
  ph.add_con(x <= 4);
  c_ps.add_disjunct(std::move(ph));

  ph = MyDisj{1};
  ph.add_con(x >= 5);
  ph.add_con(x <= 7);
  c_ps.add_disjunct(std::move(ph));

  c_ps.omega_reduce();
  bool ok = (c_ps.size() == 1);
  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
