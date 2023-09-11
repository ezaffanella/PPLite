#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Var y(1);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;

  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  MyDisj ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  print_pset(c_ps, "=== c_ps ===");

  c_ps.concatenate_assign(c_ps);

  print_pset(c_ps, "=== c_ps.concatenate_assign(c_ps) ===");

  MyPSet kr(2, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  cs.push_back(y >= 1);
  cs.push_back(y <= 3);
  ph = MyDisj(2);
  ph.add_cons(cs);
  kr.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  cs.push_back(y >= 0);
  cs.push_back(y <= 2);
  ph = MyDisj(2);
  ph.add_cons(cs);
  kr.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  cs.push_back(y >= 0);
  cs.push_back(y <= 2);
  ph = MyDisj(2);
  ph.add_cons(cs);
  kr.add_disjunct(ph);

  cs.clear();
  cs.push_back(x >= 1);
  cs.push_back(x <= 3);
  cs.push_back(y >= 1);
  cs.push_back(y <= 3);
  ph = MyDisj(2);
  ph.add_cons(cs);
  kr.add_disjunct(ph);

  print_pset(kr, "=== kr ===");

  return c_ps.check_inv() && c_ps.equals(kr);
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
