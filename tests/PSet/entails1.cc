#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs;
  cs.push_back(x >= 0);
  MyDisj ph(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  MyPSet c_ps1(1, Spec_Elem::EMPTY);
  Cons cs1;
  cs1.push_back(x >= 0);
  cs1.push_back(x <= 2);
  ph = MyDisj(1);
  ph.add_cons(cs1);
  c_ps1.add_disjunct(ph);

  bool ok = c_ps1.definitely_entails(c_ps);

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
