#include "pplite_test.hh"

bool
test01() {
  MyPSet ps1(1, Spec_Elem::EMPTY);
  MyPSet ps2(1, Spec_Elem::EMPTY);
  bool b = ps1.contains(ps2);
  bool c = ps2.contains(ps1);
  bool bs = not ps1.strictly_contains(ps2);
  bool cs = not ps2.strictly_contains(ps1);

  ps1.add_disjunct(MyDisj(1));
  bool b1 = ps1.contains(ps2);
  bool c1 = not ps2.contains(ps1);
  bool bs1 = ps1.strictly_contains(ps2);
  bool cs1 = not ps2.strictly_contains(ps1);

  ps2.add_disjunct(MyDisj(1));
  bool b2 = ps1.contains(ps2);
  bool c2 = ps2.contains(ps1);
  bool bs2 = not ps1.strictly_contains(ps2);
  bool cs2 = not ps2.strictly_contains(ps1);

  bool ok = b && c && b1 && c1 && b2 && c2;
  bool oks = bs && cs && bs1 && cs1 && bs2 && cs2;

  return ok && oks;
}

bool
test02() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);
  Cons cs { x >= 0, x <= 2 };
  MyDisj ph(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  cs = { x >= 1, x <= 4 };
  ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  cs = { x >= 1, x <= 3 };
  ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps.add_disjunct(ph);

  MyPSet c_ps1(1, Spec_Elem::EMPTY);

  cs = { x >= 1, x <= 3 };
  ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps1.add_disjunct(ph);

  bool ok = c_ps.contains(c_ps1)
    && not c_ps1.contains(c_ps)
    && c_ps.strictly_contains(c_ps1)
    && not c_ps1.strictly_contains(c_ps);

  cs = { x >= 1, x <= 4 };
  ph = MyDisj(1);
  ph.add_cons(cs);
  c_ps1.add_disjunct(ph);

  bool ok1 = c_ps.contains(c_ps1)
    && not c_ps1.contains(c_ps)
    && c_ps.strictly_contains(c_ps1)
    && not c_ps1.strictly_contains(c_ps);

  return ok && ok1;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
