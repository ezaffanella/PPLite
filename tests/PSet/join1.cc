#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  MyPSet c_ps(1, Spec_Elem::EMPTY);

  MyDisj ph(1);
  ph.add_con(x >= 0);
  ph.add_con(x <= 2);

  c_ps.add_disjunct(ph);

  ph = MyDisj(1);
  ph.add_con(x >= 1);
  ph.add_con(x <= 3);

  MyPSet c_ps1(1, Spec_Elem::EMPTY);
  c_ps1.add_disjunct(ph);
  c_ps.join_assign(c_ps1);

  ph = MyDisj(1);
  ph.add_con(x >= 0);
  ph.add_con(x <= 3);

  MyPSet c_ps2(1, Spec_Elem::EMPTY);
  c_ps2.add_disjunct(ph);

  bool ok = c_ps.definitely_entails(c_ps2);
  bool ok1 = !c_ps2.definitely_entails(c_ps);

  return ok && ok1;
}

bool
test02() {
  Var x(0);

  MyPSet ps1(1, Spec_Elem::EMPTY);
  MyPSet ps2 = ps1;
  MyPSet known_res = ps1;

  MyDisj d1(1);
  d1.add_cons({ 0 <= x, x <= 2 });

  MyDisj d2(1);
  d2.add_cons({ 4 <= x, x <= 6 });

  MyDisj d3(1);
  d3.add_cons({ 8 <= x, x <= 10 });

  ps1.add_disjunct(d1);
  ps1.add_disjunct(d2);

  ps2.add_disjunct(d2);
  ps2.add_disjunct(d3);

  known_res.add_disjunct(d1);
  known_res.add_disjunct(d2);
  known_res.add_disjunct(d3);

  MyPSet res = ps1;
  res.join_assign(ps2);

  bool ok1 = ps1.definitely_entails(res);
  bool ok2 = ps2.definitely_entails(res);
  bool ok3 = not res.definitely_entails(ps1);
  bool ok4 = not res.definitely_entails(ps2);
  bool ok5 = res.equals(known_res);

  return ok1 && ok2 && ok3 && ok4 && ok5;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
