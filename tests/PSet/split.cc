#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Var y(1);

  MyDisj ph1(2);
  ph1.add_cons( { x >= 0, x <= 3, y >= 0, y <= 3 } );

  MyDisj ph2(2);
  ph2.add_cons( { x >= 5, x <= 6, y >= 1, y <= 3 } );

  MyPSet ps(2, Spec_Elem::EMPTY);
  ps.add_disjunct(ph1);
  ps.add_disjunct(ph2);

  MyPSet ps_else = ps.integral_split(x == 2);

  MyDisj ph_eq(2);
  ph_eq.add_cons( { x == 2, 0 <= y, y <= 3 } );

  MyPSet kr_then(2, Spec_Elem::EMPTY);
  kr_then.add_disjunct(ph_eq);

  MyDisj ph_lt(2);
  ph_lt.add_cons( { 0 <= x, x <= 1, 0 <= y, y <= 3 } );
  MyDisj ph_gt(2);
  ph_gt.add_cons( { x == 3, 0 <= y, y <= 3 } );

  MyPSet kr_else(2, Spec_Elem::EMPTY);
  kr_else.add_disjunct(ph_lt);
  kr_else.add_disjunct(ph_gt);
  kr_else.add_disjunct(ph2);

  using namespace IO_Operators;
  print_pset(ps, "*** ps_then ***");
  print_pset(ps_else, "*** ps_else ***");

  return ps.equals(kr_then) && ps_else.equals(kr_else);
}


BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
