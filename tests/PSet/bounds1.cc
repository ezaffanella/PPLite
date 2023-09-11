#include "pplite_test.hh"

bool
test01() {
  MyPSet ps(0, Spec_Elem::EMPTY);
  Linear_Expr zero;
  bool ok1 = ps.is_bounded_expr(false, zero);
  bool ok2 = ps.is_bounded_expr(true, zero);

  ps.add_disjunct(MyDisj(0));
  bool ok3 = ps.is_bounded_expr(false, zero);
  bool ok4 = ps.is_bounded_expr(true, zero);

  return ok1 && ok2 && ok3 && ok4;
}

bool
test02() {
  Var x(0);
  Cons cs = { x > 0, x <= 1 };
  MyDisj ph(1, Spec_Elem::UNIVERSE, Topol::NNC);
  ph.add_cons(cs);

  MyPSet ps(1, Spec_Elem::EMPTY, Topol::NNC);
  ps.add_disjunct(ph);

  ph = MyDisj(1, Spec_Elem::UNIVERSE, Topol::NNC);
  ph.add_con(x > 1);
  ps.add_disjunct(ph);

  bool ok1 = not ps.is_bounded_expr(false, x);
  bool ok2 = ps.is_bounded_expr(true, x);

  return ok1 && ok2;
}

bool
test03() {
  Var x(0);
  MyPSet ps(1, Spec_Elem::EMPTY);

  bool ok1 = ps.is_bounded_expr(false, x);
  bool ok2 = ps.is_bounded_expr(true, x);

  if (!ok1 || !ok2)
    return false;

  ps.add_disjunct(MyDisj(1));

  ok1 = !ps.is_bounded_expr(false, x);
  ok2 = !ps.is_bounded_expr(true, x);

  return ok1 && ok2;
}

BEGIN_MAIN
//  DO_TEST(test01);
  DO_TEST(test02);
//  DO_TEST(test03);
END_MAIN
