#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Var y(1);

  MyDisj ph1(2);
  ph1.add_con(x >= 0);
  ph1.add_con(x <= 2);

  MyDisj ph2(2);
  ph2.add_con(y >= 3);
  ph2.add_con(y <= 5);

  MyPSet c_ps(2, Spec_Elem::EMPTY);
  c_ps.add_disjunct(ph1);
  c_ps.add_disjunct(ph2);

  c_ps.affine_preimage(x, x + y);

  ph1.affine_preimage(x, x + y);
  ph2.affine_preimage(x, x + y);

  auto d1 = c_ps.seq().front();
  auto d2 = c_ps.seq().back();

  bool ok1 = d1.check_inv() && d1 == ph1;

  print_cons(ph1, "*** ph1 ***");
  print_cons(d1, "*** d1 ***");

  bool ok2 = d2.check_inv() && d2 == ph2;

  print_cons(ph2, "*** ph2 ***");
  print_cons(d2, "*** d2 ***");

  return ok1 && ok2;
}

BEGIN_MAIN
  DO_TEST(test01);
END_MAIN
