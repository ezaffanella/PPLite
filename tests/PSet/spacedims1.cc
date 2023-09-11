#include "pplite_test.hh"

bool
test01() {
  Var x(0);

  MyDisj ph1(1);
  ph1.add_con(x == 1);

  MyDisj ph2(1);
  ph2.add_con(x <= 2);
  MyPSet ps(1, Spec_Elem::EMPTY);

  ps.add_disjunct(ph1);
  ps.add_disjunct(ph2);

  dim_type m = 2;

  ps.add_space_dims(m);
  bool ok = (ps.space_dim() == 3 && ps.affine_dim() == 3);

  ps.add_space_dims(m, true);
  bool ok1 = (ps.space_dim() == 5 && ps.affine_dim() == 3);

  ps.remove_higher_space_dims(4);
  bool ok2 = (ps.space_dim() == 4 && ps.affine_dim() == 3);

  MyPSet ps2(7, Spec_Elem::EMPTY);
  MyDisj ph3(7);
  ph3.add_con(x >= 1);
  ph3.add_con(x <= 0);
  ps2.add_disjunct(ph3);
  bool ok3 = (ps2.space_dim() == 7 && ps2.affine_dim() == 0);

  return ok && ok1 && ok2 && ok3;
}

bool
test02() {
  Var x(0);

  MyDisj ph1(1);
  ph1.add_con(x == 1);

  MyDisj ph2(1);
  ph2.add_con(x <= 2);
  MyPSet ps(1, Spec_Elem::EMPTY);

  ps.add_disjunct(ph1);
  ps.add_disjunct(ph2);

  dim_type m = 2;

  ps.add_space_dims(m);
  bool ok = (ps.space_dim() == 3 && ps.affine_dim() == 3);

  ps.add_space_dims(m, true);
  bool ok1 = (ps.space_dim() == 5 && ps.affine_dim() == 3);

  ps.remove_higher_space_dims(4);
  bool ok2 = (ps.space_dim() == 4 && ps.affine_dim() == 3);

  return ok && ok1 && ok2;
}

bool
test03() {
  Var x(0);
  Var y(1);
  Var z(2);
  Var w(3);

  MyDisj ph1(4);
  ph1.add_con(x == 1);
  ph1.add_con(z == 1);

  MyDisj ph2(4);
  ph2.add_con(x <= 2);
  ph2.add_con(z == 1);
  MyPSet ps(4, Spec_Elem::EMPTY);

  ps.add_disjunct(ph1);
  ps.add_disjunct(ph2);

  Vars_Set to_be_removed;
  to_be_removed.insert(y);
  to_be_removed.insert(w);

  ps.remove_space_dims(to_be_removed);
  bool ok = (ps.space_dim() == 2 && ps.affine_dim() == 1);

  return ok && ps.check_inv();
}

bool
test04() {
  Var x(0);
  Var y(1);
  Var z(2);
  Var w(3);

  MyDisj ph1(4);
  ph1.add_con(x == 1);
  ph1.add_con(z == 1);

  MyDisj ph2(4);
  ph2.add_con(x <= 2);
  ph2.add_con(z == 1);
  MyPSet ps(4, Spec_Elem::EMPTY);

  ps.add_disjunct(ph1);
  ps.add_disjunct(ph2);

  ps.expand_space_dim(y, 2);
  bool ok = (ps.space_dim() == 6 && ps.affine_dim() == 5);

  return ok && ps.check_inv();
}

bool
test05() {
  Var x(0);
  Var y(1);
  Var z(2);
  Var w(3);

  MyDisj ph1(4);
  ph1.add_con(x == 1);
  ph1.add_con(z == 1);

  MyDisj ph2(4);
  ph2.add_con(x <= 2);
  ph2.add_con(z == 1);
  MyPSet ps(4, Spec_Elem::EMPTY);

  ps.add_disjunct(ph1);
  ps.add_disjunct(ph2);

  Vars_Set to_be_folded;
  to_be_folded.insert(y);
  to_be_folded.insert(w);

  ps.fold_space_dims(to_be_folded, z);
  bool ok = (ps.space_dim() == 2 && ps.affine_dim() == 2);

  return ok && ps.check_inv();
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
//DO_TEST(test04);
  DO_TEST(test05);
END_MAIN
