#include "pplite_test.hh"

bool
test01() {
  Var x(0);
  Var y(1);
  Var z(1);
  MyPSet c_ps(3, Spec_Elem::EMPTY);
  Cons cs;
  cs.push_back(x >= 0);
  cs.push_back(x <= 2);
  cs.push_back(z <= 2);
  cs.push_back(z >= 2);
  MyDisj ph(3);
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  Cons cs1;
  cs1.push_back(y >= 3);
  cs1.push_back(y <= 5);
  cs1.push_back(x == 6);
  MyDisj ph1(3);
  ph1.add_cons(std::move(cs1));
  c_ps.add_disjunct(std::move(ph1));

  dim_type d = c_ps.affine_dim();
  bool ok = (d == 3);

  auto i = c_ps.begin();
  print_cons(*i, "*** phi ***");
  ++i;
  print_cons(*i, "*** phi1 ***");

  c_ps.add_con(z == 2);

  dim_type d1 = c_ps.affine_dim();
  bool ok1 = (d1 == 2);

  auto j = c_ps.begin();
  print_cons(*j, "*** phj ***");
  ++j;
  print_cons(*j, "*** phj1 ***");

  return ok && ok1;
}

bool
test02() {
  Var x(0);
  Var y(1);
  Var z(1);
  MyPSet c_ps(3, Spec_Elem::EMPTY, Topol::NNC);
  Cons cs;
  cs.push_back(x > 0);
  cs.push_back(x <= 2);
  cs.push_back(z <= 2);
  cs.push_back(z >= 2);
  MyDisj ph(3, Topol::NNC);
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  Cons cs1;
  cs1.push_back(y >= 3);
  cs1.push_back(y <= 5);
  cs1.push_back(x == 6);
  MyDisj ph1(3, Topol::NNC);
  ph1.add_cons(std::move(cs1));
  c_ps.add_disjunct(std::move(ph1));

  dim_type d = c_ps.affine_dim();
  bool ok = (d == 3);

  auto i = c_ps.begin();
  print_cons(*i, "*** phi ***");
  ++i;
  print_cons(*i, "*** phi1 ***");

  c_ps.add_con(z == 2);

  dim_type d1 = c_ps.affine_dim();
  bool ok1 = (d1 == 2);

  auto j = c_ps.begin();
  print_cons(*j, "*** phj ***");
  ++j;
  print_cons(*j, "*** phj1 ***");

  return ok && ok1;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
