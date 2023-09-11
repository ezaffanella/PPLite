#include "pplite_test.hh"

bool
test01() {
  MyPSet ps(0, Spec_Elem::EMPTY);
  Affine_Expr ae;
  Rational value;
  bool included = false;
  Gen g = point();
  bool ok = !ps.max(ae, value, &included, &g)
    && value.is_zero() && not included && g.is_point();

  if (!ok)
    return false;

  ps.add_disjunct(MyDisj(0));

  bool ok1 = ps.max(ae, value, &included, &g)
    && value.is_zero() && included && g.is_point();

  nout << (included ? "maximum" : "supremum") << " = " << value;
  nout << " @ ";
  using namespace pplite::IO_Operators;
  nout << g;
  nout << endl;

  return ok1;
}

bool
test02() {
  Var x(0);
  Cons cs;
  Affine_Expr ae(x);

  Rational max_value;
  bool max_included = false;
  Gen max_g = point();
  MyPSet ps(1, Spec_Elem::EMPTY, Topol::NNC);

  cs.clear();
  cs.push_back(x >= 3);
  cs.push_back(x < 14);
  MyDisj ph {1, Topol::NNC};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  bool ok = ps.max(ae, max_value, &max_included, &max_g)
    && max_value == Rational(14)
    && not max_included
    && max_g.is_closure_point();

  nout << max_value;
  nout << " @ ";
  using namespace pplite::IO_Operators;
  nout << max_g;
  nout << endl;

  if (!ok)
    return false;

  cs.clear();
  cs.push_back(x >= 3);
  cs.push_back(x <= 14);
  ph = MyDisj {1, Topol::NNC};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  ok = ps.max(ae, max_value, &max_included, &max_g)
    && max_value == Rational(14)
    && max_included && max_g.is_point();

  nout << max_value;
  nout << " @ ";
  nout << max_g;
  nout << endl;

  return ok;
}

bool
test03() {
  Var x(0);
  Var y(1);
  Cons cs;
  Affine_Expr ae(9*x + y);

  Rational value;
  bool included = false;
  Gen g = point();
  MyPSet ps(2, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 3);
  cs.push_back(3*x <= 14);
  cs.push_back(y >= 5);
  cs.push_back(11*y <= 87);
  MyDisj ph {2};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  bool ok = ps.max(ae, value, &included, &g)
    && value == Rational(549, 11)
    && included && g.is_point() && g.divisor() == 33;

  nout << value;
  nout << " @ ";
  using namespace pplite::IO_Operators;
  nout << g;
  nout << endl;

  if (!ok)
    return false;

  cs.clear();
  cs.push_back(x - 3*y >= 5);
  cs.push_back(x <= 28);
  cs.push_back(y >= 5);
  cs.push_back(4*y <= 85);
  ph = MyDisj{2};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  ok = ps.max(ae, value, &included, &g)
    && value == Rational(779, 3)
    && included && g.is_point() && g.divisor() == 3;

  nout << value;
  nout << " @ ";
  nout << g;
  nout << endl;

  return ok;
}

bool
test04() {
  Var x(0);
  Affine_Expr ae(x);

  Rational value;
  bool included = false;
  Gen g = point();
  MyPSet ps(1, Spec_Elem::EMPTY);

  bool ok = !ps.max(ae, value, &included, &g)
    && value.is_zero()
    && not included && g.is_point();

  if (!ok)
    return false;

  ps.add_disjunct(MyDisj(1));
  ok = !ps.max(ae, value, &included, &g)
    && value.is_zero()
    && not included && g.is_point();

  return ok;
}

bool
test05() {
  Affine_Expr ae;
  Rational value;
  bool included = false;
  Gen g = point();
  MyPSet ps(0, Spec_Elem::EMPTY);

  bool ok = !ps.min(ae, value, &included, &g)
    && value.is_zero()
    && not included && g.is_point();

  if (!ok)
    return false;

  ps.add_disjunct(MyDisj(0));

  ok = ps.min(ae, value, &included, &g)
    && value.is_zero()
    && included && g.is_point();

  nout << (included ? "minimum" : "supremum") << " = " << value;
  nout << " @ ";
  using namespace pplite::IO_Operators;
  nout << g;
  nout << endl;

  return ok;
}

bool
test06() {
  Var x(0);
  Cons cs;
  Affine_Expr ae(x);

  Rational value;
  bool min_included = false;
  Gen min_g = point();

  MyPSet ps(1, Spec_Elem::EMPTY, Topol::NNC);

  cs.clear();
  cs.push_back(2*x > 3);
  cs.push_back(x < 14);
  MyDisj ph {1, Topol::NNC};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  bool ok = ps.min(ae, value, &min_included, &min_g)
    && value == Rational(3, 2)
    && not min_included && min_g.is_closure_point()
    && min_g.divisor() == 2;

  nout << value;
  nout << " @ ";
  using namespace pplite::IO_Operators;
  nout << min_g;
  nout << endl;

  if (!ok)
    return false;

  cs.clear();
  cs.push_back(2*x >= 3);
  cs.push_back(x < 14);
  ph = MyDisj{1, Topol::NNC};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  ok = ps.min(ae, value, &min_included, &min_g)
    && value == Rational(3, 2)
    && min_included && min_g.is_point()
    && min_g.divisor() == 2;

  nout << value;
  nout << " @ ";
  nout << min_g;
  nout << endl;

  return ok;
}

bool
test07() {
  Var x(0);
  Var y(1);
  Cons cs;
  Affine_Expr ae(x + y);

  Rational value;
  bool included = false;
  Gen g = point();
  MyPSet ps(2, Spec_Elem::EMPTY);

  cs.clear();
  cs.push_back(x >= 3);
  cs.push_back(x <= 4);
  cs.push_back(y >= 5);
  cs.push_back(y <= 8);
  MyDisj ph {2};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  bool ok = ps.min(ae, value, &included, &g)
    && value == Rational(8)
    && included && g.is_point();

  nout << value;
  nout << " @ ";
  using namespace pplite::IO_Operators;
  nout << g;
  nout << endl;

  if (!ok)
    return false;

  cs.clear();
  cs.push_back(x - y >= 1);
  cs.push_back(x <= 8);
  cs.push_back(y >= 2);
  cs.push_back(y <= 10);
  ph = MyDisj{2};
  ph.add_cons(std::move(cs));
  ps.add_disjunct(std::move(ph));

  ok = ps.min(ae, value, &included, &g)
    && value == Rational(5)
    && included && g.is_point();

  nout << value;
  nout << " @ ";
  nout << g;
  nout << endl;

  return ok;
}

bool
test08() {
  Var x(0);
  Affine_Expr ae(x);

  Rational value;
  bool included = false;
  Gen g = point();
  MyPSet ps(1, Spec_Elem::EMPTY);

  bool ok = !ps.min(ae, value, &included, &g)
    && value.is_zero()
    && not included && g.is_point();

  if (!ok)
    return false;

  ps.add_disjunct(MyDisj(1));
  ok = !ps.min(ae, value, &included, &g)
    && value.is_zero()
    && not included && g.is_point();

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
  DO_TEST(test06);
  DO_TEST(test07);
  DO_TEST(test08);
END_MAIN
