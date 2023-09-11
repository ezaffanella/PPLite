#include "pplite_test.hh"

namespace {

bool
test01() {
  Var x{0};
  MyPSet c_ps {1, Spec_Elem::EMPTY, Topol::CLOSED};
  Cons cs;

  cs.emplace_back(x >= 0);
  cs.emplace_back(x <= 2);
  MyDisj ph { 1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  cs.clear();
  cs.emplace_back(x >= 1);
  cs.emplace_back(x <= 3);
  ph = MyDisj{1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps.add_disjunct(std::move(ph));

  c_ps.collapse(1);

  cs.clear();
  cs.emplace_back(x >= 0);
  cs.emplace_back(x <= 3);
  MyPSet c_ps_expected {1, Spec_Elem::EMPTY };
  ph = {1, Topol::CLOSED };
  ph.add_cons(std::move(cs));
  c_ps_expected.add_disjunct(std::move(ph));

  using namespace IO_Operators;
  print_pset(c_ps, "*** c_ps ***");
  print_pset(c_ps_expected, "*** c_ps_expected ***");

  bool ok = c_ps.definitely_entails(c_ps_expected);
  bool ok1 = c_ps_expected.definitely_entails(c_ps);
  bool ok2 = (c_ps.size() == 1);

  return ok && ok1 && ok2;
}

bool
test02() {
  Var x(0);
  MyDisj ph1(1, Spec_Elem::UNIVERSE);
  ph1.add_con(x >= 1);
  MyDisj ph2(1, Spec_Elem::UNIVERSE);
  ph2.add_con(x <= 0);
  MyPSet ps1 { ph1 };
  MyPSet ps2 { ph2 };
  ps1.intersection_assign(ps2);
  ps1.collapse(1);
  return ps1.size() == 0;
}

} // namespace

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
END_MAIN
