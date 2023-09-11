#include "pplite_test.hh"

// Powerset of C polyhedra: difference_assign().
bool
test01() {
  auto x = Var {0};
  auto y = Var {1};

  MyPSet c_ps1 {2, Spec_Elem::UNIVERSE, Topol::CLOSED};
  MyPSet c_ps2 {2, Topol::CLOSED};
  Cons cs;
  cs.emplace_back(x >= 0);
  cs.emplace_back(x <= 1);
  cs.emplace_back(y >= 0);
  cs.emplace_back(y <= 1);
  c_ps2.add_cons(std::move(cs));

  using namespace IO_Operators;
  c_ps1.difference_assign(c_ps2);

  print_pset(c_ps1, "*** c_ps1 ***");

  return true;
}

// Powerset of C polyhedra: difference_assign().
bool
test02() {
  auto x = Var {0};
  auto y = Var {1};

  MyPSet c_ps1 {2, Spec_Elem::UNIVERSE, Topol::CLOSED};
  MyPSet c_ps2 {2, Topol::CLOSED};
  Cons cs;
  cs.emplace_back(x >= 0);
  cs.emplace_back(x <= 1);
  cs.emplace_back(y >= 0);
  cs.emplace_back(y <= 1);
  c_ps2.add_cons(std::move(cs));

  using namespace IO_Operators;
  c_ps1.difference_assign(c_ps2);

  print_pset(c_ps1, "*** c_ps1 ***");

  return true;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
// DO_TEST(test03);
END_MAIN
