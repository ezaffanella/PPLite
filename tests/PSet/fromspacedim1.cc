#include "pplite_test.hh"

bool
test01() {
  MyPSet ps(0, Spec_Elem::EMPTY);

  bool ok = (ps.check_inv() && ps.is_empty() && ps.space_dim() == 0);

  return ok;
}

bool
test02() {
  MyPSet ps(0, Spec_Elem::UNIVERSE);

  bool ok = (ps.check_inv() && ps.is_universe() && ps.space_dim() == 0);

  return ok;
}

bool
test03() {
  MyPSet ps(4, Spec_Elem::EMPTY);

  bool ok = (ps.check_inv() && ps.is_empty() && ps.space_dim() == 4);

  return ok;
}

bool
test04() {
  MyPSet ps(4, Spec_Elem::UNIVERSE);

  bool ok = (ps.check_inv() && ps.is_universe() && ps.space_dim() == 4);

  return ok;
}

bool
test05() {
  MyPSet ps(4);

  bool ok = (ps.check_inv() && ps.is_universe() && ps.space_dim() == 4);

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
END_MAIN
