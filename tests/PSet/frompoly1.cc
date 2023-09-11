#include "pplite_test.hh"

// Constructs the powerset of polyhedra from an empty polyhedron.
bool
test01() {
  MyDisj c_ph { 0, Spec_Elem::EMPTY, Topol::CLOSED };
  MyPSet c_ps { c_ph };

  bool ok = (c_ps.check_inv() && c_ps.is_empty() &&
             c_ps.space_dim() == 0);

  MyDisj nnc_ph { 0, Spec_Elem::EMPTY, Topol::NNC };
  MyPSet nnc_ps { nnc_ph };

  ok = ok
    && (nnc_ps.check_inv() && nnc_ps.is_empty() &&
        nnc_ps.space_dim() == 0);

  return ok;
}

// Constructs the powerset of polyhedra from a closed polyhedron.
bool
test02() {
  auto x = Var{0};
  auto y = Var{1};
  auto z = Var{2};
  auto w = Var{3};

  MyDisj c_ph { 4, Topol::CLOSED };
  c_ph.add_con(x <= 2);
  c_ph.add_con(z == 1);
  MyDisj nnc_ph { c_ph };
  nnc_ph.set_topology(Topol::NNC);

  MyPSet c_pps1 { c_ph };
  c_pps1.set_topology(Topol::CLOSED);
  MyPSet c_pps2 { 4, Spec_Elem::EMPTY, Topol::CLOSED };
  c_pps2.add_disjunct(c_ph);

  MyPSet nnc_pps1 { c_ph };
  nnc_pps1.set_topology(Topol::NNC);
  MyPSet nnc_pps2 { 4, Spec_Elem::EMPTY, Topol::NNC };
  nnc_pps2.add_disjunct(nnc_ph);

  bool ok = (c_pps1 == c_pps2 && nnc_pps1 == nnc_pps2);

  auto &c_phi = *(c_pps1.begin());
  print_cons(c_phi, "*** c_phi ***");
  auto &nnc_phi = *(nnc_pps1.begin());
  print_cons(nnc_phi, "*** nnc_phi ***");

  return ok && c_pps1.check_inv() && nnc_pps1.check_inv();
}

// Constructs the powerset of polyhedra from an nnc polyhedron.
bool
test03() {
  auto x = Var{0};
  auto y = Var{1};
  auto z = Var{2};
  auto w = Var{3};

  MyDisj nnc_ph { 4, Topol::NNC };
  nnc_ph.add_con(x <= 2);
  nnc_ph.add_con(z == 1);
  MyDisj c_ph { nnc_ph };
  c_ph.set_topology(Topol::CLOSED);

  MyPSet c_pps1 { nnc_ph };
  c_pps1.set_topology(Topol::CLOSED);
  MyPSet c_pps2 { 4, Topol::CLOSED, Spec_Elem::EMPTY };
  c_pps2.add_disjunct(c_ph);

  MyPSet nnc_pps1 { nnc_ph };
  MyPSet nnc_pps2 {4, Topol::NNC, Spec_Elem::EMPTY };
  nnc_pps2.add_disjunct(nnc_ph);

  bool ok = (c_pps1 == c_pps2 && nnc_pps1 == nnc_pps2);

  auto const &c_phi = *(c_pps1.begin());
  print_cons(c_phi, "*** c_phi ***");
  auto const &nnc_phi = *(nnc_pps1.begin());
  print_cons(nnc_phi, "*** nnc_phi ***");

  return ok && c_pps1.check_inv() && nnc_pps1.check_inv();
}

// Constructs the powerset of nnc polyhedra from a powerset of
// closed polyhedra.
bool
test04() {
  auto x = Var{0};
  auto y = Var{1};

  MyDisj ph { 2, Topol::CLOSED };
  ph.add_con(2*x == 1);
  ph.add_con(y >= 0);

  MyPSet pps_c { ph };

  MyPSet pps { pps_c };
  pps.set_topology(Topol::NNC);

  MyPSet known_pps { 2, Topol::NNC };

  known_pps.add_con(2*x == 1);
  known_pps.add_con(y >= 0);

  bool ok = (pps == known_pps);

  auto const &phi = *(pps.begin());
  print_cons(phi, "*** phi ***");

  return ok;
}

// Constructs the powerset of nnc polyhedra from a powerset of
// closed polyhedra.
bool
test05() {
  auto x = Var{0};
  auto y = Var{1};

  MyPSet pps { 2, Topol::CLOSED };

  pps.add_con(x >= 1);
  pps.add_con(x <= 1);
  pps.add_con(y >= 0);

  MyPSet pps1 { pps };

  auto const &phi1 = *(pps.begin());
  print_cons(phi1, "*** phi1 ***");

  pps.check_inv();

  auto const &phi = *(pps.begin());
  if(phi.check_inv())
    print_cons(phi, "*** phi after ok check ***");
  else return false;

  bool ok = true;

  return ok;
}

BEGIN_MAIN
  DO_TEST(test01);
  DO_TEST(test02);
  DO_TEST(test03);
  DO_TEST(test04);
  DO_TEST(test05);
END_MAIN
