#include "pplite_test.hh"

bool
test01() {
  auto x = Var {0};

  MyPSet ps1 { 1, Topol::NNC, Spec_Elem::EMPTY };

  MyDisj ph { 1, Topol::NNC };
  ph.add_con(x >= 0 );
  ph.add_con(x <= 1);

  ps1.add_disjunct(ph);

  MyPSet ps2 = ps1;

  MyDisj ph2 { 1, Topol::NNC };
  ph2.add_con(x >= 2);
  ph2.add_con(x <= 3);

  ps2.add_disjunct(std::move(ph2));

  ps2.widening_assign(ps1, Widen_Impl::BHRZ03);

  MyPSet ps_expected { 1, Topol::NNC, Spec_Elem::EMPTY };
  ph2 = MyDisj{ 1, Topol::NNC };
  ph2.add_con( x >= 0 );
  ph2.add_con( x <= 1 );
  ps_expected.add_disjunct(std::move(ph2));
  ph2 = MyDisj{ 1, Topol::NNC };
  ph2.add_con( x >= 2 );
  ph2.add_con( x <= 3 );
  ps_expected.add_disjunct(std::move(ph2));
  ph2 = MyDisj{ 1, Topol::NNC };
  ph2.add_con( x >= 3 );
  ps_expected.add_disjunct(std::move(ph2));

  using namespace IO_Operators;

  print_pset(ps2, "*** ps2 widen ***");
  print_pset(ps_expected, "*** ps_expected ***");

  return ps2.geom_covers(ps_expected);
}

bool
test02() {
  auto x = Var {0};

  MyPSet ps1 { 1, Topol::CLOSED, Spec_Elem::EMPTY };

  MyDisj ph1 { 1, Topol::CLOSED };
  ph1.add_con(x == 1 );
  MyDisj ph2 { 1, Topol::CLOSED };
  ph2.add_con(x == 3);

  ps1.add_disjunct(ph1);
  ps1.add_disjunct(ph2);

  MyPSet ps2 = ps1;

  MyDisj ph3 { 1, Topol::CLOSED };
  ph3.add_con(x == 2);

  ps2.add_disjunct(std::move(ph3));

  ps2.widening_assign(ps1, Widen_Impl::BHRZ03);

  MyPSet ps_expected { 1, Topol::CLOSED, Spec_Elem::EMPTY };
  ph3 = { 1, Topol::CLOSED };
  ph3.add_con( x >= 1 );
  ph3.add_con( x <= 3 );
  ps_expected.add_disjunct(std::move(ph3));

  using namespace IO_Operators;

  print_pset(ps2, "*** ps2 widen ***");
  print_pset(ps_expected, "*** ps_expected ***");

  return ps2.geom_covers(ps_expected);
}

bool
test03() {
  Var x {0};
  Var y {1};

  MyDisj p1 {2};
  p1.add_con(x >= 2);
  p1.add_con(x <= 8);
  p1.add_con(y >= 6);
  p1.add_con(y <= 9);
  p1.add_con(x - y >= -4);
  p1.add_con(x + y <= 14);

  MyDisj p2 {2};
  p2.add_con(x >= 2);
  p2.add_con(x <= 4);
  p2.add_con(y >= 2);
  p2.add_con(y <= 4);

  MyDisj p3 {2};
  p3.add_con(x >= 6);
  p3.add_con(x <= 13);
  p3.add_con(y >= 2);
  p3.add_con(y <= 4);

  MyDisj p4 {2};
  p4.add_con(x >= 11);
  p4.add_con(x <= 13);
  p4.add_con(y >= 7);
  p4.add_con(y <= 9);

  MyPSet ps1 {2, Topol::CLOSED, Spec_Elem::EMPTY};
  ps1.add_disjunct(p1);
  ps1.add_disjunct(p2);
  ps1.add_disjunct(p3);
  ps1.add_disjunct(p4);

  MyDisj q1 {2};
  q1.add_con(x >= 1);
  q1.add_con(x <= 14);
  q1.add_con(y >= 1);
  q1.add_con(y <= 10);

  MyDisj q2 {2};
  q2.add_con(x >= 17);
  q2.add_con(x <= 19);
  q2.add_con(y >= 9);
  q2.add_con(y <= 10);

  MyDisj q3 {2};
  q3.add_con(x >= 17);
  q3.add_con(x <= 19);
  q3.add_con(y >= 2);
  q3.add_con(y <= 5);

  MyPSet ps2 {2, Topol::CLOSED, Spec_Elem::EMPTY};
  ps2.add_disjunct(q1);
  ps2.add_disjunct(q2);
  ps2.add_disjunct(q3);

  auto ps_expected = ps2;

  ps2.widening_assign(ps1, Widen_Impl::BHRZ03);

  print_pset(ps1, "*** ps2 widen ***");
  print_pset(ps_expected, "*** expected ***");

  return ps2.geom_covers(ps_expected);
}

bool
test04() {
  Var x {0};
  Var y {1};

  MyDisj p1 {2};
  p1.add_con(x + y <= 10);
  p1.add_con(-x + y >= -4);
  p1.add_con(x >= 2);
  p1.add_con(x <= 12);
  p1.add_con(x >= 3);
  p1.add_con(x <= 11);
  p1.add_con(y >= 7);
  p1.add_con(y <= 11);

  MyDisj p2 {2};
  p2.add_con(x >= 2);
  p2.add_con(x <= 4);
  p2.add_con(y >= 2);
  p2.add_con(y <= 4);

  MyDisj p3 {2};
  p3.add_con(x >= 10);
  p3.add_con(x <= 12);
  p3.add_con(y >= 2);
  p3.add_con(y <= 4);

  MyPSet ps1 {2, Topol::CLOSED, Spec_Elem::EMPTY};
  ps1.add_disjunct(p1);
  ps1.add_disjunct(p2);
  ps1.add_disjunct(p3);

  MyDisj q1 {2};
  q1.add_con(x >= 1);
  q1.add_con(x <= 13);
  q1.add_con(y >= 1);
  q1.add_con(y >= 12);

  MyPSet qs1 {2, Topol::CLOSED, Spec_Elem::EMPTY};
  qs1.add_disjunct(q1);

  auto ps_expected = qs1;

  qs1.widening_assign(ps1, Widen_Impl::BHRZ03);

  print_pset(qs1, "*** qs1 widen ***");
  print_pset(ps_expected, "*** expected ***");

  return qs1.geom_covers(ps_expected);
}

bool
test05() {
  Var x {0};
  Var y {1};

  MyDisj p1 {2};
  p1.add_con(x >= 2);
  p1.add_con(x <= 12);
  p1.add_con(x >= 3);
  p1.add_con(x <= 11);
  p1.add_con(y >= 7);
  p1.add_con(y <= 11);
  p1.add_con(x + y <= 10);
  p1.add_con(-x + y >= -4);

  MyDisj p2 {2};
  p2.add_con(x >= 2);
  p2.add_con(x <= 4);
  p2.add_con(y >= 2);
  p2.add_con(y <= 4);

  MyDisj p3 {2};
  p3.add_con(x >= 10);
  p3.add_con(x <= 12);
  p3.add_con(y >= 2);
  p3.add_con(y <= 4);

  MyPSet ps1 {2, Topol::CLOSED, Spec_Elem::EMPTY};
  ps1.add_disjunct(p1);
  ps1.add_disjunct(p2);
  ps1.add_disjunct(p3);

  MyDisj q1 {2};
  q1.add_con(x >= 1);
  q1.add_con(x <= 13);
  q1.add_con(y >= 1);
  q1.add_con(y <= 12);

  MyDisj q2 {2};
  q2.add_con(x >= 17);
  q2.add_con(x <= 19);
  q2.add_con(y >= 1);
  q2.add_con(y <= 6);

  MyDisj q3 {2};
  q3.add_con(x >= 17);
  q3.add_con(x <= 19);
  q3.add_con(y >= 10);
  q3.add_con(y <= 12);

  MyPSet qs1 {2, Topol::CLOSED, Spec_Elem::EMPTY};
  qs1.add_disjunct(q1);
  qs1.add_disjunct(q2);
  qs1.add_disjunct(q3);

  auto ps_expected = qs1;

  qs1.widening_assign(ps1, Widen_Impl::BHRZ03);

  print_pset(qs1, "*** qs1 widen ***");
  print_pset(ps_expected, "*** expected ***");

  return qs1.geom_covers(ps_expected);
}

bool
test06() {
  Var x {0};
  Var y {1};

  MyDisj p1 { 2, Topol::NNC };
  p1.add_con(x >= 1);
  p1.add_con(x <= 4);
  p1.add_con(y >= 1);
  p1.add_con(y <= 2);

  MyDisj p2 { 2, Topol::NNC };
  p2.add_con(x >= 1);
  p2.add_con(x <= 3);
  p2.add_con(y >= 3);
  p2.add_con(y <= 5);

  MyDisj p3 { 2, Topol::NNC };
  p3.add_con(x >= 6);
  p3.add_con(x <= 8);
  p3.add_con(y >= 2);
  p3.add_con(y <= 5);

  MyDisj p4 { 2, Topol::NNC };
  p4.add_con(x >= 1);
  p4.add_con(x <= 4);
  p4.add_con(y >= 6);
  p4.add_con(y <= 7);

  MyPSet ps1 { 2, Topol::NNC, Spec_Elem::EMPTY};
  ps1.add_disjunct(p1);
  ps1.add_disjunct(p2);
  ps1.add_disjunct(p3);

  MyPSet ps2 = ps1;
  ps2.add_disjunct(p4);

  ps2.widening_assign(ps1, Widen_Impl::BHRZ03);

  MyDisj p_d { 2, Topol::NNC };
  p_d.add_con(x >= 1);
  p_d.add_con(x <= 8);
  p_d.add_con(y >= 5);
  p_d.add_con(2*x + 7*y >= 51);

  MyPSet ps_expected { 2, Topol::NNC, Spec_Elem::EMPTY};
  ps_expected.add_disjunct(p1);
  ps_expected.add_disjunct(p2);
  ps_expected.add_disjunct(p3);
  ps_expected.add_disjunct(p4);
  ps_expected.add_disjunct(p_d);

  print_pset(ps2, "*** ps2 widen ***");
  print_pset(ps_expected, "*** expected ***");

  // return ps2.geom_covers(ps_expected);
  return ps_expected.geom_covers(ps2);
}
bool
test07() {
  Var A{0};
  Var B{1};

  MyDisj p{2, Topol::CLOSED};
  MyDisj q{2, Topol::CLOSED};
  MyDisj r{2, Topol::CLOSED};
  MyDisj s{2, Topol::CLOSED};
  p.add_con(A >= 1);
  p.add_con(B == 0);
  q.add_con(A >= 2);
  q.add_con(A <= 7);
  q.add_con(B == 1);
  r.add_con(A >= 3);
  r.add_con(A <= 8);
  r.add_con(B == 1);
  s.add_con(A >= 1);
  s.add_con(A <= 6);
  s.add_con(B == 1);
  MyPSet P{2, Topol::CLOSED, Spec_Elem::EMPTY};
  P.add_disjunct(p);
  P.add_disjunct(q);
  P.add_disjunct(r);
  P.add_disjunct(s);
  MyPSet Q{2, Topol::CLOSED, Spec_Elem::EMPTY};
  Q.add_disjunct(p);
  Q.add_disjunct(q);
  Q.add_disjunct(s);

  print_pset(P, "*** P widen ***");
  print_pset(Q, "*** Q widen ***");

  MyPSet old_P = P;
  P.widening_assign(Q, Widen_Impl::BHRZ03);

  print_pset(P, "*** P widen ***");
  return P.geom_covers(old_P) && P.geom_covers(Q);
}

bool
test08() {
  Var X{0};
  Var Y{1};

  MyDisj p1{2, Topol::CLOSED};
  p1.add_con(X >= 0);
  p1.add_con(Y >= 0);
  p1.add_con(X <= 2);
  p1.add_con(Y <= 1);

  MyDisj p2{2, Topol::CLOSED};
  p2.add_con(X >= 0);
  p2.add_con(Y >= 2);
  p2.add_con(X <= 1);
  p2.add_con(Y <= 3);

  MyDisj p3{2, Topol::CLOSED};
  p3.add_con(X >= 3);
  p3.add_con(Y >= 1);
  p3.add_con(X <= 4);
  p3.add_con(Y <= 3);

  MyPSet T1{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T1.add_disjunct(p1);
  T1.add_disjunct(p2);
  T1.add_disjunct(p3);

  MyDisj p4{2, Topol::CLOSED};
  p4.add_con(X >= 0);
  p4.add_con(Y >= 4);
  p4.add_con(X <= 2);
  p4.add_con(Y <= 5);

  MyPSet T2{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T2.add_disjunct(p1);
  T2.add_disjunct(p2);
  T2.add_disjunct(p3);
  T2.add_disjunct(p4);

  print_pset(T1, "*** T1 ***");
  print_pset(T2, "*** T2 ***");

  MyPSet old_T2 = T2;
  T2.widening_assign(T1, Widen_Impl::BHRZ03);

  MyDisj pd{2, Topol::CLOSED};
  pd.add_con(X >= 0);
  pd.add_con(X <= 4);
  pd.add_con(X + 2*Y >= 10);

  MyPSet known_result = old_T2;
  known_result.add_disjunct(pd);

  print_pset(T2, "*** T2 widen ***");
  print_pset(known_result, "*** known_result ***");

  return T2 == known_result
    && T2.geom_covers(old_T2)
    && T2.geom_covers(T1);
}

bool
test09() {
  Var A{0};
  Var B{1};

  MyDisj p{2, Topol::CLOSED};
  MyDisj q{2, Topol::CLOSED};
  MyDisj r{2, Topol::CLOSED};
  MyDisj s{2, Topol::CLOSED};
  p.add_con(A >= 1);
  p.add_con(B == 0);
  q.add_con(A >= 2);
  q.add_con(A <= 7);
  q.add_con(B == 1);
  r.add_con(A >= 3);
  r.add_con(A <= 8);
  r.add_con(B == 1);
  s.add_con(A >= 1);
  s.add_con(A <= 6);
  s.add_con(B == 1);
  MyPSet P{2, Topol::CLOSED, Spec_Elem::EMPTY};
  P.add_disjunct(p);
  P.add_disjunct(q);
  P.add_disjunct(r);
  P.add_disjunct(s);
  MyPSet Q{2, Topol::CLOSED, Spec_Elem::EMPTY};
  Q.add_disjunct(p);
  Q.add_disjunct(q);
  Q.add_disjunct(s);

  print_pset(P, "*** P ***");
  print_pset(Q, "*** Q ***");

  MyPSet old_P = P;
  P.widening_assign(Q, Widen_Impl::BHRZ03);

  print_pset(P, "*** P widen ***");

  return P.geom_covers(old_P) && P.geom_covers(Q);
}

bool
test10() {
  Var X{0};
  Var Y{1};

  MyDisj p1{2, Topol::CLOSED};
  p1.add_con(X >= 0);
  p1.add_con(Y >= 0);
  p1.add_con(X <= 2);
  p1.add_con(Y <= 1);

  MyDisj p3{2, Topol::CLOSED};
  p3.add_con(X >= 3);
  p3.add_con(Y >= 1);
  p3.add_con(X <= 4);
  p3.add_con(Y <= 3);

  MyDisj p4{2, Topol::CLOSED};
  p4.add_con(X >= 0);
  p4.add_con(Y >= 4);
  p4.add_con(X <= 2);
  p4.add_con(Y <= 5);

  MyPSet T1{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T1.add_disjunct(p1);
  T1.add_disjunct(p3);
  T1.add_disjunct(p4);

  MyDisj p2{2, Topol::CLOSED};
  p2.add_con(X >= 0);
  p2.add_con(Y >= 2);
  p2.add_con(X <= 1);
  p2.add_con(Y <= 3);

  MyPSet T2{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T2.add_disjunct(p1);
  T2.add_disjunct(p2);
  T2.add_disjunct(p3);
  T2.add_disjunct(p4);

  print_pset(T1, "*** T1 ***");
  print_pset(T2, "*** T2 ***");

  MyPSet old_T2 = T2;
  T2.widening_assign(T1, Widen_Impl::BHRZ03);

  MyDisj phull_T2{2, Topol::CLOSED};
  phull_T2.add_con(X >= 0);
  phull_T2.add_con(X <= 4);
  phull_T2.add_con(Y >= 0);
  phull_T2.add_con(Y <= 5);
  phull_T2.add_con(X - 2*Y <= 2);
  phull_T2.add_con(X + Y <= 7);

  MyPSet known_result{2, Topol::CLOSED, Spec_Elem::EMPTY};
  known_result.add_disjunct(phull_T2);

  print_pset(T2, "*** T2 widen ***");
  print_pset(known_result, "*** known result ***");

  return T2 == known_result
    && T2.geom_covers(old_T2)
    && T2.geom_covers(T1);
}

// This tests the first case of the widening definition when the widening
// does nothing as the lgo for the poly-hull is decreasing.
bool
test11() {
  Var X{0};
  Var Y{1};

  MyDisj p1{2, Topol::CLOSED};
  p1.add_con(Y >= 2);
  p1.add_con(Y - X <= 2);
  p1.add_con(X + Y <= 4);

  MyDisj p2{2, Topol::CLOSED};
  p2.add_con(X >= 0);
  p2.add_con(Y >= 0);
  p2.add_con(X <= 1);
  p2.add_con(Y <= 1);

  MyDisj p3{2, Topol::CLOSED};
  p3.add_con(X >= 2);
  p3.add_con(Y >= 0);
  p3.add_con(X <= 4);
  p3.add_con(Y <= 1);

  MyDisj p4{2, Topol::CLOSED};
  p4.add_con(X >= 3);
  p4.add_con(Y >= 2);
  p4.add_con(X <= 4);
  p4.add_con(Y <= 3);

  MyPSet T1{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T1.add_disjunct(p1);
  T1.add_disjunct(p2);
  T1.add_disjunct(p3);
  T1.add_disjunct(p4);

  MyDisj q1{2, Topol::CLOSED};
  q1.add_con(X >= 0);
  q1.add_con(Y >= 0);
  q1.add_con(X <= 4);
  q1.add_con(Y <= 4);

  MyDisj q2{2, Topol::CLOSED};
  q2.add_con(X >= 5);
  q2.add_con(Y >= 3);
  q2.add_con(X <= 6);
  q2.add_con(Y <= 4);

  MyDisj q3{2, Topol::CLOSED};
  q3.add_con(X >= 5);
  q3.add_con(Y >= 0);
  q3.add_con(X <= 6);
  q3.add_con(Y <= 2);

  MyPSet T2{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T2.add_disjunct(q1);
  T2.add_disjunct(q2);
  T2.add_disjunct(q3);

  print_pset(T1, "*** T1 ***");
  print_pset(T2, "*** T2 ***");

  MyPSet old_T2 = T2;
  T2.widening_assign(T1, Widen_Impl::BHRZ03);

  print_pset(T2, "*** T2 widen ***");

  return T2 == old_T2
    && T2.geom_covers(old_T2)
    && T2.geom_covers(T1);
}

// This tests the first case of the widening definition when the widening
// does nothing; the poly-hull is stable with respect to the certificate
// and the multiset ordering for this certificate is decreasing.
bool
test12() {
  Var X{0};
  Var Y{1};

  MyDisj p1{2, Topol::CLOSED};
  p1.add_con(X >= 1);
  p1.add_con(Y >= 4);
  p1.add_con(X <= 7);
  p1.add_con(Y <= 7);
  p1.add_con(X - Y <= 2);
  p1.add_con(X + Y >= 6);

  MyDisj p2{2, Topol::CLOSED};
  p2.add_con(X >= 1);
  p2.add_con(Y >= 1);
  p2.add_con(X <= 3);
  p2.add_con(Y <= 3);

  MyDisj p3{2, Topol::CLOSED};
  p3.add_con(X >= 5);
  p3.add_con(Y >= 1);
  p3.add_con(X <= 7);
  p3.add_con(Y <= 3);

  MyPSet T1{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T1.add_disjunct(p1);
  T1.add_disjunct(p2);
  T1.add_disjunct(p3);

  MyDisj q1{2, Topol::CLOSED};
  q1.add_con(X >= 0);
  q1.add_con(Y >= 0);
  q1.add_con(X <= 8);
  q1.add_con(Y <= 8);

  MyDisj q2{2, Topol::CLOSED};
  q2.add_con(X >= 10);
  q2.add_con(Y >= 6);
  q2.add_con(X <= 12);
  q2.add_con(Y <= 8);

  MyDisj q3{2, Topol::CLOSED};
  q3.add_con(X >= 10);
  q3.add_con(Y >= 0);
  q3.add_con(X <= 12);
  q3.add_con(Y <= 4);

  MyPSet T2{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T2.add_disjunct(q1);
  T2.add_disjunct(q2);
  T2.add_disjunct(q3);

  print_pset(T1, "*** T1 ***");
  print_pset(T2, "*** T2 ***");

  MyPSet old_T2 = T2;
  T2.widening_assign(T1, Widen_Impl::BHRZ03);

  print_pset(T2, "*** T2 widen ***");

  return T2 == old_T2
    && T2.geom_covers(old_T2)
    && T2.geom_covers(T1);
}

// This tests the first case of the widening definition when the widening
// of the elements of the set reduces the multiset ordering.
bool
test13() {
  Var X{0};
  Var Y{1};

  MyDisj p1{2, Topol::CLOSED};
  p1.add_con(Y >= 2);
  p1.add_con(Y <= 3);
  p1.add_con(Y - X <= 2);
  p1.add_con(X + Y <= 8);

  MyDisj p2{2, Topol::CLOSED};
  p2.add_con(X >= 0);
  p2.add_con(Y >= 0);
  p2.add_con(X <= 1);
  p2.add_con(Y <= 1);

  MyDisj p3{2, Topol::CLOSED};
  p3.add_con(X >= 5);
  p3.add_con(Y >= 0);
  p3.add_con(X <= 8);
  p3.add_con(Y <= 1);

  MyDisj p4{2, Topol::CLOSED};
  p4.add_con(X >= 7);
  p4.add_con(Y >= 4);
  p4.add_con(X <= 8);
  p4.add_con(Y <= 5);

  MyPSet T1{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T1.add_disjunct(p1);
  T1.add_disjunct(p2);
  T1.add_disjunct(p3);
  T1.add_disjunct(p4);

  MyDisj q1{2, Topol::CLOSED};
  q1.add_con(Y >= 2);
  q1.add_con(Y <= 4);
  q1.add_con(Y - X <= 2);
  q1.add_con(X + Y <= 8);

  MyPSet T2{2, Topol::CLOSED, Spec_Elem::EMPTY};
  T2.add_disjunct(q1);
  T2.add_disjunct(p2);
  T2.add_disjunct(p3);
  T2.add_disjunct(p4);

  print_gens(p1, "=== p1 ===");
  print_gens(q1, "=== q1 ===");

  print_pset(T1, "*** T1 ***");
  print_pset(T2, "*** T2 ***");

  MyPSet old_T2 = T2;

  T2.widening_assign(T1, Widen_Impl::BHRZ03);

  MyDisj r1{2, Topol::CLOSED};
  r1.add_con(Y >= 2);
  r1.add_con(Y - X <= 2);
  r1.add_con(X + Y <= 8);

  MyPSet known_result{2, Topol::CLOSED, Spec_Elem::EMPTY};
  known_result.add_disjunct(r1);
  known_result.add_disjunct(p2);
  known_result.add_disjunct(p3);
  known_result.add_disjunct(p4);

  print_pset(T2, "*** T2 widen ***");
  print_pset(known_result, "*** known_result ***");

  return T2 == known_result
    && T2.geom_covers(old_T2)
    && T2.geom_covers(T1);
}

const MyDisj&
aux1_test14(unsigned n) {
  Var x(0);
  Var y(1);

  static std::vector<MyDisj>p;
  if (p.size() < 5) {
    p.resize(5, MyDisj{2, Topol::CLOSED});
    p[2].add_con(0 <= x);
    p[2].add_con(x <= 4);
    p[2].add_con(0 <= y);
    p[2].add_con(y <= 4);
    p[1] = p[2];
    p[1].add_con(x-y <= 3);
    p[0] = p[1];
    p[0].add_con(x+y >= 1);

    p[3].add_con(0 <= x);
    p[3].add_con(x <= 8);
    p[3].add_con(0 <= y);
    p[3].add_con(y <= 8);
    p[3].add_con(x+y <= 14);
    p[3].add_con(x-y >= -6);
    p[4] = p[3];
    p[3].add_con(5*x-y >= -2);
    p[3].add_con(x+3*y >= 3);
    p[4].add_con(4*x-y >= -3);
    p[4].add_con(x+2*y >= 2);
  }

  if (n >= p.size()) {
    unsigned new_size = p.size();
    while (n >= new_size)
      new_size *= 2;
    p.resize(p.size()*2);
  }

  if (p[n].is_universe()) {
    p[n] = aux1_test14(n-4);
    p[n].affine_image(x, 2*x);
    p[n].affine_image(y, -2*y, 8);
  }

  return p[n];
}

MyPSet
aux2_test14(unsigned n) {
  MyPSet s{2, Topol::CLOSED, Spec_Elem::EMPTY};
  if (n == 0) {

    nout << "S0 = { P0 }\n";

    s.add_disjunct(aux1_test14(0));
    return s;
  }

  const int p_base = (n-1)/3*4;

  switch (n % 3) {
  case 1:

    nout << "S" << n << " = { "
         << "P" << p_base + 1 << ", "
         << "P" << p_base + 3 << " }" << endl;

    s.add_disjunct(aux1_test14(p_base + 1));
    s.add_disjunct(aux1_test14(p_base + 3));
    break;
  case 2:

    nout << "S" << n << " = { "
         << "P" << p_base + 2 << ", "
         << "P" << p_base + 3 << " }" << endl;

    s.add_disjunct(aux1_test14(p_base + 2));
    s.add_disjunct(aux1_test14(p_base + 3));
    break;
  case 0:

    nout << "S" << n << " = { "
         << "P" << p_base + 2 << ", "
         << "P" << p_base + 4 << " }" << endl;

    s.add_disjunct(aux1_test14(p_base + 2));
    s.add_disjunct(aux1_test14(p_base + 4));
    break;
  }
  return s;
}

void
aux3_test14(std::ostream& s, const Var v) {
  s << char('x' + v.id());
}

bool
test14() {
  // Install the alternate output function.
  Var::output_function.set(aux3_test14);

  MyPSet T = aux2_test14(0);

  print_pset(T, "*** T0 ***");
  // nout << "T0 = " << T << endl;

  bool converged = false;
  for (unsigned n = 1; !converged && n <= 20; ++n) {
    MyPSet Sn = aux2_test14(n);

    nout << "S" << n << " = ";
    print_pset(Sn, "");

    Sn.join_assign(T);
    Sn.widening_assign(T, Widen_Impl::BHRZ03);

    nout << "T" << n << " = ";
    print_pset(Sn, "");

    if (Sn.definitely_entails(T))
      converged = true;
    else
      swap(Sn, T);
  }

  return converged;
}

const MyDisj&
aux1_test15(unsigned n) {
  Var x(0);
  Var y(1);
  static std::vector<MyDisj> p;
  if (p.size() < 5) {
    p.resize(5, MyDisj{2, Topol::CLOSED});
    p[2].add_con(0 <= x);
    p[2].add_con(x <= 4);
    p[2].add_con(0 <= y);
    p[2].add_con(y <= 4);
    p[1] = p[2];
    p[1].add_con(x-y <= 3);
    p[0] = p[1];
    p[0].add_con(x+y >= 1);

    p[3].add_con(0 <= x);
    p[3].add_con(x <= 8);
    p[3].add_con(0 <= y);
    p[3].add_con(y <= 8);
    p[3].add_con(x+y <= 14);
    p[3].add_con(x-y >= -6);
    p[4] = p[3];
    p[3].add_con(5*x-y >= -2);
    p[3].add_con(x+3*y >= 3);
    p[4].add_con(4*x-y >= -3);
    p[4].add_con(x+2*y >= 2);
  }

  if (n >= p.size()) {
    unsigned new_size = p.size();
    while (n >= new_size)
      new_size *= 2;
    p.resize(p.size()*2);
  }

  if (p[n].is_universe()) {
    p[n] = aux1_test15(n-4);
    p[n].affine_image(x, 2*x);
    p[n].affine_image(y, -2*y, 8);
  }

  return p[n];
}

MyPSet
aux2_test15(unsigned n) {
  MyPSet s{2, Topol::CLOSED, Spec_Elem::EMPTY};
  if (n == 0) {

    nout << "S0 = { P0 }" << endl;

    s.add_disjunct(aux1_test15(0));
    return s;
  }

  const int p_base = (n-1)/3*4;

  switch (n % 3) {
  case 1:

    nout << "S" << n << " = { "
         << "P" << p_base + 1 << ", "
         << "P" << p_base + 3 << " }" << endl;

    s.add_disjunct(aux1_test15(p_base + 1));
    s.add_disjunct(aux1_test15(p_base + 3));
    break;
  case 2:

    nout << "S" << n << " = { "
         << "P" << p_base + 2 << ", "
         << "P" << p_base + 3 << " }" << endl;

    s.add_disjunct(aux1_test15(p_base + 2));
    s.add_disjunct(aux1_test15(p_base + 3));
    break;
  case 0:

    nout << "S" << n << " = { "
         << "P" << p_base + 2 << ", "
         << "P" << p_base + 4 << " }" << endl;

    s.add_disjunct(aux1_test15(p_base + 2));
    s.add_disjunct(aux1_test15(p_base + 4));
    break;
  }
  return s;
}

void
aux3_test15(std::ostream& s, const Var v) {
  s << char('x' + v.id());
}

bool
test15() {
  // Install the alternate output function.
  Var::output_function.set(aux3_test15);

  MyPSet T = aux2_test15(0);

  nout << "T0 = ";
  print_pset(T, "");

  bool converged = false;
  for (unsigned n = 1; !converged && n <= 20; ++n) {
    MyPSet Sn = aux2_test15(n);

    nout << "S" << n << " = ";
    print_pset(Sn, "");

    Sn.join_assign(T);
    Sn.widening_assign(T, Widen_Impl::BHRZ03);

    nout << "T" << n << " = ";
    print_pset(Sn, "");

    if (Sn.definitely_entails(T))
      converged = true;
    else
      swap(Sn, T);
  }

  return converged;
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
  DO_TEST(test09);
  DO_TEST(test10);
  DO_TEST(test11);
  DO_TEST(test12);
  // DO_TEST(test13);
  DO_TEST(test14);
  DO_TEST(test15);
END_MAIN
