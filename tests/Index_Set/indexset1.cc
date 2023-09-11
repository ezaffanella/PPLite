/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2023 Enea Zaffanella <enea.zaffanella@unipr.it>

This file is part of PPLite.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pplite_test.hh"

void print(const Index_Set& iset) {
  for (auto i : iset)
    std::cout << i << " ";
  std::cout << "\n";
}

bool test01() {
  Index_Set iset;
  iset.set(4); iset.set(18);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(4); known_iset.set(18);

  bool ok = (iset == known_iset);

  std::cout<<"row 0\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test02() {
  Index_Set iset;
  iset.set(21); iset.set(39);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(21); known_iset.set(37);

  bool ok = (iset == known_iset);

  std::cout<<"row 1\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test03() {
  Index_Set iset;
  iset.set(20); iset.set(55);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(20); known_iset.set(53);

  bool ok = (iset == known_iset);

  std::cout<<"row 2\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test04() {
  Index_Set iset;
  iset.set(0); iset.set(3); iset.set(33);
  iset.set(66);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(0); known_iset.set(3); known_iset.set(31);
  known_iset.set(64);

  bool ok = (iset == known_iset);

  std::cout<<"row 3\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test05() {
  Index_Set iset;
  iset.set(3); iset.set(6); iset.set(66);
  iset.set(71);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(3); known_iset.set(6); known_iset.set(64);
  known_iset.set(69);

  bool ok = (iset == known_iset);

  std::cout<<"row 4\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test06() {
  Index_Set iset;
  iset.set(3); iset.set(66); iset.set(67);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(3); known_iset.set(64); known_iset.set(65);

  bool ok = (iset == known_iset);

  std::cout<<"row 5\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test07() {
  Index_Set iset;
  iset.set(3); iset.set(17); iset.set(66);
  iset.set(134);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(3); known_iset.set(17); known_iset.set(64);
  known_iset.set(130);

  bool ok = (iset == known_iset);

  std::cout<<"row 6\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test08() {
  Index_Set iset;
  iset.set(3); iset.set(17); iset.set(66);
  iset.set(114);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(3); known_iset.set(17); known_iset.set(64);
  known_iset.set(112);

  bool ok = (iset == known_iset);

  std::cout<<"fake row 6\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test09() {
  Index_Set iset;
  iset.set(3); iset.set(5); iset.set(11);
  iset.set(66); iset.set(68); iset.set(133);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(3); known_iset.set(5); known_iset.set(11);
  known_iset.set(64); known_iset.set(66); known_iset.set(129);

  bool ok = (iset == known_iset);

  std::cout<<"row 7\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test10() {
  Index_Set iset;
  iset.set(0); iset.set(9); iset.set(33);
  iset.set(53); iset.set(56); iset.set(146);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(0); known_iset.set(9); known_iset.set(31);
  known_iset.set(51); known_iset.set(54); known_iset.set(134);

  bool ok = (iset == known_iset);

  std::cout<<"row 8\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test11() {
  Index_Set iset;
  iset.set(0); iset.set(33); iset.set(42);
  iset.set(53); iset.set(56); iset.set(154);

  Index_Set neg;
  neg.set(24); neg.set(25); neg.set(131); neg.set(132);
  neg.set(138); neg.set(139); neg.set(140);  neg.set(141);
  neg.set(142); neg.set(143); neg.set(144); neg.set(145);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(0); known_iset.set(31); known_iset.set(40);
  known_iset.set(51); known_iset.set(54); known_iset.set(142);

  bool ok = (iset == known_iset);

  std::cout<<"row 9\n";
  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test12() {
  Index_Set iset;
  iset.set(65); iset.set(134);

  Index_Set neg;
  neg.set(0);
  neg.set(129);

  std::cout<<"iset:  ";
  print(iset);

  std::cout<<"neg:  ";
  print(neg);

  iset.remove_all(neg.begin(), neg.end());

  Index_Set known_iset;
  known_iset.set(64); known_iset.set(132);

  bool ok = (iset == known_iset);

  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(known_iset);
  std::cout<<"\n";

  return ok;
}

bool test13() {
  const ulong in[] = { 2, 36, 130, 199, 206 };

  const ulong neg[] = {
    13, 61, 64, 65, 136, 138, 139, 141, 145, 146, 147, 148, 149, 150,
    156, 157, 158, 159, 160, 161, 162, 175, 179
  };

  const ulong known[] = { 2, 35, 126, 176, 183 };

  Index_Set iset;
  for (auto i : in)
    iset.set(i);

  Index_Set iset_neg;
  for (auto i : neg)
    iset_neg.set(i);

  Index_Set iset_known;
  for (auto i : known)
    iset_known.set(i);

  std::cout<<"iset:  ";
  print(iset);

  std::cout<<"neg:   ";
  print(iset_neg);

  iset.remove_all(iset_neg.begin(), iset_neg.end());

  bool ok = (iset == iset_known);

  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(iset_known);
  std::cout<<"\n";

  return ok;
}

bool test14() {
  const ulong in[] = { 130 };

  const ulong neg[] = {
    0, 1, 64
  };

  const ulong known[] = { 127 };

  Index_Set iset;
  for (auto i : in)
    iset.set(i);

  Index_Set iset_neg;
  for (auto i : neg)
    iset_neg.set(i);

  Index_Set iset_known;
  for (auto i : known)
    iset_known.set(i);

  std::cout<<"iset:  ";
  print(iset);

  std::cout<<"neg:   ";
  print(iset_neg);

  iset.remove_all(iset_neg.begin(), iset_neg.end());

  bool ok = (iset == iset_known);

  std::cout<<"iset:  ";
  print(iset);
  std::cout<<"known: ";
  print(iset_known);
  std::cout<<"\n";

  return ok;
}

bool test15() {
  Bits bs;
  bs.set(1024);
  bs.reset_from(1024);
  return bs.check_inv();
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
  DO_TEST(test13);
  DO_TEST(test14);
  DO_TEST(test15);
END_MAIN
