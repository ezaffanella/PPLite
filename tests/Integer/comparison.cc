/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2026 Enea Zaffanella <enea.zaffanella@unipr.it>

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

// test Integer comparison operators

bool
test01() {
  // equal signed integer
  Integer i(-313);
  bool t1 = (i == -313  );
  bool t2 = (i == -313L );
  bool t4 = (-313   == i);
  bool t5 = (-313L  == i);

  return t1 && t2 && t4 && t5;
}

bool
test02() {
  // equal unsigned integer
  Integer i(2017);
  bool t1 = (i == 2017U  );
  bool t2 = (i == 2017UL );
  bool t4 = (2017U   == i);
  bool t5 = (2017UL  == i);

  return t1 && t2 && t4 && t5;
}

bool
test03() {
  // equal another Integer
  Integer i(20917);
  Integer j(20917);
  return (i == j);
}

bool
test04() {
  // not equal signed integer
  Integer i(-313);
  bool t1 = ((i != -313  ) == false);
  bool t2 = ((i != -313L ) == false);
  bool t4 = ((-313   != i) == false);
  bool t5 = ((-313L  != i) == false);

  return t1 && t2 && t4 && t5;
}

bool
test05() {
  // not equal unsigned integer
  Integer i(21884);
  bool t1 = ((i != 21884U  ) == false);
  bool t2 = ((i != 21884UL ) == false);
  bool t4 = ((21884U   != i) == false);
  bool t5 = ((21884UL  != i) == false);

  return t1 && t2 && t4 && t5;
}

bool
test06() {
  // not equal another Integer
  Integer i(20917);
  Integer j(20917);

  return ((i != j) == false);
}

bool
test07() {
  // less than, less or equal with signed integer, unsigned integer, Integer
  Integer i(100);
  Integer j(100);
  Integer k(101);
  bool t1 = (((i <  100 ) == false) && ((i <  100L ) == false));
  bool t2 = (((i <= 100 ) == true ) && ((i <= 100L ) == true )
          && ((i <= 101 ) == true ) && ((i <= 101L ) == true ));
  bool t3 = (((i <  101 ) == true ) && ((i <  101L ) == true ));
  bool t4 = (((i <  100U) == false) && ((i <  100UL) == false));
  bool t5 = (((i <= 100U) == true ) && ((i <= 100UL) == true )
          && ((i <= 101U) == true ) && ((i <= 101UL) == true ));
  bool t6 = (((i <  101U) == true ) && ((i <  101UL) == true ));
  bool t7 =  ((i <  j   ) == false);
  bool t8 = (((i <= j   ) == true ) && ((i <= j + 1) == true ));
  bool t9 =  ((i <  k   ) == true );
  return t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9;
}

bool
test08() {
  // greater than, greater or equal with signed integer, unsigned integer, Integer
  Integer i(100);
  Integer j(100);
  Integer k( 99);
  bool t1 = (((i >  100 ) == false) && ((i >  100L ) == false));
  bool t2 = (((i >= 100 ) == true ) && ((i >= 100L ) == true )
          && ((i >= 100 ) == true ) && ((i >= 100L ) == true ));
  bool t3 = (((i >   99 ) == true ) && ((i >   99L ) == true ));
  bool t4 = (((i >  100U) == false) && ((i >  100UL) == false));
  bool t5 = (((i >= 100U) == true ) && ((i >= 100UL) == true )
          && ((i >=  99U) == true ) && ((i >=  99UL) == true ));
  bool t6 = (((i >   99U) == true ) && ((i >   99UL) == true ));
  bool t7 =  ((i >  j   ) == false);
  bool t8 = (((i >= j   ) == true ) && ((i >= j - 1) == true ));
  bool t9 =  ((i >  k   ) == true );
  return t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9;
}

bool
test09() {
  // testing the compare method
  Integer i(8);
  Integer j(-8);
  Integer k;

  bool t1 = (compare(i, j) >  0);
  bool t2 = (compare(j, i) <  0);
  bool t3 = (compare(i, k) >  0);
  bool t4 = (compare(k, i) <  0);
  bool t5 = (compare(j, k) <  0);
  bool t6 = (compare(k, j) >  0);
  bool t7 = (compare(i, i) == 0);
  bool t8 = (compare(j, j) == 0);
  bool t9 = (compare(k, k) == 0);

  return t1 && t2 && t3 && t4 && t5 && t6 && t7 && t8 && t9;
}

bool
test10() {
  // testing the sgn method
  Integer i(8);
  Integer j(-8);
  Integer k;

  bool t1 = (sgn(i) ==  1);
  bool t2 = (sgn(j) == -1);
  bool t3 = (sgn(k) ==  0);

  return t1 && t2 && t3;
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
END_MAIN
