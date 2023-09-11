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

#include <sstream>

// test stream input and output

bool
test01() {
  // input / output from stream with stream operators
  std::stringstream is("123");
  std::stringstream os;

  Integer i;
  is >> i;
  os << i;

  return check_print(os.str(), "123");
}

bool
test02() {
  // input / output from stream with stream operators
  std::stringstream is;
  std::stringstream os;

  Integer i;
  is >> i;
  os << i;

  return check_print(os.str(), "0");
}

bool
test03() {
  // input / output from stream with stream operators
  std::stringstream is("123 not a number");
  std::stringstream os;

  Integer i;
  is >> i;
  os << i;

  return check_print(os.str(), "123");
}

bool
test04() {
  // input / output from stream with stream operators
  std::stringstream is("not a number 123");
  std::stringstream os;

  Integer i;
  is >> i;
  os << i;

  return check_print(os.str(), "0");
}

bool
test05() {
  // input / output from stream with ascii_load and ascii_dump
  std::stringstream is("3.14 is pi");
  std::stringstream os;

  Integer i;
  is >> i;
  os << i;

  return check_print(os.str(), "3");
}

bool
test06() {
  // input / output from stream with ascii_load and ascii_dump
  std::stringstream is("123");
  std::stringstream os;

  Integer i;
  i.ascii_load(is);
  i.ascii_dump(os);

  return check_print(os.str(), "123");
}

bool
test07() {
  // input / output from stream with ascii_load and ascii_dump
  std::stringstream is;
  std::stringstream os;

  Integer i;
  i.ascii_load(is);
  i.ascii_dump(os);

  return check_print(os.str(), "0");
}

bool
test08() {
  // input / output from stream with ascii_load and ascii_dump
  std::stringstream is("123 not a number");
  std::stringstream os;

  Integer i;
  i.ascii_load(is);
  i.ascii_dump(os);

  return check_print(os.str(), "123");
}

bool
test09() {
  // input / output from stream with ascii_load and ascii_dump
  std::stringstream is("not a number 123");
  std::stringstream os;

  Integer i;
  i.ascii_load(is);
  i.ascii_dump(os);

  return check_print(os.str(), "0");
}

bool
test10() {
  // input / output from stream with ascii_load and ascii_dump
  std::stringstream is("3.14 is pi");
  std::stringstream os;

  Integer i;
  i.ascii_load(is);
  i.ascii_dump(os);

  return check_print(os.str(), "3");
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
