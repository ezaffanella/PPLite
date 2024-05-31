/* PPLite: a lightweight library for convex polyhedra derived from PPL.
   Copyright (C) 2018-2024 Enea Zaffanella <enea.zaffanella@unipr.it>

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

#include "pplite-config.h"
#include "Bits.hh"
#include "ascii_dump_load.hh"

#include <cassert>
#include <iostream>

namespace pplite {

namespace {

inline Bits::Word
shifted_copy(Bits::Word src1, Bits::Word src2, dim_type start) {
  assert(0 < start && start < Bits::word_size);
  /*
    Example (8 bit word for clarity):
    src1 = [a7 a6 a5 a4 a3 a2 a1 a0]
    src2 = [b7 b6 b5 b4 b3 b2 b1 b0]
    start = 3
  */
  // right shift the leftmost `word_size - start' bits from src1
  src1 = (src1 >> start);
  // now src1 = [0 0 0 a7 a6 a5 a4 a4]
  // left shift the rightmost `start' bits from src2
  src2 = (src2 << (Bits::word_size - start));
  // now src2 = [b2 b1 b0 0 0 0 0 0]
  return src1 | src2; // [b2 b1 b0 a7 a6 a5 a4 a3]
}

inline Bits::Word
partial_copy(Bits::Word src, dim_type src_start,
             Bits::Word dst, dim_type dst_start,
             dim_type count) {
  assert(count > 0);
  constexpr auto w_size = Bits::word_size;
  constexpr auto w_ones = Bits::word_ones();
  assert(src_start + count <= w_size);
  assert(dst_start + count <= w_size);
  Bits::Word src_mask = (w_ones << src_start);
  if (src_start + count < w_size)
    src_mask &= ~(w_ones << (src_start + count));
  Bits::Word dst_mask = ~(w_ones << dst_start);
  if (dst_start + count < w_size)
    dst_mask |= (w_ones << (dst_start + count));
  Bits::Word src_masked_shifted = (src_start < dst_start)
    ? ((src & src_mask) << (dst_start - src_start))
    : ((src & src_mask) >> (src_start - dst_start));
  return (dst & dst_mask) | src_masked_shifted;
}

} // namespace

void
Bits::copy_from(dim_type dst_k,
                const Bits& src, dim_type src_k, dim_type count) {
  auto [dst_wk, dst_bk] = split(dst_k);
  auto [src_wk, src_bk] = split(src_k);

#ifndef NDEBUG
  auto min_bits = dst_k + count;
  auto min_words = min_bits / word_size;
  if (min_bits % word_size != 0)
    ++min_words;
  assert(num_words() >= min_words);
#endif // NDEBUG

  auto& dst_ws = words;
  const auto& src_ws = src.words;

  // First (partial) word.
  if (dst_bk < src_bk) {
    auto partial_count = std::min(word_size - src_bk, count);
    dst_ws[dst_wk] = partial_copy(src_ws[src_wk], src_bk,
                                  dst_ws[dst_wk], dst_bk, partial_count);
    if (count <= partial_count)
      return;
    count -= partial_count;
    ++src_wk;
    assert(src_bk + partial_count == word_size);
    src_bk = 0;
    dst_bk += partial_count;
    partial_count = std::min(word_size - dst_bk, count);
    dst_ws[dst_wk] = partial_copy(src_ws[src_wk], src_bk,
                                  dst_ws[dst_wk], dst_bk, partial_count);
    if (count <= partial_count)
      return;
    count -= partial_count;
    ++dst_wk;
    src_bk = partial_count;
    dst_bk = 0;
  } else {
    assert(dst_bk >= src_bk);
    auto partial_count = std::min(word_size - dst_bk, count);
    dst_ws[dst_wk] = partial_copy(src_ws[src_wk], src_bk,
                                  dst_ws[dst_wk], dst_bk, partial_count);
    if (count <= partial_count)
      return;
    count -= partial_count;
    ++dst_wk;
    if (dst_bk == src_bk) {
      ++src_wk;
      src_bk = 0;
      dst_bk = 0;
    } else {
      src_bk += partial_count;
      dst_bk = 0;
    }
  }
  // Now dst word is aligned.
  assert(dst_bk == 0);

  // Copy other (complete) words.
  auto count_div = count / word_size;
  if (src_bk == 0) {
    // Lucky case: also src is aligned.
    while (count_div > 0) {
      dst_ws[dst_wk] = src_ws[src_wk];
      ++src_wk;
      ++dst_wk;
      --count_div;
    }
    // Last (partial) word.
    count %= word_size;
    if (count > 0)
      dst_ws[dst_wk] = partial_copy(src_ws[src_wk], 0,
                                    dst_ws[dst_wk], 0, count);
  } else {
    assert(src_bk > 0);
    // src is misaligned: perform shifted copies.
    while (count_div > 0) {
      dst_ws[dst_wk] = shifted_copy(src_ws[src_wk],
                                    src_ws[src_wk+1], src_bk);
      ++src_wk;
      ++dst_wk;
      --count_div;
    }
    // Last word.
    count %= word_size;
    if (count > 0) {
      auto partial_count = std::min(word_size - src_bk, count);
      dst_ws[dst_wk] = partial_copy(src_ws[src_wk], src_bk,
                                    dst_ws[dst_wk], 0, partial_count);
      if (count <= partial_count)
        return;
      count -= partial_count;
      ++src_wk;
      assert(src_bk + partial_count == word_size);
      src_bk = 0;
      dst_bk += partial_count;
      assert(count <= word_size - dst_bk);
      dst_ws[dst_wk] = partial_copy(src_ws[src_wk], src_bk,
                                    dst_ws[dst_wk], dst_bk, count);
    }
  }
}

void
Bits::print_as_bits(std::ostream& os) const {
  const auto& x = *this;
  os << x[0];
  if (!x.empty()) {
    for (auto i : range(1, last()+1))
      os << ' ' << x[i];
  }
  os << "\n";
}

void
Bits::print_as_index_set(std::ostream& os) const {
  const auto& x = *this;
  auto sz = x.size();
  os << sz << " : ";
  if (sz == 0) {
    os << "{}\n";
    return;
  }
  bool first = true;
  for (auto i : x) {
    os << (first ? "{ " : ", ") << i;
    first = false;
  }
  os << " }\n";
}

bool
Bits::ascii_load(std::istream& is) {
  clear();
  dim_type num = 0;
  if (not ((is >> num) && (num >= 0)))
    return false;
  if (not ascii_load_string(is, ":"))
    return false;
  if (num == 0)
    return ascii_load_string(is, "{}");

  if (not ascii_load_string(is, "{"))
    return false;
  dim_type idx = not_a_dim();
  if (not ((is >> idx) && (idx >= 0)))
    return false;
  set(idx);
  dim_type old_idx = idx;
  for (auto idx : range(1, num)) {
    if (not ascii_load_string(is, ","))
      return false;
    if (not ((is >> idx) && (idx >= 0) && (idx > old_idx)))
      return false;
    set(idx);
    old_idx = idx;
  }
  return ascii_load_string(is, "}");
}

} // namespace pplite
