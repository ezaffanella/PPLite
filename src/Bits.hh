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

#ifndef pplite_Bits_hh
#define pplite_Bits_hh 1

#include "globals.hh"
#include "utils.hh"

#include "Local_Stats.hh"

#if PPLITE_LOW_LEVEL_COUNTERS
#include "Low_Level_Stats.hh"
#endif

#include <cassert>
#include <iosfwd>
#include <limits>
#include <vector>

namespace pplite {

class Bits {
public:
  Bits() = default;
  Bits(const Bits& y) = default;
  Bits(Bits&& y) = default;
  Bits& operator=(const Bits& y) = default;
  Bits& operator=(Bits&& y) = default;
  ~Bits() = default;

  explicit Bits(dim_type i)
    : Bits() { set(i); }
  explicit Bits(std::pair<dim_type, dim_type> p)
    : Bits(std::max(p.first, p.second)) {
    set(std::min(p.first, p.second));
  }
  template <typename Iter>
  Bits(Iter first, Iter last) : Bits() {
    for ( ; first != last; ++first)
      set(*first);
  }

  bool check_inv() const {
    return num_words() == 0 || words.back() != 0;
  }


  using Word = unsigned long;
  // These depend on chosen Word implementation type.
  static dim_type ctz(Word w) { return __builtin_ctzl(w); }
  static dim_type clz(Word w) { return __builtin_clzl(w); }
  static dim_type popcount(Word w) { return __builtin_popcountl(w); }
  // End of things depending on Word implementation type.
  static constexpr dim_type word_size = std::numeric_limits<Word>::digits;
  static constexpr dim_type end_pos = std::numeric_limits<dim_type>::max();
  static constexpr Word word_ones() {
    return std::numeric_limits<Word>::max();
  }

  using Words = std::vector<Word>;
  Words& impl() { return words; }
  const Words& impl() const { return words; }
  dim_type num_words() const { return num_rows(words); }
  dim_type num_bits() const { return num_words() * word_size; }

  // Provide a read-only bidirectional iterator on the bits that are set
  class const_iterator {
  public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = dim_type;
    using pointer = const value_type*;
    using reference = const value_type&;
    using difference_type = std::ptrdiff_t;

  private:
    const Bits* bits_ptr;
    value_type pos;
    static constexpr value_type end_pos = Bits::end_pos;

  public:
    const_iterator() = default;
    const_iterator(const const_iterator&) = default;
    const_iterator(const_iterator&&) = default;
    const_iterator& operator=(const const_iterator&) = default;
    const_iterator& operator=(const_iterator&&) = default;
    ~const_iterator() = default;

    const_iterator(const Bits* is, bool get_first)
      : bits_ptr(is),
        pos(get_first ? is->first() : end_pos) { }

    value_type operator*() const {
      assert(pos != end_pos);
      return pos;
    }

    const_iterator& operator++() {
      assert(pos != end_pos);
      pos = bits_ptr->next(pos);
      return *this;
    }

    const_iterator& operator--() {
      assert(pos != bits_ptr->first());
      pos = (pos == end_pos)
        ? bits_ptr->last()
        : bits_ptr->prev(pos);
      return *this;
    }

    bool operator==(const const_iterator& j) const {
      assert(bits_ptr == j.bits_ptr);
      return pos == j.pos;
    }

    bool operator!=(const const_iterator& j) const {
      return !(*this == j);
    }

  }; // class const_iterator

  NOTHROW_DEFAULT_AND_MOVES(const_iterator);

  const_iterator begin() const { return const_iterator(this, true); }
  const_iterator end() const { return const_iterator(this, false); }

  // Single bit getters.
  bool test(dim_type k) const {
    const auto ids = split(k);
    if (ids.wk >= num_words())
      return false;
    return words[ids.wk] & get_mask(ids.bk);
  }
  bool operator[](dim_type k) const { return test(k); }

  bool empty() const {
    assert(check_inv());
    return num_words() == 0;
  }

  bool intersects(const Bits& y) const {
#if PPLITE_LOW_LEVEL_COUNTERS
    LLOp_Stats::bump_sat_count();
#endif
    const auto& x = *this;
    auto min_nw = std::min(x.num_words(), y.num_words());
    for (auto i : range(min_nw)) {
      if (x.words[i] & y.words[i])
        return true;
    }
    return false;
  }

  dim_type count_ones() const {
#if PPLITE_LOW_LEVEL_COUNTERS
    LLOp_Stats::bump_sat_count();
#endif
    auto res = 0;
    for (auto w : words)
      res += popcount(w);
    return res;
  }

  dim_type count_ones_in_union(const Bits& y) const {
#if PPLITE_LOW_LEVEL_COUNTERS
    LLOp_Stats::bump_sat_count();
#endif
    const Bits& x = *this;
    auto res = 0;
    auto x_nw = x.num_words();
    auto y_nw = y.num_words();
    auto min_nw = std::min(x_nw, y_nw);
    auto max_nw = std::max(x_nw, y_nw);
    for (auto i : range(min_nw)) {
      auto w = x.words[i] | y.words[i];
      res += popcount(w);
    }
    const auto& tail_ws = (x_nw == max_nw) ? x.words : y.words;
    for (auto i : range(min_nw, max_nw))
      res += popcount(tail_ws[i]);
    return res;
  }

  dim_type size() const { return count_ones(); }

  void print_as_bits(std::ostream& os) const;
  void print_as_index_set(std::ostream& os) const;

  void print(std::ostream& os, bool as_index_set = false) const {
    if (as_index_set)
      print_as_index_set(os);
    else
      print_as_bits(os);
  }

  void ascii_dump(std::ostream& os) const { print(os, true); }
  bool ascii_load(std::istream& is);

  // Modifiers
  void swap(Bits& y) noexcept {
    using std::swap;
    swap(words, y.words);
  }
  void clear() { words.clear(); }

  // Single bit setters.
  void set(dim_type k) {
    const auto ids = split(k);
    if (ids.wk < num_words())
      words[ids.wk] |= get_mask(ids.bk);
    else {
      words.resize(ids.wk + 1, 0);
      words.back() = get_mask(ids.bk);
    }
    assert(check_inv());
  }
  void reset(dim_type k) {
    const auto ids = split(k);
    const auto nw = num_words();
    if (ids.wk >= nw)
      return;
    words[ids.wk] &= ~get_mask(ids.bk);
    if (ids.wk == nw - 1)
      rtrim();
    assert(check_inv());
  }

  // Prefix setter (k excluded).
  void set_until(dim_type k) {
    assert(k >= 0);
    if (k == 0)
      return;
    // Decrement, so k is *included*
    --k;
    const auto ids = split(k);
    const auto old_nw = num_words();
    if (old_nw > ids.wk) {
      for (auto i : range(ids.wk))
        words[i] = word_ones();
      words[ids.wk] |= get_ones(ids.bk + 1);
    } else {
      assert(old_nw <= ids.wk);
      words.assign(ids.wk + 1, word_ones());
      words.back() = get_ones(ids.bk + 1);
    }
    assert(check_inv());
  }

  // Suffix clearer (k included).
  void reset_from(dim_type k) {
    const auto nw = num_words();
    const auto ids = split(k);
    if (ids.wk >= nw)
      return;
    if (ids.bk == 0)
      words.resize(ids.wk);
    else {
      words[ids.wk] &= get_ones(ids.bk);
      words.resize(ids.wk + 1);
    }
    rtrim();
    assert(check_inv());
  }

  // Prefix setter (k excluded).
  void complement_until(dim_type k) {
    assert(k >= 0);
    if (k == 0)
      return;
    // Decrement, so k is *included*
    --k;
    const auto ids = split(k);
    const auto old_nw = num_words();
    const auto min_nw = std::min(ids.wk, old_nw);
    for (auto i : range(min_nw))
      words[i] = ~words[i];
    if (old_nw > ids.wk)
      words[ids.wk] ^= get_ones(ids.bk + 1);
    else {
      words.resize(ids.wk + 1, word_ones());
      words.back() = get_ones(ids.bk + 1);
    }
    rtrim();
    assert(check_inv());
  }

  static Bits get_union(const Bits& x, const Bits& y) {
#if PPLITE_LOW_LEVEL_COUNTERS
    LLOp_Stats::bump_sat_union();
#endif
    assert(x.check_inv() && y.check_inv());
    const auto x_nw = x.num_words();
    const auto y_nw = y.num_words();
    const auto min_nw = std::min(x_nw, y_nw);
    const auto max_nw = std::max(x_nw, y_nw);
    Bits res;
    res.words.reserve(max_nw);
    for (auto i : range(min_nw))
      res.words.push_back(x.words[i] | y.words[i]);
    const auto& tail_ws = (x_nw == max_nw) ? x.words : y.words;
    for (auto i : range(min_nw, max_nw))
      res.words.push_back(tail_ws[i]);
    assert(res.check_inv());
    return res;
  }

  void operator|=(const Bits& y) {
#if PPLITE_LOW_LEVEL_COUNTERS
    LLOp_Stats::bump_sat_union();
#endif
    auto& x = *this;
    assert(x.check_inv() && y.check_inv());
    const auto x_nw = x.num_words();
    const auto y_nw = y.num_words();
    const auto min_nw = std::min(x_nw, y_nw);
    for (auto i : range(min_nw))
      x.words[i] |= y.words[i];
    if (x_nw < y_nw) {
      x.words.reserve(y_nw);
      for (auto i : range(x_nw, y_nw))
        x.words.push_back(y.words[i]);
    }
    assert(x.check_inv());
  }

  void operator&=(const Bits& y) {
    auto& x = *this;
    const auto nw = std::min(x.num_words(), y.num_words());
    x.words.resize(nw);
    for (auto i : range(nw))
      x.words[i] &= y.words[i];
    x.rtrim();
  }

  void operator-=(const Bits& y) {
    auto& x = *this;
    const auto nw = std::min(x.num_words(), y.num_words());
    for (auto i : range(nw))
      x.words[i] &= ~y.words[i];
    x.rtrim();
  }

  void operator>>=(dim_type shift) {
    assert(shift >= 0);
    if (shift == 0 || empty())
      return;
    // Push front (whole) words
    auto ids = split(shift);
    if (ids.bk > 0)
      ++ids.wk;
    words.insert(words.begin(), ids.wk, 0);
    if (ids.bk == 0)
      return;
    // Shift back some bits
    auto& x = *this;
    auto dst_k = 0;
    auto src_k = word_size - ids.bk;
    auto count = num_bits() - src_k;
    x.copy_from(dst_k, x, src_k, count);
    // Clear trailing bits
    x.reset_from(dst_k + count);
  }

  template <typename Iter>
  void remove_all(Iter rem_first, Iter rem_last);
  template <typename Indices>
  void remove_all(const Indices& indices) {
    remove_all(std::begin(indices), std::end(indices));
  }
  // For debugging purposes.
  template <typename Iter>
  Bits remove_all_aux(Iter rem_first, Iter rem_last) const;

  // Implementation helpers

  void copy_from(dim_type dst_pos,
                 const Bits& src, dim_type src_pos, dim_type count);

  dim_type find_from_word(dim_type wk) const {
    assert(wk >= 0);
    const auto nw = num_words();
    for (auto i : range(wk, nw)) {
      if (words[i] != 0)
        return (i * word_size) + lowest_bit(words[i]);
    }
    return end_pos;
  }

  dim_type first() const { return find_from_word(0); }

  dim_type next(dim_type k) const {
    ++k;
    const auto ids = split(k);
    const auto nw = num_words();
    if (ids.wk >= nw)
      return end_pos;
    auto rest = (words[ids.wk] >> ids.bk);
    if (rest == 0)
      return find_from_word(ids.wk + 1);
    else
      return k + lowest_bit(rest);
  }

  dim_type last() const {
    auto wk = num_words();
    while (wk > 0) {
      --wk;
      if (words[wk] != 0)
        return wk * word_size + highest_bit(words[wk]);
    }
    return end_pos;
  }

  dim_type prev(dim_type k) const {
    if (k == 0)
      return end_pos;
    --k;
    const auto nw = num_words();
    if (nw == 0)
      return end_pos;
    auto wk = word_idx(k);
    if (wk >= nw)
      wk = nw - 1;
    else {
      auto mask = (word_ones() >> (word_size - 1 - bit_idx(k)));
      auto word = (words[wk] & mask);
      if (word != 0)
        return wk * word_size + highest_bit(word);
    }

    for (auto i : bwd_range(wk)) {
      if (words[i] != 0)
        return i * word_size + highest_bit(words[i]);
    }
    return end_pos;
  }

private:
  Words words;

  static dim_type lowest_bit(Word b) { return ctz(b); }
  static dim_type highest_bit(Word b) { return word_size - 1 - clz(b); }

  static dim_type word_idx(dim_type k) { return k / word_size; }
  static dim_type bit_idx(dim_type k) { return k % word_size; }
  struct Ids { dim_type wk; dim_type bk; };
  static Ids split(dim_type k) { return { word_idx(k), bit_idx(k) }; }

  static Word get_mask(dim_type k) {
    assert(0 <= k && k < word_size);
    return Word{1} << k;
  }
  static Word get_ones(dim_type k) {
    assert(0 <= k && k <= word_size);
    return word_ones() >> (word_size - k);
  }

  void rtrim() {
    auto nw = num_words();
    while (nw > 0 && words[nw - 1] == 0)
      --nw;
    words.resize(nw);
  }

}; // class Bits

NOTHROW_DEFAULT_AND_MOVES(Bits);

inline void
swap(Bits& x, Bits& y) noexcept { x.swap(y); }

inline Bits
operator&(const Bits& x, const Bits& y) {
  auto res = x;
  res &= y;
  return res;
}

inline Bits
operator|(const Bits& x, const Bits& y) {
  auto res = x;
  res |= y;
  return res;
}

inline bool
operator==(const Bits& x, const Bits& y) {
  assert(x.check_inv() && y.check_inv());
  return x.impl() == y.impl();
}

inline bool
operator!=(const Bits& x, const Bits& y) {
  return !(x == y);
}

inline bool
operator<(const Bits& x, const Bits& y) {
  assert(x.check_inv() && y.check_inv());
  return x.impl() < y.impl();
}

inline bool
subset_eq(const Bits& x, const Bits& y) {
#if PPLITE_LOW_LEVEL_COUNTERS
  LLOp_Stats::bump_sat_incl();
#endif
  assert(x.check_inv() && y.check_inv());
  const auto x_nw = x.num_words();
  const auto y_nw = y.num_words();
  if (x_nw > y_nw)
    return false;
  const auto& x_ws = x.impl();
  const auto& y_ws = y.impl();
  for (auto i : range(x_nw)) {
    if ((x_ws[i] & ~y_ws[i]) != 0)
      return false;
  }
  return true;
}

inline bool
subset_ne(const Bits& x, const Bits& y) {
#if PPLITE_LOW_LEVEL_COUNTERS
  LLOp_Stats::bump_sat_incl();
#endif
  assert(x.check_inv() && y.check_inv());
  const auto x_nw = x.num_words();
  const auto y_nw = y.num_words();
  if (x_nw > y_nw)
    return false;
  const auto& x_ws = x.impl();
  const auto& y_ws = y.impl();
  bool not_equal = (y_nw > x_nw);
  for (auto i : range(x_nw)) {
    if ((x_ws[i] & ~y_ws[i]) != 0)
      return false;
    if (!not_equal && x_ws[i] != y_ws[i])
      not_equal = true;
  }
  return not_equal;
}

template <typename Iter>
Bits
Bits::remove_all_aux(Iter rem_first, Iter rem_last) const {
  const auto& x = *this;
  Bits res;
  if (x.empty())
    return res;
  if (rem_first == rem_last) {
    res = x;
    return res;
  }
  auto i = x.begin(), i_end = x.end();
  auto j = rem_first, j_end = rem_last;
  dim_type removed = 0;
  while (i != i_end && j != j_end) {
    auto idx_i = *i;
    auto idx_j = *j;
    if (idx_i < idx_j) {
      res.set(idx_i - removed);
      ++i;
      continue;
    } else if (idx_i == idx_j) {
      ++i;
      ++j;
      ++removed;
    } else {
      ++j;
      ++removed;
    }
  }
  for ( ; i != i_end; ++i) {
    res.set(*i - removed);
  }
  return res;
}

template <typename Iter>
void
Bits::remove_all(Iter rem_first, Iter rem_last) {
  auto& x = *this;
  if (x.empty())
    return;
  if (rem_first == rem_last)
    return;
  dim_type x_last = x.last();
  assert(x_last != end_pos);
  // No need to copy first chunk.
  dim_type dst_pos = *rem_first;
  dim_type src_pos = dst_pos + 1;
  if (src_pos > x_last) {
    x.reset_from(dst_pos);
    return;
  }
  ++rem_first;
  for ( ; rem_first != rem_last; ++rem_first) {
    auto rem_pos = std::min(*rem_first, x_last + 1);
    assert(rem_pos >= src_pos);
    auto count = rem_pos - src_pos;
    if (count > 0) {
      x.copy_from(dst_pos, x, src_pos, count);
      dst_pos += count;
      src_pos += count;
    }
    ++src_pos;
    if (src_pos > x_last) {
      x.reset_from(dst_pos);
      return;
    }
  }
  // Copy suffix (if any).
  assert(x_last >= src_pos);
  auto count = (x_last + 1) - src_pos;
  if (count > 0)
    x.copy_from(dst_pos, x, src_pos, count);
  // Clear trailing bits.
  x.reset_from(dst_pos + count);
}

// Index_Set and NS_Rows are just aliases for Bits and a vector of them
using Index_Set = Bits;
using NS_Rows = std::vector<Index_Set>;

} // namespace pplite

#endif // !defined(pplite_Bits_hh)

