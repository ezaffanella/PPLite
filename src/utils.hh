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

#ifndef pplite_utils_hh
#define pplite_utils_hh 1

#include "globals.hh"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <functional>
#include <list>
#include <vector>

namespace pplite {

// This helper function returns a dim_type value:
// used to avoid warnings about signed/unsigned comparisons.
template <typename Container>
inline dim_type
num_rows(const Container& rows) { return rows.size(); }

// range: enables for-range loops on a range of dim_type values.
// Typical usages:
//   for (auto i : range(20)) ... iterates from 0 to 19
//   for (auto i : range(10, 20)) ... iterates from 10 to 19
struct range {
  const dim_type first;
  const dim_type last;

  range(dim_type first_, dim_type last_) noexcept
    : first(first_), last(last_) {
    assert(0 <= first && first <= last);
  }

  explicit range(dim_type last) noexcept : range(0, last) {}

  // A *pseudo* (read-only, forward) iterator.
  // Note: it deliberately does NOT implement the iterator concept
  // (e.g., it has no post-increment operator).
  // Only meant to be used *implicitly*, in code such as
  //   for (auto i : range(10)) ...
  struct pseudo_iter {
    dim_type pos;
    explicit pseudo_iter(dim_type p = 0) noexcept : pos(p) {}
    dim_type operator*() const { return pos; }
    void operator++() { ++pos; } // void return type is meant.
    bool operator==(const pseudo_iter& y) const { return pos == y.pos; }
    bool operator!=(const pseudo_iter& y) const { return pos != y.pos; }
  }; // struct pseudo_iter

  pseudo_iter begin() const { return pseudo_iter(first); }
  pseudo_iter end() const { return pseudo_iter(last); }
}; // struct range

// Same as range, but iterations go backwards from last-1 to first.
struct bwd_range : public range {
  using range::range; // inheriting ctors

  // A *pseudo* (read-only, backward) iterator.
  struct bwd_pseudo_iter : public range::pseudo_iter {
    using pseudo_iter::pseudo_iter; // inheriting ctors
    dim_type operator*() { return pos-1; }
    void operator++() { --pos; } // void return type is meant
  }; // struct bwd_pseudo_iter

  bwd_pseudo_iter begin() const { return bwd_pseudo_iter(last); }
  bwd_pseudo_iter end() const { return bwd_pseudo_iter(first); }
}; // struct bwd_range

// dim_range: enables for-range loops on a type having space_dim().
// Typical usage:
//   for (auto i : dim_range(c)) ... iterates from 0 to c.space_dim() - 1
struct dim_range : public range {
  template <typename T>
  explicit dim_range(const T& t) noexcept
    : range(0, t.space_dim()) {}
}; // struct dim_range

// Same as above, but iterating backwards.
struct bwd_dim_range : public bwd_range {
  template <typename T>
  explicit bwd_dim_range(const T& t) noexcept
    : bwd_range(0, t.space_dim()) {}
}; // struct bwd_dim_range

// Enables for-range loops on the *indices* of a random-access container
// (e.g., Cons, Gens, Sat, Vars, ...) having method size()
// Usage: for (auto i : index_range(cs)) ...
struct index_range : public range {
  template <typename Container>
  explicit index_range(const Container& c) noexcept
    : range(0, num_rows(c)) {}
}; // struct index_range

// Same as above, but iterating backwards.
struct bwd_index_range : public bwd_range {
  template <typename Container>
  explicit bwd_index_range(const Container& c) noexcept
    : bwd_range(0, num_rows(c)) {}
}; // struct bwd_index_range


template <typename Row, typename Iter>
inline void
erase_using_sorted_indices(std::vector<Row>& rows,
                           Iter first, Iter last) {
  if (first == last)
    return;

  // `idx_unused' and `idx' start with value `*first' (instead of `0'),
  // because the loop would just increment them in the preceding iterations.
  using Index = typename std::iterator_traits<Iter>::value_type;
  const Index idx_end = rows.size();
  auto idx_unused = *first;
  assert(0 <= idx_unused);
  if (idx_unused >= idx_end)
    return;
  auto idx = idx_unused + 1;
  ++first;
  while (first != last && idx != idx_end) {
    auto new_idx = *first;
    assert(idx_unused < idx && idx <= new_idx);
    if (idx == new_idx) {
      // Don't swap row idx.
      ++first;
      ++idx;
    } else {
      // Keep (i.e., swap to its correct place) row idx.
      using std::swap;
      swap(rows[idx_unused], rows[idx]);
      ++idx_unused;
      ++idx;
    }
  }
  // Move up the remaining rows, if any.
  for ( ; idx < idx_end; ++idx) {
    using std::swap;
    swap(rows[idx_unused], rows[idx]);
    ++idx_unused;
  }
  // The trailing rows have to be erased.
  rows.resize(idx_unused);
}

template <typename Row, typename Indices>
inline void
erase_using_sorted_indices(std::vector<Row>& rows, const Indices& indices) {
  erase_using_sorted_indices(rows, std::begin(indices), std::end(indices));
}

template <typename Range, typename Pred>
inline bool
all_of(const Range& range, Pred pred) {
  return std::all_of(std::begin(range), std::end(range), pred);
}

template <typename Range, typename Pred>
inline bool
any_of(const Range& range, Pred pred) {
  return std::any_of(std::begin(range), std::end(range), pred);
}

template <typename Range, typename Pred>
inline bool
none_of(const Range& range, Pred pred) {
  return std::none_of(std::begin(range), std::end(range), pred);
}

inline size_t
hash_size(size_t val) {
  std::hash<size_t> hash_fn;
  return hash_fn(val);
}

/* Boost code (GPL compatible) */
inline void
hash_combine(size_t& lhs, size_t rhs) {
  lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
}

/* Adapter for a sequence of elements into a sequence of Val_Type. */
template <typename Dst_Type>
struct SequenceAdapter {
  using value_type = Dst_Type;
  using List = typename std::list<const value_type*>;
  using list_const_it = typename List::const_iterator;

  struct const_iterator {
    using value_type = Dst_Type;
    using pointer = const value_type*;
    using reference = const value_type&;
    using difference_type = typename list_const_it::difference_type;
    using iterator_category = typename list_const_it::iterator_category;

    list_const_it it;

    const_iterator(list_const_it i) : it(i) {}
    reference operator*() const         { return *(*it); }
    pointer operator->() const          { return (*it); }
    bool operator==(const_iterator y)   { return it == y.it; }
    bool operator!=(const_iterator y)   { return it != y.it; }
    const_iterator operator++()         { ++it; return *this; }
  };

  struct TrueFunc {
    template <typename Src_Type>
    bool operator()(const Src_Type&) { return true; }
  };
  struct IdAddr {
    const Dst_Type* operator()(const Dst_Type& x) { return &x; }
  };

  template <typename Iter,
            typename Src_Type = typename std::iterator_traits<Iter>::value_type,
            typename CondFunc = std::function< bool (const Src_Type&)>,
            typename Reinterpreter = std::function< const Dst_Type* (const Src_Type&)>>
  void
  append(Iter first, Iter last, CondFunc cond,
         Reinterpreter get_ptr = IdAddr()) {
    for (auto i = first; i != last; ++i) {
      if (cond(*i))
        ptrs.push_back(get_ptr(*i));
    }
  }

  bool empty() const                    { return ptrs.empty(); }
  dim_type size() const                 { return ptrs.size(); }
  const_iterator begin() const          { return ptrs.begin(); }
  const_iterator end() const            { return ptrs.end(); }
  void push_back(const value_type* p)   { ptrs.push_back(p); }

private:
  List ptrs;
};

} // namespace pplite

#endif // !defined(pplite_utils_hh)
