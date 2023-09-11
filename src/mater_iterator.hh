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

#ifndef pplite_mater_iterator_hh
#define pplite_mater_iterator_hh 1

#include "globals.hh"
#include "Con.hh"
#include "Gen.hh"
#include "support_utils.hh"

#include <cassert>
#include <iterator>
#include <memory>
#include <type_traits>

namespace pplite {

inline bool
skippable(const Con* c) {
  return c->is_strict_inequality() && c->is_tautological();
}
inline bool
skippable(const Gen*) { return false; }

template <typename> class Proxy_Iter;

template <typename Sys, typename Impl>
struct Mater_Sys {
  // Types.
  using sys_type = Sys;
  using container_type = typename sys_type::rows_type;
  using size_type = typename container_type::size_type;
  using value_type = typename container_type::value_type;
  using cache_type = std::vector<std::unique_ptr<value_type>>;
  using const_iterator = Proxy_Iter<Mater_Sys>;

  // Data members.
  const sys_type* ptr1;
  const sys_type* ptr2;
  dim_type sk_offset;
  dim_type ns_offset;
  dim_type end_offset;
  mutable cache_type cache;
  const Impl* impl;

  // Constructor.
  explicit Mater_Sys(const Impl* p_impl)
    : impl(p_impl) {
    std::tie(ptr1, ptr2) = impl->get_sys_ptrs(ptr1);
    sk_offset = num_rows(ptr1->sing_rows) + num_rows(ptr2->sing_rows);
    ns_offset = sk_offset + num_rows(ptr1->sk_rows) + num_rows(ptr2->sk_rows);
    end_offset = ns_offset + num_rows(ptr1->ns_rows) + num_rows(ptr2->ns_rows);
    cache.resize(end_offset - ns_offset);
  }

  // Movable.
  Mater_Sys(Mater_Sys&&) noexcept = default;
  Mater_Sys& operator=(Mater_Sys&&) noexcept = default;
  ~Mater_Sys() = default;
  // Non-copiable.
  Mater_Sys() = delete;
  Mater_Sys(const Mater_Sys&) = delete;
  Mater_Sys& operator=(const Mater_Sys&) = delete;

  dim_type end_pos() const { return end_offset; }

  bool is_skippable(dim_type pos) const {
    assert(pos < end_pos());
    // Generators are never skippable.
    if (std::is_same<value_type, Gen>::value)
      return false;
    // Equations are not skippable.
    if (pos < sk_offset)
      return false;
    if (pos < ns_offset)
      return skippable(get_value_ptr(pos));
    if (impl->is_necessarily_closed())
      return true;
    // Non-skel element: we only skip the efc in ptr1
    // (the non-pending constraints).
    if (impl->sat_g.empty())
      return false;
    pos -= ns_offset;
    return pos < num_rows(ptr1->ns_rows) &&
      detail::is_empty_face_cutter(ptr1->ns_rows[pos],
                                   impl->sat_g, impl->gs.sk_rows);
  }

  void maybe_cache(dim_type cache_pos,
                   const sys_type* ptr,
                   dim_type ns_pos) const {
    if (cache[cache_pos] == nullptr) {
      auto m = detail::materialize(ptr->ns_rows[ns_pos], ptr->sk_rows);
      cache[cache_pos].reset(new value_type(std::move(m)));
    }
  }

  const value_type* get_value_ptr(dim_type pos) const {
    if (pos < sk_offset) {
      // Singular.
      if (pos < num_rows(ptr1->sing_rows))
        return &ptr1->sing_rows[pos];
      pos -= num_rows(ptr1->sing_rows);
      return &ptr2->sing_rows[pos];
    }
    if (pos < ns_offset) {
      // Skeleton.
      pos -= sk_offset;
      if (pos < num_rows(ptr1->sk_rows))
        return &ptr1->sk_rows[pos];
      pos -= num_rows(ptr1->sk_rows);
      return &ptr2->sk_rows[pos];
    }
    // Non-skeleton.
    assert(pos < end_pos());
    pos -= ns_offset;
    if (pos < num_rows(ptr1->ns_rows)) {
      maybe_cache(pos, ptr1, pos);
      return cache[pos].get();
    }
    maybe_cache(pos, ptr2, pos - num_rows(ptr1->ns_rows));
    return cache[pos].get();
  }

  const_iterator cbegin() const { return const_iterator(this, false); }
  const_iterator cend() const { return const_iterator(this, true); }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }

  dim_type size() const { return std::distance(cbegin(), cend()); }

}; // Mater_Sys

template <typename Cont>
struct Cont_Proxy {
  using container_type = Cont;
  using value_type = typename container_type::value_type;
  using const_iterator = Proxy_Iter<Cont_Proxy>;

  container_type c;

  explicit Cont_Proxy(const container_type& cont) : c(cont) {}
  explicit Cont_Proxy(container_type&& cont) : c(std::move(cont)) {}

  const value_type* get_value_ptr(dim_type pos) const {
    return &(c[pos]);
  }
  bool is_skippable(dim_type pos) const {
    return skippable(get_value_ptr(pos));
  }
  dim_type end_pos() const { return num_rows(c); }
  dim_type size() const { return std::distance(cbegin(), cend()); }

  const_iterator cbegin() const { return const_iterator(this, false); }
  const_iterator cend() const { return const_iterator(this, true); }
  const_iterator begin() const { return cbegin(); }
  const_iterator end() const { return cend(); }
}; // Cont_Proxy


// A bidirectional const_iterator based on Proxy.
template <typename Proxy>
class Proxy_Iter {
public:
  using proxy_type = Proxy;
  const proxy_type* ptr;
  dim_type pos;

  void skip_forward() {
    while (pos < ptr->end_pos() && ptr->is_skippable(pos))
      ++pos;
  }
  void skip_backward() noexcept {
    assert(pos > 0);
    --pos;
    while (pos > 0 && ptr->is_skippable(pos))
      --pos;
  }

  using iterator_category = std::bidirectional_iterator_tag;
  using value_type = typename proxy_type::value_type;
  using pointer = const value_type*;
  using reference = const value_type&;
  using difference_type = dim_type;

  Proxy_Iter() noexcept : ptr(nullptr), pos(0) {}
  Proxy_Iter(const proxy_type* msys, bool at_end)
    : ptr(msys), pos(0) {
    assert(ptr != nullptr);
    if (at_end)
      pos = ptr->end_pos();
    else
      skip_forward();
  }

  Proxy_Iter(const Proxy_Iter&) = default;
  Proxy_Iter(Proxy_Iter&&) noexcept = default;
  Proxy_Iter& operator=(const Proxy_Iter&) = default;
  Proxy_Iter& operator=(Proxy_Iter&&) noexcept = default;
  ~Proxy_Iter() = default;

  pointer operator->() const { return ptr->get_value_ptr(pos); }
  reference operator*() const { return *(operator->()); }

  Proxy_Iter& operator++() {
    assert(pos != ptr->end_pos());
    ++pos;
    skip_forward();
    return *this;
  }
  Proxy_Iter& operator--() {
    assert(pos > 0);
    --pos;
    skip_backward();
    return *this;
  }

  Proxy_Iter operator++(int) {
    auto old = *this;
    operator++();
    return old;
  }
  Proxy_Iter operator--(int) {
    auto old = *this;
    operator--();
    return old;
  }

}; // class Proxy_Iter

template <typename Proxy>
inline bool
operator==(const Proxy_Iter<Proxy>& x, const Proxy_Iter<Proxy>& y) {
  assert(x.ptr == y.ptr);
  return x.pos == y.pos;
}

template <typename Proxy>
inline bool
operator!=(const Proxy_Iter<Proxy>& x, const Proxy_Iter<Proxy>& y) {
  return !(x == y);
}

} // namespace pplite

#endif // !defined(pplite_mater_iterator_hh)
