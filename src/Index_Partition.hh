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

#ifndef pplite_Index_Partition_hh
#define pplite_Index_Partition_hh 1

#include <algorithm>
#include <iterator>
#include <vector>

namespace pplite {
namespace detail {

enum class Range { EQ, POS, NEG, POS_NEG };

using Ranges = std::vector<Range>;

template <Range policy>
class Index_Partition {
  const Ranges& labels;
public:
  Index_Partition() = delete;
  Index_Partition(const Index_Partition&) = default;
  Index_Partition& operator=(const Index_Partition&) = default;
  ~Index_Partition() = default;

  // Not movable
  Index_Partition(Index_Partition&&) = delete;
  Index_Partition& operator=(Index_Partition&&) = delete;

  explicit Index_Partition(const Ranges& ls) : labels(ls) {}

  dim_type tot_size() const { return labels.size(); }
  Range operator[](dim_type i) const { return labels[i]; }

  // read-only forward iterator for `policy' labels
  template <Range policy_it>
  class iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = dim_type;
    using pointer = const value_type*;
    using reference = const value_type&;
    using difference_type = std::ptrdiff_t;

  private:
    const Index_Partition& part;
    dim_type pos;
    dim_type end_pos;

  public:
    iterator(const Index_Partition& partition, bool get_first)
      : part(partition), end_pos(part.tot_size()) {
      if (!get_first){
        pos = end_pos;
        return;
      }
      // get_first
      for (pos = 0; pos < end_pos; ++pos){
        if (part[pos] == policy_it)
          break;
      }
    }

    value_type operator*() const { return pos; }
    iterator operator++() {
      if (pos == end_pos) return *this;
      do { ++pos; }
      while (pos != end_pos && part[pos] != policy_it);
      return *this;
    }
    bool operator==(const iterator& i) const { return pos == i.pos; }
    bool operator!=(const iterator& i) const { return pos != i.pos; }
  }; // iterator

  iterator<policy> begin() const {
    return iterator<policy>(*this, true);
  }
  iterator<policy> end() const {
    return iterator<policy>(*this, false);
  }
  bool empty() const {
    return std::find(labels.begin(), labels.end(), policy) == labels.end();
  }
}; // class Index_Partition

} // namespace detail
} // namespace pplite

#endif // not defined pplite_Index_Partition_hh
