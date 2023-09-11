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

#include "pplite-config.h"

#include "F_Poly.hh"
#include "Poly_widen.hh"

#include <algorithm>
#include <map>
#include <memory>

namespace pplite {

namespace detail {

using Block = F_Poly_Impl::Block;
using Blocks = F_Poly_Impl::Blocks;
using Factor = F_Poly_Impl::Factor;
using Factors = F_Poly_Impl::Factors;

template <typename T>
inline void
splice(std::vector<T>& xs, std::vector<T>& ys) {
  xs.reserve(xs.size() + ys.size());
  for (auto& y : ys)
    xs.push_back(std::move(y));
  ys.clear();
}

template <typename T, typename U, typename Pred>
void
synchronous_sort(std::vector<T>& xs, std::vector<U>& ys,
                 Pred x_lesser) {
  assert(num_rows(xs) == num_rows(ys));

  // Compute the order of the elements according to xs elems.
  auto pred = [&xs, &x_lesser] (dim_type i, dim_type j)
              { return x_lesser(xs[i], xs[j]); };
  std::list<dim_type> order;
  for (auto i : bwd_index_range(xs))
    order.insert(std::lower_bound(order.begin(), order.end(), i, pred), i);

  assert(num_rows(xs) == num_rows(order));

  // Sort xs and ys according to the computed order.
  std::vector<T> xs_res;
  xs_res.reserve(xs.size());
  std::vector<U> ys_res;
  ys_res.reserve(ys.size());
  for (auto i : order) {
    xs_res.push_back(std::move(xs[i]));
    ys_res.push_back(std::move(ys[i]));
  }

  xs = std::move(xs_res);
  ys = std::move(ys_res);
}

template <typename T, typename U>
inline void
synchronous_sort(std::vector<T>& xs, std::vector<U>& ys) {
  assert(num_rows(xs) == num_rows(ys));
  synchronous_sort(xs, ys, std::less<T>());
}

inline void
add_var(Block& b, Var v) {
  if (std::find(b.begin(), b.end(), v.id()) == b.end())
    b.push_back(v.id());
}

// Sort according to the positions in ext_b,
// so that a convert_back yields a sorted block.
inline void
sort_back(Block& int_b, const Block& ext_b) {
  sort(int_b.begin(), int_b.end(),
       [&ext_b] (dim_type i, dim_type j) { return (ext_b[i] < ext_b[j]); });
}

inline void
sort_back(Blocks& int_bs, const Block& ext_b) {
  for (auto& int_b : int_bs)
    sort_back(int_b, ext_b);
}

inline bool
contains(const Block& b, dim_type i) {
  return std::find(b.begin(), b.end(), i) != b.end();
}

inline bool
contains(const Block& b1, const Block& b2) {
  return all_of(b2, [&b1](dim_type i) { return contains(b1, i); });
}

inline Var
convert(Var v, const Block& b) {
  for (auto i : bwd_index_range(b)) {
    if (b[i] == v.id())
      return Var(i);
  }
  assert(false);
  abort();
}

inline Linear_Expr
convert(const Linear_Expr& expr, const Block& b) {
  Linear_Expr res;
  for (auto i : bwd_index_range(b))
    add_mul_assign(res, expr.get(b[i]), Var(i));
  return res;
}

inline Con
convert(const Con& c, const Block& b) {
  return Con(convert(c.linear_expr(), b), c.inhomo_term(), c.type());
}

inline Gen
convert(const Gen& g, const Block& b) {
  return Gen(g.type(), convert(g.linear_expr(), b), g.impl().inhomo);
}

inline Block
convert_back(const Block& int_b, const Block& ext_b) {
  auto res = int_b;
  for (auto& i : res)
    i = ext_b[i];
  return res;
}

inline Linear_Expr
convert_back(const Linear_Expr& expr, const Block& b) {
  Linear_Expr res;
  for (auto i : bwd_index_range(b))
    add_mul_assign(res, expr.get(i), Var(b[i]));
  return res;
}

inline Con
convert_back(const Con& c, const Block& b) {
  return Con(convert_back(c.linear_expr(), b), c.inhomo_term(), c.type());
}

inline Gen
convert_back(const Gen& g, const Block& b) {
  return Gen(g.type(), convert_back(g.linear_expr(), b), g.impl().inhomo);
}

inline Linear_Expr
convert_back_and_forth(const Linear_Expr& expr,
                       const Block& b_back, const Block& b_forth) {
  assert(contains(b_forth, b_back));
  Linear_Expr res;
  for (auto i : bwd_index_range(b_back) ) {
    if (expr.get(i).is_zero())
      continue;
    Var v_forth = convert(Var(b_back[i]), b_forth);
    add_mul_assign(res, expr.get(i), v_forth);
  }
  return res;
}

inline Con
convert_back_and_forth(const Con& c,
                       const Block& b_back, const Block& b_forth) {
  return Con(convert_back_and_forth(c.linear_expr(), b_back, b_forth),
             c.inhomo_term(), c.type());
}

inline dim_type
has_block_index(const Blocks& bs, Var v) {
  for (auto i : bwd_index_range(bs) ) {
    if (contains(bs[i], v.id()))
      return i;
  }
  return not_a_dim();
}

inline dim_type
find_block_index(const Blocks& bs, Var v) {
  auto res = has_block_index(bs, v);
  assert(res != not_a_dim());
  return res;
}

inline dim_type
find_block_index(const Blocks& bs, const Linear_Expr& expr) {
  // Note: assumes that bs contains a block b containing all dims
  // having a non-zero coeff in expr.
  auto i = expr.last_nonzero();
  assert(i < expr.space_dim());
  return find_block_index(bs, Var(i));
}

inline dim_type
find_block_index(const Blocks& bs, const Con& c) {
  return find_block_index(bs, c.linear_expr());
}

inline dim_type
get_containing_block_index(const Blocks& bs, const Block& b) {
  for (auto i : bwd_index_range(bs)) {
    if (!are_disjoint(bs[i], b)) {
      assert(contains(bs[i], b));
      return i;
    }
  }
  assert(false);
  abort();
}

inline dim_type
has_containing_block_index(const Blocks& bs, const Block& b) {
  for (auto i : bwd_index_range(bs)) {
    if (are_disjoint(bs[i], b))
      continue;
    return contains(bs[i], b) ? i : not_a_dim();
  }
  return not_a_dim();
}

inline void
reduce_blocks(Blocks& bs, Var v) {
  for (auto& b : bs) {
    for (auto& i : b)
      if (i > v.id())
        --i;
  }
}

inline void
reduce_blocks(Blocks& bs, const Index_Set& vars) {
  // Iterating backward is meant
  for (auto it = vars.end(), it_end = vars.begin(); it != it_end; ) {
    --it;
    detail::reduce_blocks(bs, Var(*it));
  }
}

inline void
concat(Block& b1, const Block& b2) {
  b1.insert(b1.end(), b2.begin(), b2.end());
}

inline void
concat(Block& b1, Factor& f1, const Block& b2, const Factor& f2) {
  concat(b1, b2);
  f1.concatenate_assign(f2);
}

inline Factor
itv_to_factor(const Itv& itv, Topol topol) {
  assert(not itv.is_empty());
  Factor f(1, topol);
  Var var(0);
  if (itv.is_singleton())
    f.add_con(get_eq_con(var, itv));
  else {
    if (itv.has_lb())
      f.add_con(get_lb_con(var, itv));
    if (itv.has_ub())
      f.add_con(get_ub_con(var, itv));
  }
  return f;
}

inline void
concat(Block& b, Factor& f, dim_type i, const Itv& itv) {
  assert(not itv.is_empty());
  b.push_back(i);
  Factor f_itv = itv_to_factor(itv, f.topology());
  f.concatenate_assign(f_itv);
}

inline Block
get_proper_itvs_block(const Itvs& itvs) {
  Block res;
  for (auto i : index_range(itvs)) {
    if (not itvs[i].is_empty())
      res.push_back(i);
  }
  return res;
}

struct Join_Info {
  Block itvs;
  Index_Set blocks;
};

Join_Info
joining_itvs_and_blocks(const Itvs& itvs, const Blocks& bs,
                        const Block& joiner) {
  Join_Info j_info;
  for (auto i : joiner) {
    if (not itvs[i].is_empty())
      j_info.itvs.push_back(i);
  }
  for (auto i : bwd_index_range(bs)) {
    if (!are_disjoint(bs[i], joiner))
      j_info.blocks.set(i);
  }
  return j_info;
}

Index_Set
joining_blocks(const Blocks& bs, const Block& joiner) {
  Index_Set joining;
  for (auto i : bwd_index_range(bs)) {
    if (!are_disjoint(bs[i], joiner))
      joining.set(i);
  }
  return joining;
}

void
join_blocks(Blocks& bs, const Block& joiner) {
  auto joining = joining_blocks(bs, joiner);
  assert(joining.size() > 0);
  if (joining.size() == 1)
    return;
  auto join_idx = joining.first();
  joining.reset(join_idx);
  for (auto i : joining)
    concat(bs[join_idx], bs[i]);
  erase_using_sorted_indices(bs, joining);
}

std::pair<Block, Factor>
join_factors(const Blocks& bs, const Factors& fs,
             const Index_Set& joining) {
  auto i = joining.first();
  auto res = std::make_pair(bs[i], fs[i]);
  for (auto j : joining) {
    if (i == j) continue;
    concat(res.first, bs[j]);
    res.second.concatenate_assign(fs[j]);
  }
  return res;
}

void
remap_factor(Factor& f, const Block& src_b, const Block& dst_b) {
  const auto sz = num_rows(src_b);
  assert(sz == num_rows(dst_b));
  Dims pf(sz);
  for (auto j : range(sz)) {
    for (auto k : range(sz)) {
      if (src_b[j] == dst_b[k]) {
        pf[j] = k;
        break;
      }
    }
  }
  f.map_space_dims(pf);
}

std::pair<Block, Factor>
merge(Topol topol, const Itvs& itvs, const Blocks& bs, const Factors& fs,
      const Join_Info& j_info) {
  auto res = std::make_pair(Block(), Factor(0, topol));
  auto& res_b = res.first;
  auto& res_f = res.second;

  if (not j_info.blocks.empty()) {
    auto first = j_info.blocks.first();
    for (auto i : j_info.blocks) {
      if (i != first)
        concat(res_b, res_f, bs[i], fs[i]);
      else {
        res_b = bs[i];
        res_f = fs[i];
      }
    }
  }
  for (auto i : j_info.itvs)
    concat(res_b, res_f, i, itvs[i]);
  return res;
}

inline Factor
merge_and_remap(Topol topol, const Itvs& itvs,
                const Blocks& bs, const Factors& fs,
                const Join_Info& j_info, const Block& dst_b) {
  auto p = merge(topol, itvs, bs, fs, j_info);
  auto& src_b = p.first;
  auto& src_f = p.second;
  if (src_b != dst_b)
    remap_factor(src_f, src_b, dst_b);
  return src_f;
}

template <typename Cons_Proxy>
Blocks
cons_based_partition(dim_type sdim, const Cons_Proxy& cs) {
  assert(sdim > 1);
  Blocks res(sdim);
  for (auto i : range(sdim))
    res[i].push_back(i);
  for (const auto& c : cs) {
    // check for (syntactic) inconsistency
    if (c.is_inconsistent())
      return Blocks();
    Block b_c = extract_block(c);
    // no join needed if c is an interval constraint
    if (num_rows(b_c) > 1) {
      join_blocks(res, b_c);
      if (num_rows(res) == 1)
        return res;
    }
  }
  return res;
}

Factors
project_on_partition(const Factor& f, const Blocks& bs) {
  const auto nbs = num_rows(bs);
  Factors res;
  res.reserve(nbs);
  // init with universe of correct dim and topology
  for (const auto& b : bs)
    res.emplace_back(num_rows(b), f.topology());
  // partition (and convert) the constraints
  for (const auto& c : f.cons()) {
    auto i = find_block_index(bs, c);
    assert(0 <= i && i < nbs);
    res[i].add_con(convert(c, bs[i]));
  }
  // minimize // FIXME: needed? worth?
  for (auto& ri : res)
    ri.minimize();
  return res;
}

// Sorts dims in (external) ascending order and computes normalized form.
inline void
normalize_factor(Block& b, Factor& f) {
  assert(num_rows(b) == f.space_dim());
  assert(f.space_dim() > 0);
  if (not std::is_sorted(b.begin(), b.end())) {
    auto nb = num_rows(b);
    Dims perm(nb), inv_perm(nb);
    std::iota(inv_perm.begin(), inv_perm.end(), 0);
    synchronous_sort(b, inv_perm);
    for (auto i : range(nb))
      perm[inv_perm[i]] = i;
    f.map_space_dims(perm);
  }
  (void) f.normalized_cons();
}

// Factorizes non-trivial block b/f:
// this puts the factorization in new_bs/new_fs (which are empty on entry);
// if f is detected inconsistent, new_bs and new_fs are left empty.
// Note: b and f are *consumed*, since they *might* be moved into
// the output parameters new_bs/new_fs.
// Note: we do not force minimization of f, so as to be lazy if possible.
inline void
factorize(Block& b, Factor& f, bool normalize,
          Blocks& new_bs, Factors& new_fs) {
  assert(num_rows(b) == f.space_dim());
  assert(f.space_dim() > 1);
  assert(new_bs.empty() && new_fs.empty());

  struct Counter {
    long calls = 0;
    long empty = 0;
    long non_fact = 0;
    long fact = 0;
    long blocks_out = 0;

    void check(const Blocks& new_bs) {
      ++calls;
      if (new_bs.size() == 0)
        ++empty;
      else if (new_bs.size() == 1) {
        ++non_fact;
        ++blocks_out;
      } else {
        ++fact;
        blocks_out += new_bs.size();
        std::vector<dim_type> vs;
        for (const auto& b : new_bs)
          vs.push_back(b.size());
        sort(vs.begin(), vs.end());
        auto sdim = std::accumulate(vs.begin(), vs.end(), 0);
        std::cerr << "=== fact: " << sdim
                  << " -> ";
        for (auto s : vs) { std::cerr << s << " , "; }
        std::cerr << "\n";
      }
    }

    ~Counter() {
      std::cerr << "\n=== detail::factorize counter:\n";
      std::cerr << "=== num_calls = " << calls << "\n";
      std::cerr << "=== found empty = " << empty << "\n";
      std::cerr << "=== non-factorizable = " << non_fact << "\n";
      std::cerr << "=== factorizable = " << fact << "\n";
      std::cerr << "=== blocks out = " << blocks_out << "\n";
      std::cerr << "=== avg blocks out = "
                << ((double) blocks_out) / calls << "\n";
    }
  };

#define FACT_STATS 0
#if FACT_STATS
  static Counter cc;
  #define UPDATE_STATS do { cc.check(new_bs); } while (false)
#else
  #define UPDATE_STATS
#endif

  // try to be lazy: partition on syntactic constraints
  Blocks part = cons_based_partition(f.space_dim(), f.cons());
  auto part_size = num_rows(part);
  if (part_size == 0) {
    UPDATE_STATS;
    return;
  }
  if (part_size == 1) {
    if (normalize) {
      normalize_factor(b, f);
      if (f.is_empty()) {
        UPDATE_STATS;
        return;
      }
    }
    new_bs.push_back(std::move(b));
    new_fs.push_back(std::move(f));
    UPDATE_STATS;
    return;
  }

  assert(part_size > 1);
  // force ascending (external) order in new blocks.
  sort_back(part, b);
  // actually compute partition projections
  new_fs = project_on_partition(f, part);
  // convert from internal to external new blocks
  new_bs = std::move(part);
  for (auto i : range(part_size)) {
    auto& nb = new_bs[i];
    auto& nf = new_fs[i];
    nb = convert_back(nb, b);
    if (normalize) {
      normalize_factor(nb, nf);
      if (nf.is_empty()) {
        new_bs.clear();
        new_fs.clear();
        UPDATE_STATS;
        return;
      }
    }
  }
  UPDATE_STATS;
}

// Factorizes block having index idx; returns true if it is empty.
inline bool
factorize(dim_type idx, Blocks& bs, Factors& fs, bool normalize) {
  auto& b = bs[idx];
  auto& f = fs[idx];
  const auto sdim = f.space_dim();
  if (sdim == 1) {
    f.minimize();
    return f.is_empty();
  }
  assert(sdim > 1);
  const bool f_was_minimized = f.is_minimized();
  Blocks new_bs;
  Factors new_fs;
  factorize(b, f, normalize, new_bs, new_fs);
  if (new_bs.empty())
    return true;
  if (f_was_minimized) {
    b = std::move(new_bs.back());
    f = std::move(new_fs.back());
    new_bs.pop_back();
    new_fs.pop_back();
    splice(bs, new_bs);
    splice(fs, new_fs);
    return false;
  }

  assert(not f_was_minimized);
  // Blocks have been factorized according to syntactic constraints,
  // hence some of them might be empty or further factorizable.
  // Now really force minimization and factorize them.
  Blocks tmp_bs, tmp_bs_i;
  Factors tmp_fs, tmp_fs_i;
  for (auto i : index_range(new_fs)) {
    auto& bi = new_bs[i];
    auto& fi = new_fs[i];
    fi.minimize();
    if (fi.is_empty())
      return true;
    if (fi.space_dim() == 1) {
      tmp_bs.push_back(std::move(bi));
      tmp_fs.push_back(std::move(fi));
    } else {
      factorize(bi, fi, normalize, tmp_bs_i, tmp_fs_i);
      assert(not tmp_bs_i.empty());
      // accumulate partial results
      splice(tmp_bs, tmp_bs_i);
      splice(tmp_fs, tmp_fs_i);
    }
  }
  // Now we have maximal factorization (of b/f) in tmp_bs/tmp_fs.
  // Move last element in b/f and splice the remaining elements in bs/fs.
  b = std::move(tmp_bs.back());
  f = std::move(tmp_fs.back());
  tmp_bs.pop_back();
  tmp_fs.pop_back();
  splice(bs, tmp_bs);
  splice(fs, tmp_fs);
  return false;
}

inline Block
block_union(const Block& b1, const Block& b2) {
  Block res = b1;
  for (auto i2 : b2) {
    if (!contains(b1, i2))
      res.push_back(i2);
  }
  return res;
}

Blocks
blocks_lub(const Itvs& itvs1, const Blocks& bs1,
           const Itvs& itvs2, const Blocks& bs2) {
  auto add_itvs_blocks
    = [](const Blocks& bs1, const Itvs& itvs2, Blocks& lub) {
        for (const auto& b1 : bs1) {
          for (auto i1 : b1) {
            if (itvs2[i1].is_empty())
              continue;
          }
          lub.push_back(b1);
        }
      };

  Blocks lub;
  add_itvs_blocks(bs1, itvs2, lub);
  add_itvs_blocks(bs2, itvs1, lub);

  for (const auto& b1 : bs1) {
    bool first = true;
    for (const auto& b2 : bs2) {
      if (are_disjoint(b1, b2))
        continue;
      if (first) {
        lub.push_back(block_union(b1, b2));
        first = false;
      } else
        lub.back() = block_union(lub.back(), b2);
    }
  }

  Index_Set to_be_rem;
  const auto nb = num_rows(lub);
  for (auto i : range(nb)) {
    if (to_be_rem.test(i))
      continue;
    auto j1 = i;
    for (auto j2 : range(j1+1, nb)) {
      if (to_be_rem.test(j2))
        continue;
      if (are_disjoint(lub[j1], lub[j2]))
        continue;
      lub[j2] = block_union(lub[j1], lub[j2]);
      to_be_rem.set(j1);
      j1 = j2;
    }
  }
  erase_using_sorted_indices(lub, to_be_rem);
  return lub;
}

struct Refactor_Proxy {
  const Topol topol;
  const Itvs& itvs;
  const Blocks& bs;
  const Factors& fs;
  const Blocks& lub;
  mutable std::vector<std::unique_ptr<Factor>> cache;

  const Itv& get_itv(dim_type i) const {
    assert(0 <= i && i <= num_rows(itvs) && !itvs[i].is_empty());
    return itvs[i];
  }

  const Factor& operator[](dim_type i) const {
    assert(0 <= i && i <= num_rows(lub));
    if (cache[i])
      return *cache[i];
    auto j_info = joining_itvs_and_blocks(itvs, bs, lub[i]);
    if (j_info.blocks.size() == 1
        && bs[j_info.blocks.first()] == lub[i]) {
      assert(j_info.itvs.empty());
      // No need to join&remap
      auto idx = j_info.blocks.first();
      return fs[idx];
    } else {
      auto f = merge_and_remap(topol, itvs, bs, fs, j_info, lub[i]);
      cache[i].reset(new Factor(std::move(f)));
      return *cache[i];
    }
  }

  Refactor_Proxy(const F_Poly::Impl& impl, const Blocks& lub_)
    : topol(impl.topol), itvs(impl.itvs),
      bs(impl.blocks), fs(impl.factors),
      lub(lub_), cache(num_rows(lub_)) {
    assert(not impl.empty);
  }
  ~Refactor_Proxy() = default;

  Refactor_Proxy(const Refactor_Proxy&) = delete;
  Refactor_Proxy& operator=(const Refactor_Proxy&) = delete;
  Refactor_Proxy(Refactor_Proxy&&) = delete;
  Refactor_Proxy& operator=(Refactor_Proxy&&) = delete;
};

} // namespace detail

void
F_Poly::erase_empty_block(dim_type block_idx) {
  assert(0 <= block_idx && block_idx < num_rows(blocks));
  assert(blocks[block_idx].empty());
  blocks.erase(blocks.begin() + block_idx);
  factors.erase(factors.begin() + block_idx);
}

void
F_Poly::factorize(bool normalize) const {
  if (empty || blocks.empty())
    return;

  // Process blocks in increasing size order (i.e., factor space dim).
  Block proc_order(num_rows(blocks));
  std::iota(proc_order.begin(), proc_order.end(), 0);
  std::sort(proc_order.begin(), proc_order.end(),
            [this](dim_type i, dim_type j) {
              return blocks[i].size() < blocks[j].size();
            });

  // logically const
  auto& x = const_cast<F_Poly&>(*this);
  for (auto idx : proc_order) {
    if (detail::factorize(idx, x.blocks, x.factors, normalize)) {
      x.set_empty();
      return;
    }
  }
  x.blocks_to_itvs();
  assert(x.check_inv());
}

void
F_Poly::block_to_itv(dim_type block_idx) {
  assert(block_idx < num_rows(blocks));
  auto& b = blocks[block_idx];
  auto& f = factors[block_idx];
  assert(is_boxable(f));
  auto i = b.back();
  assert(is_block_dim(i));
  itvs[i] = f.get_bounds(Var(0));
  b.clear();
  erase_empty_block(block_idx);
}

void
F_Poly::blocks_to_itvs() {
  assert(!empty);
  for (auto i : bwd_index_range(blocks)) {
    if (is_boxable(factors[i]))
      block_to_itv(i);
  }
  assert(check_inv());
}

void
F_Poly::itv_to_block(dim_type i, dim_type block_idx) {
  assert(!empty);
  assert(is_itv_dim(i) && block_idx < num_rows(blocks));
  detail::concat(blocks[block_idx], factors[block_idx], i, itvs[i]);
  itvs[i].set_empty();
  assert(check_inv());
}

F_Poly
F_Poly::from_poly(const Poly& ph) {
  auto res = F_Poly(ph.space_dim(), ph.topology());
  const auto& cs = ph.cons();
  res.add_cons(cs.begin(), cs.end());
  return res;
}

Poly
F_Poly::to_poly() const {
  if (space_dim() == 0 || is_empty())
    return Poly(space_dim(),
                (is_empty()) ? Spec_Elem::EMPTY : Spec_Elem::UNIVERSE,
                topology());

  using namespace detail;
  Join_Info j_info;
  j_info.itvs = get_proper_itvs_block(itvs);
  j_info.blocks.set_until(num_rows(blocks));
  Block b(space_dim());
  std::iota(b.begin(), b.end(), 0);
  return merge_and_remap(topol, itvs, blocks, factors, j_info, b);
}

void
F_Poly::m_swap(F_Poly& y) noexcept {
  using std::swap;
  swap(dim, y.dim);
  swap(topol, y.topol);
  swap(empty, y.empty);
  swap(is_normalized, y.is_normalized);
  swap(itvs, y.itvs);
  swap(blocks, y.blocks);
  swap(factors, y.factors);
}

dim_type
F_Poly::merge(const Block& b) {
  assert(b.size() > 0);

  // Check if no merge is needed
  auto block_idx = detail::has_containing_block_index(blocks, b);
  if (block_idx != not_a_dim())
    return block_idx;

  // CHECKME: too preemptive?
  is_normalized = false;

  // Check if merging only itv dims
  if (are_itv_dims(b)) {
    block_idx = num_rows(blocks);
    blocks.push_back(Block());
    factors.push_back(Factor(0, topol));
    for (auto i : b)
      itv_to_block(i, block_idx);
    assert(check_inv());
    return block_idx;
  }

  const auto nb = num_rows(blocks);
  dim_type merge_idx = not_a_dim();
  Blocks merge_bs;
  merge_bs.reserve(nb);
  Factors merge_fs;
  merge_fs.reserve(nb);
  for (auto idx : range(nb)) {
    auto& b_idx = blocks[idx];
    auto& f_idx = factors[idx];
    if (detail::are_disjoint(b_idx, b)) {
      merge_bs.push_back(std::move(b_idx));
      merge_fs.push_back(std::move(f_idx));
    } else if (merge_idx == not_a_dim()) {
      assert(idx == num_rows(merge_bs));
      merge_idx = idx;
      merge_bs.push_back(std::move(b_idx));
      merge_fs.push_back(std::move(f_idx));
    } else {
      detail::concat(merge_bs[merge_idx], merge_fs[merge_idx],
                     b_idx, f_idx);
    }
  }
  assert(merge_idx != not_a_dim());
  blocks = std::move(merge_bs);
  factors = std::move(merge_fs);
  // Now merge itv dims, if any.
  for (auto i : b) {
    if (is_itv_dim(i))
      itv_to_block(i, merge_idx);
  }
  return merge_idx;
}

void
F_Poly::sync(const Blocks& bs) const {
  auto& x = const_cast<F_Poly&>(*this);
  for (auto i : index_range(bs)) {
    auto j = x.merge(bs[i]);
    assert(j >= i);
    if (x.blocks[j] != bs[i]) {
      detail::remap_factor(x.factors[j], x.blocks[j], bs[i]);
      x.blocks[j] = bs[i];
    }
    if (i != j) {
      using std::swap;
      swap(x.blocks[i], x.blocks[j]);
      swap(x.factors[i], x.factors[j]);
    }
  }
  x.is_normalized = false;
  assert(x.check_inv());
}

bool
F_Poly::check_inv() const {
  std::string reason;
#ifdef NDEBUG
  auto maybe_dump = []() {};
#else // In debugging mode, be noisy.
  auto maybe_dump = [this, &reason]() {
    std::cerr << reason << std::endl;
    this->ascii_dump(std::cerr);
    std::cerr << reason << std::endl;
  };
#endif

  if (is_normalized) {

    // No factor is empty.
    if (any_of(factors, std::mem_fn(&Factor::is_empty))) {
      reason = "F_Poly broken: empty factor when normalized";
      maybe_dump();
      return false;
    }

    // If a block has space dim 1, then it must be NOT topologically closed
    for (auto i : index_range(factors)) {
      if (is_boxable(factors[i])) {
        reason = "F_Poly broken: space dim 1 block is topologically closed";
        maybe_dump();
        return false;
      }
    }

    auto copy = *this;
    copy.is_normalized = false;
    if (copy.hash() != hash()) {
      reason = "F_Poly broken: marked normalized when it is not.";
      maybe_dump();
      return false;
    }
  }

  if (dim < 0) {
    reason = "F_Poly broken: invalid space dimension";
    maybe_dump();
    return false;
  }

  // If it is marked empty, itvs, blocks and factors must be empty too.
  if (empty) {
    if (not (itvs.empty() && blocks.empty() && factors.empty())) {
      reason = "F_Poly broken: empty poly has non-empty components";
      maybe_dump();
      return false;
    }
    // No other check for an empty polyhedron.
    return true;
  }

  // Here the polyhedron is NOT marked empty

  // itvs size should equals dim
  if (num_rows(itvs) != dim) {
    reason = "F_Poly broken: info size != dim";
    maybe_dump();
    return false;
  }

  // itvs should satisfy their invariant
  if (not all_of(itvs, std::mem_fn(&Itv::check_inv))) {
    reason = "F_Poly broken: an element of itvs is broken";
    maybe_dump();
    return false;
  }

  // Blocks should not be empty.
  for (const auto& block : blocks) {
    if (block.empty()) {
      reason = "F_Poly broken: a block has space dim 0";
      maybe_dump();
      return false;
    }
  }

  // As many factors as blocks.
  if (factors.size() != blocks.size()) {
    reason = "F_Poly broken: #factors != #blocks";
    maybe_dump();
    return false;
  }

  // Each factor space dim should match its block.
  for (auto i : index_range(factors)) {
    if (factors[i].space_dim() != num_rows(blocks[i])) {
      reason = "F_Poly broken: factor vs block space dim mismatch";
      maybe_dump();
      return false;
    }
  }

  // Checking for dimensions
  Index_Set itv_dims;
  for (auto i : bwd_index_range(itvs)) {
    if (is_itv_dim(i))
      itv_dims.set(i);
  }

  // Each block dimension appears once in blocks.
  Index_Set block_dims;
  for (const auto& block : blocks) {
    for (const auto i : block) {
      if (i < 0 || i >= dim) {
        reason = "F_Poly broken: block contains an illegal space dim";
        maybe_dump();
        return false;
      }
      if (block_dims.test(i)) {
        reason = "F_Poly broken: repeated space dim in blocks";
        maybe_dump();
        return false;
      }
      block_dims.set(i);
    }
  }

  // The sum of itvs and blocks dims should equal space dim
  if (dim != itv_dims.size() + block_dims.size()) {
    reason = "F_Poly broken: itvs + blocks dimension mismatch";
    maybe_dump();
    return false;
  }

  // No dim should appear in both ivs and blocks dims
  if (itv_dims.intersects(block_dims)) {
    reason = "F_Poly broken: itv dim occurs in blocks";
    maybe_dump();
    return false;
  }

  // All tests passed.
  return true;
}

void
F_Poly::add_con(const Con& c) {
  if (empty)
    return;
  Block b = detail::extract_block(c);
  if (b.size() == 0) {
    if (c.is_tautological())
      return;
    else {
      assert(c.is_inconsistent());
      set_empty();
      return;
    }
  }
  if (num_rows(b) == 1 && are_itv_dims(b) && !c.is_strict_inequality()) {
    assert(is_proper_interval_con(c));
    auto c_dim = b.back();
    const auto& c_den = c.coeff(Var(c_dim));
    Rational b_bound { c.inhomo_term(), c_den };
    // Move inhomo on rhs
    neg_assign(b_bound);

    Itv b_itv;
    if (c.is_equality())
      b_itv.set_singleton(std::move(b_bound));
    else if (c_den > 0)
      b_itv.set_lb(std::move(b_bound));
    else
      b_itv.set_ub(std::move(b_bound));
    if (itvs[c_dim].glb_assign(b_itv))
      set_empty();
    assert(check_inv());
    return;
  }

  auto idx = merge(b);
  Con c_int = detail::convert(c, blocks[idx]);
  factors[idx].add_con(std::move(c_int));
  is_normalized = false;
  // FIXME, CHECKME: why being so eager?
  if (factors[idx].is_empty())
    set_empty();
  assert(check_inv());
}

void
F_Poly::add_line_or_ray(const Gen& g) {
  assert(!empty && g.is_line_or_ray());
  Block b = detail::extract_block(g);
  if (num_rows(b) == 1 && are_itv_dims(b)) {
    auto g_dim = b.back();
    auto& itv = itvs[g_dim];
    if (g.is_line())
      itv.set_universe();
    else {
      assert(g.is_ray() && sgn(g.coeff(Var(g_dim))) != 0);
      if (sgn(g.coeff(Var(g_dim))) == 1)
        itv.unset_ub();
      else
        itv.unset_lb();
    }
    return;
  }
  auto idx = merge(b);
  auto g_int = detail::convert(g, blocks[idx]);
  factors[idx].add_gen(std::move(g_int));
  assert(check_inv());
}

Poly_Con_Rel
F_Poly::relation_with(const Con& c) const {
  if (is_empty())
    return Poly_Con_Rel::saturates()
      && Poly_Con_Rel::is_included()
      && Poly_Con_Rel::is_disjoint();

  if (c.is_inconsistent())
    return (c.is_strict_inequality() && c.inhomo_term() == 0)
      ? Poly_Con_Rel::saturates() && Poly_Con_Rel::is_disjoint()
      : Poly_Con_Rel::is_disjoint();

  if (c.is_tautological())
    return (c.is_equality() || c.inhomo_term() == 0)
      ? Poly_Con_Rel::saturates() && Poly_Con_Rel::is_included()
      : Poly_Con_Rel::is_included();

  assert(space_dim() > 0);
  Block b = detail::extract_block(c);
  assert(b.size() > 0);
  auto& x = const_cast<F_Poly&>(*this);
  auto i = x.merge(b);
  Con c_int = detail::convert(c, blocks[i]);
  return factors[i].relation_with(c_int);
}

Poly_Gen_Rel
F_Poly::relation_with(const Gen& g) const {
  assert(space_dim() >= g.space_dim());
  if (is_empty())
    return Poly_Gen_Rel::nothing();
  if (space_dim() == 0)
    return Poly_Gen_Rel::subsumes();

  if (g.is_line() || g.is_ray()) {
    Block b_g = detail::extract_block(g);
    for (auto i : index_range(b_g)) {
      if (is_itv_dim(i)) {
        if (g.is_line() && !itvs[i].is_universe())
          return Poly_Gen_Rel::nothing();
        if (g.is_ray()) {
          if (sgn(g.coeff(Var(i))) > 0 && itvs[i].has_ub())
            return Poly_Gen_Rel::nothing();
          else if (sgn(g.coeff(Var(i))) < 0 && itvs[i].has_lb())
            return Poly_Gen_Rel::nothing();
        }
      }
    }
    for (auto i : index_range(blocks)) {
      const auto& b = blocks[i];
      if (detail::are_disjoint(b, b_g))
        continue;
      Gen g_int = detail::convert(g, b);
      const auto& f = factors[i];
      if (f.relation_with(g_int) == Poly_Gen_Rel::nothing())
        return Poly_Gen_Rel::nothing();
    }
    return Poly_Gen_Rel::subsumes();
  }

  assert(g.is_point() || g.is_closure_point());
  for (auto i : index_range(itvs)) {
    if (is_itv_dim(i) && !itvs[i].contains(g.coeff(Var(i)), g.divisor()))
      return Poly_Gen_Rel::nothing();
  }
  for (auto i : index_range(blocks)) {
    const auto& b = blocks[i];
    Gen g_int = detail::convert(g, b);
    const auto& f = factors[i];
    if (f.relation_with(g_int) == Poly_Gen_Rel::nothing())
      return Poly_Gen_Rel::nothing();
  }
  return Poly_Gen_Rel::subsumes();
}

bool
F_Poly::min(const Affine_Expr& ae, Rational& value,
            bool* included_ptr, Gen* g_ptr) const {
  if (is_empty())
    return false;

  if (space_dim() == 0) {
    value = Rational(ae.inhomo);
    if (included_ptr)
      *included_ptr = true;
    if (g_ptr)
      *g_ptr = point();
    return true;
  }

  // Maybe needed to update g_ptr.
  Linear_Expr expr_g;
  Integer div_g = 1;
  const Integer zero_i;
  Integer lcm, mult_g, mult_i;

  auto update_gen = [&](const Linear_Expr& expr_g_i,
                        const Integer& div_g_i) {
    assert(g_ptr);
    lcm_assign(lcm, div_g, div_g_i);
    exact_div_assign(mult_g, lcm, div_g);
    expr_g *= mult_g;
    exact_div_assign(mult_i, lcm, div_g_i);
    add_mul_assign(expr_g, mult_i, expr_g_i);
    div_g = std::move(lcm);
  };

  // Init value and included.
  value = Rational::zero();
  bool included = true;
  // Note: we keep and update the boolean `included'
  // and (maybe) copy it into the optional output parameter.

  // Go through itv dims
  for (auto i : dim_range(*this)) {
    if (is_block_dim(i))
      continue;
    const auto& itv_i = itvs[i];
    auto sgn_i = sgn(ae.expr.get(i));
    if (sgn_i != 0) {
      if (sgn_i > 0 && itv_i.inf_lb())
        return false;
      if (sgn_i < 0 && itv_i.inf_ub())
        return false;
      const auto& bound = (sgn_i > 0) ? itv_i.lb : itv_i.ub;
      add_mul_assign(value, Rational(ae.expr.get(i)), bound);
      if (g_ptr) {
        const auto& expr_g_i = bound.get_num() * Var(i);
        const auto& div_g_i = bound.get_den();
        update_gen(expr_g_i, div_g_i);
      }
    } else {
      assert(sgn_i == 0);
      // value is not affected: this mess in only needed to produce gen.
      if (g_ptr) {
        if (itv_i.is_universe()) {
          // Equivalent to choosing 0 as g_i coordinate
          continue;
        }
        const auto& bound = itv_i.has_lb() ? itv_i.lb : itv_i.ub;
        const auto& expr_g_i = bound.get_num() * Var(i);
        const auto& div_g_i = bound.get_den();
        update_gen(expr_g_i, div_g_i);
      }
    }
  }

  Rational value_i;
  bool included_i;
  // Note: included_i is needed even when included_ptr is false.
  bool* included_i_ptr = (included_ptr || g_ptr) ? &included_i : nullptr;
  Gen g_i = point();
  Gen* g_i_ptr = g_ptr ? &g_i : nullptr;

  // Go through block dims.
  for (auto i : index_range(blocks)) {
    const auto& b_i = blocks[i];
    auto expr_i = detail::convert(ae.expr, b_i);
    if (expr_i.is_zero())
      continue;
    auto ae_i = Affine_Expr(std::move(expr_i));
    if (!factors[i].min(ae_i, value_i, included_i_ptr, g_i_ptr))
      return false;
    value += value_i;
    if (included_i_ptr) {
      assert(included_i_ptr == &included_i);
      included = included && included_i;
    }
    if (g_ptr) {
      assert(g_i_ptr && (g_i_ptr == &g_i));
      auto expr_g_i = detail::convert_back(g_i.linear_expr(), b_i);
      const auto& div_g_i = g_i.divisor();
      update_gen(expr_g_i, div_g_i);
    }
  }

  value += Rational(ae.inhomo);
  if (included_ptr)
    *included_ptr = included;
  if (g_ptr)
    *g_ptr = included ? point(expr_g, div_g) : closure_point(expr_g, div_g);
  return true;
}

Itv
F_Poly::get_bounds(Var var) const {
  if (is_empty())
    return Itv(Spec_Elem::EMPTY);
  if (is_itv_dim(var.id()))
    return itvs[var.id()];
  else {
    auto idx = detail::find_block_index(blocks, var);
    Var b_var = detail::convert(var, blocks[idx]);
    return factors[idx].get_bounds(b_var);
  }
}

Itv
F_Poly::get_bounds(const Affine_Expr& ae) const {
  Itv res(Spec_Elem::EMPTY);
  if (is_empty())
    return res;

  res.set_singleton(Rational(ae.inhomo));

  Index_Set non_zeroes = ae.expr.non_zeroes();
  for (auto i : non_zeroes) {
    if (is_block_dim(i))
      continue;
    auto itv_i = itvs[i];
    itv_i.mul_assign(Rational(ae.expr[i]));
    res.add_assign(itv_i);
    if (res.is_universe())
      return res;
  }

  for (auto i : index_range(blocks)) {
    const auto& b_i = blocks[i];
    // check if block b_i is irrelevant wrt ae
    if (none_of(b_i, [&non_zeroes](dim_type d)
                     { return non_zeroes.test(d); }))
      continue;
    Linear_Expr expr_i = detail::convert(ae.expr, b_i);
    assert(not expr_i.is_zero());
    auto itv_i = factors[i].get_bounds(Affine_Expr(expr_i));
    res.add_assign(itv_i);
    if (res.is_universe())
      return res;
  }
  return res;
}

Itv
F_Poly::get_bounds(const Itv_Expr& ie) const {
  if (is_empty())
    return Itv(Spec_Elem::EMPTY);

  const auto& ie_vars = ie.first;
  const auto& ie_itvs = ie.second;
  assert(num_rows(ie_vars) == num_rows(ie_itvs));

  Itv res = Itv::zero();

  if (ie_vars.empty())
    return res;

  // Process intervals
  for (auto i : index_range(ie_vars)) {
    auto dim = ie_vars[i].id();
    const auto& ie_itv = ie_itvs[i];
    if (is_block_dim(dim) || ie_itv.is_zero())
      continue;
    auto prod_ub = itvs[dim];
    prod_ub.mul_assign(ie_itv.ub);
    auto prod_lb = itvs[dim];
    prod_lb.mul_assign(ie_itv.lb);
    prod_ub.lub_assign(prod_lb);
    res.add_assign(prod_ub);
    if (res.is_universe())
      return res;
  }

  // Process blocks
  for (auto j : index_range(blocks)) {
    const auto& bj = blocks[j];
    const auto& fj = factors[j];

    Itv_Expr bj_ie;
    auto& bj_vars = bj_ie.first;
    auto& bj_itvs = bj_ie.second;

    for (auto i : index_range(ie_vars)) {
      auto ie_var = ie_vars[i];
      const auto& ie_itv = ie_itvs[i];
      auto dim = ie_var.id();
      if (is_itv_dim(dim) || not detail::contains(bj, dim) || ie_itv.is_zero())
        continue;
      Var bj_var = detail::convert(ie_var, bj);
      bj_vars.push_back(bj_var);
      bj_itvs.push_back(ie_itv);
    }

    assert(num_rows(bj_vars) == num_rows(bj_itvs));
    if (bj_vars.empty())
      continue;
    auto bj_res = fj.get_bounds(bj_ie);
    res.add_assign(bj_res);
    if (res.is_universe())
      return res;
  }
  // Final (not universe) result
  return res;
}

Index_Set
F_Poly::get_unconstrained() const {
  const auto& x = *this;
  Index_Set res;
  if (x.is_empty())
    return res;
  // Process intervals
  for (auto i : bwd_dim_range(x)) {
    if (x.is_itv_dim(i) && x.itvs[i].is_universe())
      res.set(i);
  }
  if (x.is_normalized)
    return res;
  // Process blocks/factors
  for (auto i : index_range(x.blocks)) {
    const auto& fi = x.factors[i];
    auto fi_res = fi.get_unconstrained();
    const auto& bi = blocks[i];
    for (auto i : fi_res)
      res.set(bi[i]);
  }
  return res;
}

bool
F_Poly::is_bounded_expr(bool from_below, const Linear_Expr& expr) const {
  if (space_dim() == 0 || is_empty())
    return true;

  for (auto i : dim_range(expr)) {
    if (expr.get(i).is_zero() || is_block_dim(i))
      continue;
    const auto& itv = itvs[i];
    if (itv.is_universe())
      return false;
    bool positive = (expr.get(i) > 0);
    if (itv.inf_lb() && (positive == from_below))
      return false;
    if (itv.inf_ub() && (positive != from_below))
      return false;
  }

  for (auto i : index_range(blocks)) {
    const auto& b_i = blocks[i];
    auto expr_i = detail::convert(expr, b_i);
    if (!factors[i].is_bounded_expr(from_below, expr_i))
      return false;
  }
  return true;
}

Cons
F_Poly::itvs_to_cons() const {
  Cons res;
  for (auto i : index_range(itvs)) {
    if (is_block_dim(i))
      continue;
    Var var(i);
    const auto& itv = itvs[i];
    if (itv.is_singleton())
      res.push_back(get_eq_con(var, itv));
    else {
      if (itv.has_lb())
        res.push_back(get_lb_con(var, itv));
      if (itv.has_ub())
        res.push_back(get_ub_con(var, itv));
    }
  }
  return res;
}

Cons
F_Poly::copy_cons() const {
  if (empty)
    return { Con::zero_dim_false() };
  if (space_dim() == 0)
    return Cons();

  Cons res = itvs_to_cons();
  for (auto i : index_range(factors)) {
    const auto& bi = blocks[i];
    for (const auto& c : factors[i].cons())
      res.push_back(detail::convert_back(c, bi));
  }
  return res;
}

Gens
F_Poly::copy_gens() const {
  if (empty)
    return Gens();

  if (space_dim() == 0)
    return { point() };

  Block b;
  Factor f(0, topology());
  for (auto i : index_range(blocks))
    detail::concat(b, f, blocks[i], factors[i]);
  for (auto i : dim_range(*this))
    if (is_itv_dim(i))
      detail::concat(b, f, i, itvs[i]);
  assert(num_rows(b) == space_dim());
  assert(f.space_dim() == space_dim());
  f.map_space_dims(b);
  return f.copy_gens();
}

dim_type
F_Poly::num_min_cons() const {
  if (is_empty())
    return 1;
  dim_type res = 0;
  for (const auto& itv : proper_itvs())
    res += itv.num_min_cons();
  for (const auto& f : factors)
    res += f.num_min_cons();
  return res;
}

Gens_Info
F_Poly::gens_info() const {
  if (is_empty())
    return { 0, 0, 0, 0, 0 };

  dim_type ln = 0, r = 0, cp = 0, skp = 1, ns = 0;

  for (const auto& itv : proper_itvs()) {
    assert(not itv.is_empty());
    if (itv.is_universe())
      ++ln;
    else if (not itv.is_bounded())
      ++r;
    else if (not itv.is_singleton())
      skp *= 2;
  }

  dim_type f_ln, f_r, f_cp, f_skp, f_ns;
  for (const auto& f : factors) {
    std::tie(f_ln, f_r, f_cp, f_skp, f_ns) = f.gens_info();
    ln += f_ln;
    r += f_r;
    cp = cp * (f_cp + f_skp) + skp * f_cp;
    ns = ns * (f_ns + f_skp) + skp * f_ns;
    skp *= f_skp;
  }
  return { ln, r, cp, skp, ns };
}

dim_type
F_Poly::num_min_gens() const {
  if (is_empty())
    return 0;
  if (space_dim() == 0)
    return 1;
  dim_type ln, r, cp, skp, ns;
  std::tie(ln, r, cp, skp, ns) = gens_info();
  return ln + r + cp + skp + ns;
}

bool
F_Poly::equals(const F_Poly& y) const {
  const auto& x = *this;
  if (x.space_dim() != y.space_dim())
    return false;
  if (x.is_empty())
    return y.is_empty();
  if (y.is_empty())
    return false;
  if (x.space_dim() == 0)
    return true;

  // check for luck
  if (x.num_min_cons() != y.num_min_cons())
    return false;
  if (x.gens_info() != y.gens_info())
    return false;

  x.normalize();
  y.normalize();
  return (x.itvs == y.itvs) && (x.blocks == y.blocks) &&
    std::equal(x.factors.begin(), x.factors.end(),
               y.factors.begin(),
               std::mem_fn(&Factor::equals));
}

dim_type
F_Poly::affine_dim() const {
  dim_type res = 0;
  for (const auto& itv : proper_itvs())
    if (not itv.is_singleton())
      ++res;
  for (const auto& f : factors)
    res += f.affine_dim();
  return res;
}

void
F_Poly::set_topology(Topol t) {
  if (topology() == t)
    return;
  topol = t;
  // TODO: if Itv will support NNC, adjust their topology here.
  for (auto& f : factors) {
    assert(f.is_topologically_closed());
    f.set_topology(t);
  }
}

void
F_Poly::set_empty() {
  empty = true;
  itvs.clear();
  blocks.clear();
  factors.clear();
  is_normalized = true;
  assert(check_inv());
}

void
F_Poly::set_universe() {
  empty = false;
  itvs.assign(dim, Itv());
  blocks.clear();
  factors.clear();
  is_normalized = true;
  assert(check_inv());
}

namespace detail {

dim_type
hull_factors(Topol topol, Itvs& itvs1, Blocks& bs1, Factors& fs1,
             const Refactor_Proxy& fs2) {
  dim_type join_idx = not_a_dim();
  Index_Set joined_blocks;
  Factor fs2_join(0, topol);

  // First join blocks component.
  const auto nb = num_rows(bs1);
  assert(num_rows(fs1) == nb);
  for (auto i : range(nb)) {
    if (fs1[i] == fs2[i])
      continue;
    if (join_idx == not_a_dim()) {
      join_idx = i;
      fs2_join = fs2[i];
    } else {
      concat(bs1[join_idx], fs1[join_idx], bs1[i], fs1[i]);
      fs2_join.concatenate_assign(fs2[i]);
      joined_blocks.set(i);
    }
  }

  // Then join itvs component.
  for (auto i : index_range(itvs1)) {
    auto& itv1 = itvs1[i];
    if (itv1.is_empty() || itv1.is_universe())
      continue;
    if (not has_block_index(fs2.lub, Var(i)))
      continue;
    const auto& itv2 = fs2.get_itv(i);
    if (itv1 == itv2)
      continue;
    if (join_idx == not_a_dim()) {
      // Have to create new block/factor in bs1/fs1.
      join_idx = num_rows(bs1);
      bs1.push_back(Block());
      fs1.push_back(Factor(0, topol));
    }
    auto b2_dummy = bs1[join_idx];
    concat(bs1[join_idx], fs1[join_idx], i, itv1);
    concat(b2_dummy, fs2_join, i, itv2);
    // the i-th interval has been merged in a block
    itv1.set_empty();
  }

  if (join_idx != not_a_dim()) {
    erase_using_sorted_indices(bs1, joined_blocks);
    erase_using_sorted_indices(fs1, joined_blocks);
    // Safe: join_idx < all indices in joined_blocks
    fs1[join_idx].poly_hull_assign(fs2_join);
  }
  return join_idx;
}

} // namespace detail

void
F_Poly::poly_hull_assign(const F_Poly& y) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  if (y.is_empty())
    return;
  if (x.is_empty() || x.space_dim() == 0) {
    x = y;
    return;
  }

  // Check for unconstrained dims (speculative optimization)
  // FIXME, CHECKME: is factorization worth?
  x.factorize();
  y.factorize();
  for (auto i : dim_range(x)) {
    if (y.is_itv_dim(i)) {
      auto& y_itv = y.itvs[i];
      if (y_itv.is_universe())
        x.unconstrain(Var(i));
      else if (x.is_itv_dim(i)) {
        auto& x_itv = x.itvs[i];
        if ((x_itv.inf_lb() && y_itv.inf_ub())
            || (x_itv.inf_ub() && y_itv.inf_lb()))
          x_itv.set_universe();
      }
    }
  }

  using namespace detail;
  Blocks lub = blocks_lub(x.itvs, x.blocks, y.itvs, y.blocks);

  // FIXME: speculative optimization.
  // we could drop from lub's blocks dims that are unconstrained in x.
  // This however would require changes in the refactor_proxy.
  // Index_Set x_uncon;
  // for (auto i : dim_range(x)) {
  //   if (x.itvs[i].is_universe())
  //     x_uncon.set(i);
  // }
  // // Filter from lub those dims that are unconstrained.
  // erase_from_blocks(lub, x_uncon);

  x.sync(lub);
  Refactor_Proxy y_fs(y.impl(), lub);

  auto join_idx = hull_factors(x.topol, x.itvs, x.blocks, x.factors, y_fs);
  if (join_idx != not_a_dim()) {
    // check for refactoring (only for join_idx)
    detail::factorize(join_idx, x.blocks, x.factors, false);
    x.blocks_to_itvs();
  }
  assert(x.check_inv());
}

bool
F_Poly::constrains(Var v) const {
  if (is_empty())
    return true;
  auto i = v.id();
  if (is_itv_dim(i))
    return not itvs[i].is_universe();
  else {
    auto idx = detail::find_block_index(blocks, v);
    Var b_var = detail::convert(v, blocks[idx]);
    return factors[idx].constrains(b_var);
  }
}

bool
F_Poly::contains(const F_Poly& y) const {
  const auto& x = *this;
  assert(x.check_inv() && y.check_inv());

  if (x.is_empty())
    return y.is_empty();
  if (y.is_empty())
    return true;
  if (x.space_dim() == 0)
    return true;

  // Check common itv dims
  for (auto i : bwd_dim_range(x)) {
    if (x.is_itv_dim(i) && y.is_itv_dim(i)
        && not x.itvs[i].contains(y.itvs[i]))
      return false;
  }

  // Check x itv dims that are not itv dims in y (speculative)
  for (auto i : bwd_dim_range(x)) {
    if (x.is_itv_dim(i) && y.is_block_dim(i)
        && not x.itvs[i].is_universe()) {
      auto y_idx = detail::find_block_index(y.blocks, Var(i));
      Var v = detail::convert(Var(i), y.blocks[y_idx]);
      auto y_itv = y.factors[y_idx].get_bounds(v);
      if (not x.itvs[i].contains(y_itv))
        return false;
    }
  }

  Blocks lub = detail::blocks_lub(x.itvs, x.blocks, y.itvs, y.blocks);
  detail::Refactor_Proxy y_fs(y.impl(), lub);
  for (auto i : bwd_index_range(x.blocks)) {
    const auto& x_b = x.blocks[i];
    const auto& x_f = x.factors[i];
    for (const auto& x_c : x_f.cons()) {
      auto j = detail::get_containing_block_index(lub, x_b);
      const auto& y_b = lub[j];
      const auto& y_f = y_fs[j];
      auto y_c = detail::convert_back_and_forth(x_c, x_b, y_b);
      if (not y_f.relation_with(y_c).implies(Poly_Con_Rel::is_included()))
        return false;
    }
  }
  return true;
}

bool
F_Poly::boxed_contains(const F_Poly& y) const {
  const auto& x = *this;
  assert(x.check_inv() && y.check_inv());
  assert(x.get_bounding_box().contains(y.get_bounding_box()));

  if (x.space_dim() == 0 || x.is_empty() || y.is_empty())
    return true;

  // No need to (re-) check x itv dims, only check its block dims
  Blocks lub = detail::blocks_lub(x.itvs, x.blocks, y.itvs, y.blocks);
  detail::Refactor_Proxy y_fs(y.impl(), lub);
  for (auto i : bwd_index_range(x.blocks)) {
    const auto& x_b = x.blocks[i];
    const auto& x_f = x.factors[i];
    for (const auto& x_c : x_f.cons()) {
      if (!x_c.is_strict_inequality() && is_proper_interval_con(x_c))
        continue;
      auto j = detail::get_containing_block_index(lub, x_b);
      const auto& y_b = lub[j];
      const auto& y_f = y_fs[j];
      auto y_c = detail::convert_back_and_forth(x_c, x_b, y_b);
      if (not y_f.relation_with(y_c).implies(Poly_Con_Rel::is_included()))
        return false;
    }
  }
  return true;
}

bool
F_Poly::strictly_contains(const F_Poly& y) const {
  const auto& x = *this;
  assert(x.check_inv() && y.check_inv());
  return not x.equals(y) && x.contains(y);
}

bool
F_Poly::is_disjoint_from(const F_Poly& y) const {
  const auto& x = *this;
  assert(x.check_inv() && y.check_inv());
  if (x.is_empty() || y.is_empty())
    return true;
  if (x.space_dim() == 0)
    return false;

  // Speculative optimizations.
  auto nx = x.num_min_cons();
  if (nx < 2)
    return detail::simple_poly_is_disjoint_from_poly(nx, x, y);
  auto ny = y.num_min_cons();
  if (ny < 2)
    return detail::simple_poly_is_disjoint_from_poly(ny, y, x);

  // Check common itv dims
  for (auto i : bwd_dim_range(x)) {
    if (x.is_itv_dim(i) && y.is_itv_dim(i)
        && x.itvs[i].is_disjoint_from(y.itvs[i]))
      return true;
  }

  // Helper to check itv dims that are not in common
  // (speculative and only safe in the closed case)
  if (x.is_necessarily_closed()) {
    auto check_uncommon = [](const F_Poly& x, const F_Poly& y) {
      assert(x.is_necessarily_closed());
      for (auto i : bwd_dim_range(x)) {
        if (x.is_itv_dim(i) && y.is_block_dim(i)
            && not x.itvs[i].is_universe()) {
          auto y_idx = detail::find_block_index(y.blocks, Var(i));
          Var v = detail::convert(Var(i), y.blocks[y_idx]);
          auto y_itv = y.factors[y_idx].get_bounds(v);
          if (x.itvs[i].is_disjoint_from(y_itv))
            return true;
        }
      }
      return false;
    };
    // Check itv dims in x that are not itv dims in y (speculative)
    if (check_uncommon(x, y))
      return true;
    // Check itv dims in y that are not itv dims in x (speculative)
    if (check_uncommon(y, x))
      return true;
  }
  // General (expensive) case: intersection.
  F_Poly z = x;
  assert(z.check_inv() && y.check_inv());
  z.intersection_assign(y);
  return z.is_empty();
}

void
F_Poly::minimize() const {
  auto& x = const_cast<F_Poly&>(*this);
  for (auto& f : x.factors) {
    f.minimize();
    if (f.is_empty()) {
      x.set_empty();
      return;
    }
  }
}

namespace detail {

void ascii_dump_itv(const Itv& itv, std::ostream& os) {
  if (itv.is_empty()) {
    os << "block-dim";
    return;
  }
  os << "[ ";
  if (itv.inf_lb())
    os << "-inf";
  else {
    os << "lb ";
    itv.lb.ascii_dump(os);
  }
  os << " , ";
  if (itv.inf_ub())
    os << "+inf";
  else {
    os << "ub ";
    itv.ub.ascii_dump(os);
  }
  os << " ]";
}

bool ascii_load_itv(Itv& itv, std::istream& is) {
  std::string str;
  itv.set_universe();

  if (not (is >> str))
    return false;
  if (str == "block-dim") {
    itv.set_empty();
    return true;
  } else if (str != "[")
    return false;

  if (not (is >> str))
    return false;
  if (str == "lb") {
    Rational lb;
    if (not lb.ascii_load(is))
      return false;
    itv.set_lb(lb);
  } else if (str != "-inf")
    return false;

  if (not ascii_load_string(is, ","))
    return false;

  if (not (is >> str))
    return false;
  if (str == "ub") {
    Rational ub;
    if (not ub.ascii_load(is))
      return false;
    itv.set_ub(ub);
  } else if (str != "+inf")
    return false;

  if (not ascii_load_string(is, "]"))
    return false;

  return itv.check_inv();
}

void ascii_dump_block(const Block& b, std::ostream& os) {
  const auto dim = num_rows(b);
  os << dim << " : {";
  for (auto j : range(dim))
    os << " " << b[j];
  os << " }";
}

bool ascii_load_block(Block& b, std::istream& is) {
  dim_type dim = not_a_dim();
  if (not ((is >> dim) && (dim >= 0)))
    return false;
  b.resize(dim);
  std::string str;
  if (not (ascii_load_string(is, ":") && ascii_load_string(is, "{")))
    return false;
  for (auto j : range(dim)) {
    if (not ((is >> b[j]) && (b[j] >= 0)))
      return false;
  }
  if (not ascii_load_string(is, "}"))
    return false;
  return true;
}

} // namespace detail

void
F_Poly::ascii_dump(std::ostream& os) const {
  using namespace IO_Operators;
  os << "dim " << dim << "\n";
  os << "empty " << empty << "\n";
  os << "topol " << topol << "\n";
  os << "is_normalized " << is_normalized << "\n";

  os << "===itvs-component===\n";
  for (auto i : index_range(itvs)) {
    os << i << " : ";
    detail::ascii_dump_itv(itvs[i], os);
    os << "\n";
  }

  os << "===blocks-component===\n";
  os << "blocks " << num_rows(blocks) << std::endl;
  for (const auto& b : blocks) {
    detail::ascii_dump_block(b, os);
    os << "\n";
  }

  os << "factors " << num_rows(factors) << std::endl;
  for (auto i : index_range(factors)) {
    os << "===start-of-factor " << i << "\n";
    if (i >= num_rows(blocks))
      os << "missing block -- broken F_Poly ";
    else {
      os << "block = ";
      detail::ascii_dump_block(blocks[i], os);
    }
    os << "\n";
    factors[i].ascii_dump(os);
    os << "\n";
    os << "=== end of factor " << i << " ===\n";
  }
}

bool
F_Poly::ascii_load(std::istream& is) {
  std::string str;
  if (not (ascii_load_string(is, "dim") && (is >> dim) && (dim >= 0)))
    return false;
  if (not (ascii_load_string(is, "empty") && (is >> empty)))
    return false;
  if (not (ascii_load_string(is, "topol") && (is >> str)))
    return false;
  if (str == "CLOSED")
    topol = Topol::CLOSED;
  else if (str == "NNC")
    topol = Topol::NNC;
  else
    return false;
  if (not (ascii_load_string(is, "is_normalized") && (is >> is_normalized)))
    return false;

  if (not ascii_load_string(is, "===itvs-component==="))
      return false;
  itvs.resize(dim);
  for (auto i : range(dim)) {
    dim_type idx;
    if (not ((is >> idx) && (idx == i)))
      return false;
    if (not ascii_load_string(is, ":"))
      return false;
    if (not detail::ascii_load_itv(itvs[i], is))
      return false;
  }

  if (not ascii_load_string(is, "===blocks-component==="))
    return false;
  dim_type num_blocks = not_a_dim();
  if (not (ascii_load_string(is, "blocks") && (is >> num_blocks)))
    return false;
  if (num_blocks < 0)
    return false;
  blocks.resize(num_blocks);
  for (auto i : index_range(blocks)) {
    if (not detail::ascii_load_block(blocks[i], is))
      return false;
  }

  dim_type num_factors = not_a_dim();
  if (not (ascii_load_string(is, "factors") && (is >> num_factors)))
    return false;
  if (num_factors < 0)
    return false;
  factors.resize(num_factors);
  for (auto i : index_range(factors)) {
    dim_type index = not_a_dim();
    if (not (ascii_load_string(is, "===start-of-factor")
             && (is >> index) && (index == i)))
      return false;
    Block dummy;
    if (not (ascii_load_string(is, "(")
             && ascii_load_string(is, "block")
             && ascii_load_string(is, "=")
             && detail::ascii_load_block(dummy, is)
             && ascii_load_string(is, ")")))
      return false;
    if (dummy != blocks[i])
      return false;
    if (not factors[i].ascii_load(is))
      return false;
    if (not (ascii_load_string(is, "===end-of-factor")
             && (is >> index) && (index == i)))
      return false;
  }

  return check_inv();
}

void
F_Poly::add_space_dims(dim_type m, bool project) {
  assert(m >= 0);
  if (m == 0)
    return;
  if (empty) {
    dim += m;
    return;
  }
  Itv itv;
  if (project)
    itv.set_zero();
  itvs.insert(itvs.end(), m, itv);
  dim += m;
  assert(check_inv());
}

void
F_Poly::remove_higher_space_dims(dim_type new_dim) {
  assert(new_dim <= dim);
  if (empty) {
    dim = new_dim;
    return;
  }
  Index_Set vars;
  for (auto i : bwd_range(new_dim, space_dim()))
    vars.set(i);
  remove_space_dims(vars);
}

void
F_Poly::map_space_dims(const Dims& pfunc) {
  assert(space_dim() == num_rows(pfunc));
  if (space_dim() == 0)
    return;

  auto permute_dims = [this](const Dims& perm) {
    auto sd = space_dim();
    // apply perm to each itv
    Itvs new_itvs(sd);
    for (auto i : bwd_range(sd))
      new_itvs[perm[i]] = std::move(itvs[i]);
    std::swap(itvs, new_itvs);
    // apply perm to each block
    for (auto& b : blocks)
      for (auto& i : b)
        i = perm[i];
  };

  dim_type rem_count = std::count(pfunc.begin(), pfunc.end(), not_a_dim());

  if (empty) {
    dim -= rem_count;
    return;
  }

  if (rem_count == 0) {
    permute_dims(pfunc);
    assert(check_inv());
    return;
  }

  Index_Set tbr;
  for (auto i : bwd_dim_range(*this)) {
    if (pfunc[i] == not_a_dim())
      tbr.set(i);
  }
  remove_space_dims(tbr);
  // Recompute partial function after removal.
  Dims perm = pfunc;
  erase_using_sorted_indices(perm, tbr);
  permute_dims(perm);
  assert(check_inv());
}

void
F_Poly::fold_space_dims(const Index_Set& vars, Var dest) {
  assert(!vars.test(dest.id()));

  if (empty) {
    dim -= vars.size();
    return;
  }

  Block b(vars.begin(), vars.end());
  b.push_back(dest.id());

  if (are_itv_dims(b)) {
    Itv& itv_dest = itvs[dest.id()];
    for (auto v_i : vars)
      itv_dest.lub_assign(itvs[v_i]);
    remove_space_dims(vars);
    assert(check_inv());
    return;
  }

  auto i = merge(b);
  auto& b_i = blocks[i];
  auto& f_i = factors[i];

  Var c_dest = detail::convert(dest, b_i);

  Index_Set c_vars;
  for (auto v_j : vars) {
    Var c_var = detail::convert(Var(v_j), b_i);
    c_vars.set(c_var.id());
  }

  // Fold dimensions in factor f_i
  f_i.fold_space_dims(c_vars, c_dest);
  // Remove folded dimensions from b_i
  erase_using_sorted_indices(b_i, c_vars);

  // Adjust blocks
  detail::reduce_blocks(blocks, vars);

  // Remove folded dimensions from itvs
  erase_using_sorted_indices(itvs, vars);

  // Update space dimension
  dim -= vars.size();

  assert(check_inv());
}

void
F_Poly::expand_space_dim(Var v, dim_type m) {
  if (m == 0)
    return;
  if (empty) {
    dim += m;
    return;
  }

  itvs.reserve(dim + m);
  if (is_itv_dim(v.id())) {
    dim += m;
    const auto& itv = itvs[v.id()];
    itvs.resize(dim, itv);
    assert(check_inv());
    return;
  }

  auto idx = detail::find_block_index(blocks, v);
  const auto old_dim = dim;
  dim += m;
  itvs.resize(dim, Itv(Spec_Elem::EMPTY));
  if (blocks[idx].size() == 1) {
    // Add m copies of factor idx (and corresponding blocks).
    for (auto d : range(m))
      blocks.emplace_back(1, old_dim + d);
    // Note: this reserve avoids UB.
    factors.reserve(factors.size() + m);
    factors.insert(factors.end(), m, factors[idx]);
  } else {
    Var v_int = detail::convert(v, blocks[idx]);
    for (auto d : range(m))
      blocks[idx].push_back(old_dim + d);
    factors[idx].expand_space_dim(v_int, m);
    // refactor modified block
    detail::factorize(idx, blocks, factors, false);
  }
  blocks_to_itvs();
  assert(check_inv());
}

void
F_Poly::topological_closure_assign() {
  if (is_necessarily_closed())
    return;
  if (empty)
    return;
  for (auto& f : factors)
    f.topological_closure_assign();
  blocks_to_itvs();
  assert(check_inv());
}

void
F_Poly::unconstrain(const Index_Set& vars) {
  assert(vars.empty() || vars.last() < space_dim());
  if (empty || vars.empty())
    return;

  // Remove space dims from blocks/factors.
  Index_Set zerodims;
  for (auto idx : bwd_index_range(blocks)) {
    auto& b = blocks[idx];
    auto& f = factors[idx];
    Index_Set tbr;
    for (auto i : bwd_dim_range(f)) {
      if (vars.test(b[i]))
        tbr.set(i);
    }
    if (tbr.empty())
      continue;
    f.remove_space_dims(tbr);
    if (f.is_empty()) {
      set_empty();
      return;
    }
    if (tbr.size() == num_rows(b))
      zerodims.set(idx);
    else
      erase_using_sorted_indices(b, tbr);
  }
  // Remove 0-dim blocks/factors
  erase_using_sorted_indices(blocks, zerodims);
  erase_using_sorted_indices(factors, zerodims);

  // Now unconstrain space dims in itvs.
  for (auto i : vars)
    itvs[i].set_universe();

  blocks_to_itvs();
  is_normalized = false;
  assert(check_inv());
}

void
F_Poly::remove_space_dims(const Index_Set& vars) {
  if (is_empty()) {
    dim -= vars.size();
    return;
  }

  // Unconstrain dims (moving them to itvs)
  unconstrain(vars);
  // Remove dims from itvs (and adjust blocks)
  erase_using_sorted_indices(itvs, vars);
  detail::reduce_blocks(blocks, vars);

  dim -= vars.size();
  blocks_to_itvs();
  is_normalized = false;
  assert(check_inv());
}

void
F_Poly::intersection_assign(const F_Poly& y) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  if (x.empty)
    return;
  if (y.empty) {
    x.set_empty();
    return;
  }

  // Update common itv dims
  for (auto i : bwd_dim_range(x)) {
    if (x.is_itv_dim(i) && y.is_itv_dim(i)) {
      if (x.itvs[i].glb_assign(y.itvs[i])) {
        x.set_empty();
        return;
      }
    }
  }

  using namespace detail;
  Blocks lub = blocks_lub(x.itvs, x.blocks, y.itvs, y.blocks);
  x.sync(lub);
  Refactor_Proxy y_fs(y.impl(), lub);
  for (auto i : bwd_index_range(lub)) {
    auto& xi = x.factors[i];
    xi.intersection_assign(y_fs[i]);
    if (xi.is_empty()) {
      x.set_empty();
      return;
    }
  }
}

void
F_Poly::concatenate_assign(const F_Poly& y) {
  auto& x = *this;
  assert(x.topol == y.topol);
  if (x.empty || y.empty) {
    x.dim += y.dim;
    x.set_empty();
    return;
  }
  if (y.space_dim() == 0)
    return;
  if (x.space_dim() == 0) {
    x = y;
    return;
  }

  // Copy before modifying.
  const auto x_dim = x.dim;
  const auto x_nb = num_rows(x.blocks);

  x.dim += y.dim;
  x.itvs.insert(x.itvs.end(),
                y.itvs.begin(), y.itvs.end());
  x.blocks.insert(x.blocks.end(),
                  y.blocks.begin(), y.blocks.end());
  for (auto i : range(x_nb, num_rows(x.blocks))) {
    for (auto& d : x.blocks[i])
      d += x_dim;
  }
  x.factors.insert(x.factors.end(),
                   y.factors.begin(), y.factors.end());
  x.is_normalized &= y.is_normalized;
  assert(x.check_inv());
}

void
F_Poly::affine_image(Var var, const Linear_Expr& expr,
                     const Integer& inhomo, const Integer& den) {
  if (is_empty())
    return;

  Block b = detail::extract_block(expr);
  detail::add_var(b, var);
  auto i = merge(b);
  const auto& b_i = blocks[i];

  Var v_int = detail::convert(var, b_i);
  Linear_Expr le_int = detail::convert(expr, b_i);
  factors[i].affine_image(v_int, le_int, inhomo, den);
  assert(check_inv());
}

void
F_Poly::affine_preimage(Var var, const Linear_Expr& expr,
                        const Integer& inhomo, const Integer& den) {
  if (is_empty())
    return;

  Block b = detail::extract_block(expr);
  detail::add_var(b, var);
  auto i = merge(b);
  const auto& b_i = blocks[i];

  Var v_int = detail::convert(var, b_i);
  Linear_Expr le_int = detail::convert(expr, b_i);
  factors[i].affine_preimage(v_int, le_int, inhomo, den);
  assert(check_inv());
}

void
F_Poly::parallel_affine_image(const Vars& vars,
                              const Linear_Exprs& exprs,
                              const Integers& inhomos,
                              const Integers& dens) {
  detail::par_affine_image_aux(*this, vars, exprs, inhomos, dens);
}

void
F_Poly::widening_assign(const F_Poly& y, const Cons* upto_ptr,
                        Widen_Impl w_impl, Widen_Spec w_spec) {
  auto& x = *this;
  assert(x.topology() == y.topology());
  assert(x.space_dim() == y.space_dim());
  if (detail::widening_preamble(x, y, w_spec))
    return;

  if (w_spec == Widen_Spec::SAFE) {
    // Apply trivial lifting of risky widening.
    x.poly_hull_assign(y);
    x.minimize();
  }

  // FIXME: here we should check for certificate-based stability.
  // For now, we only check affine dimension.
  if (x.affine_dim() > y.affine_dim())
    return;

  Index_Set valid = detail::valid_upto_cons(*this, upto_ptr);

  for (auto i : index_range(itvs))
    if (x.is_itv_dim(i) && y.is_itv_dim(i))
      x.itvs[i].widen_assign(y.itvs[i]);

  // FIXME: temporary, could be improved?
  Blocks lub = detail::blocks_lub(x.itvs, x.blocks, y.itvs, y.blocks);
  x.sync(lub);
  detail::Refactor_Proxy y_fs(y.impl(), lub);
  for (auto i : bwd_index_range(lub))
    x.factors[i].widening_assign(y_fs[i], w_impl, Widen_Spec::RISKY);

  detail::add_valid_upto_cons(x, valid, upto_ptr);

  assert(x.check_inv());
}

BBox
F_Poly::get_bounding_box() const {
  if (is_empty())
    return BBox(space_dim(), Spec_Elem::EMPTY);

  BBox res(space_dim());
  res.itvs = itvs;
  for (auto i : bwd_index_range(blocks)) {
    BBox f_bbox = factors[i].get_bounding_box();
    const auto& block = blocks[i];
    // Move factor bounds into res.
    for (auto fd : bwd_index_range(block)) {
      auto sd = block[fd];
      res.itvs[sd] = std::move(f_bbox.itvs[fd]);
    }
  }
  res.maybe_update_volume_info();
  return res;
}

F_Poly
F_Poly::split_aux(const Con& c, Topol t, bool integral) {
  assert(not integral || t == Topol::CLOSED);
  assert(not c.is_equality() || integral);
  // Init result with empty poly
  auto res = F_Poly(space_dim(), Spec_Elem::EMPTY, topology());
  if (is_empty())
    return res;
  Block b = detail::extract_block(c);
  if (b.size() == 0) {
    if (c.is_tautological())
      return res;
    assert(c.is_inconsistent());
    m_swap(res);
    return res;
  }
  if (b.size() == 1) {
    assert(is_proper_interval_con(c));
    // Special case: check if it is an integral or non-strict split
    // on an interval dimension.
    const dim_type d = b[0];
    if ((integral || topol == Topol::CLOSED) && is_itv_dim(d)) {
      auto& itv_d = itvs[d];
      Itv res_d = split_itv(itv_d, c, integral);
      if (res_d.is_empty()) {
        if (itv_d.is_empty()) {
          /* both empty */
          set_empty();
        } else {
          /* only res empty: nothing to do */
        }
      } else {
        if (itv_d.is_empty()) {
          /* only this empty: swap */
          m_swap(res);
        } else {
          /* none empty: copy */
          res = *this;
        }
        res.itvs[d] = std::move(res_d);
      }
      return res;
    }
  }
  // Here we work on blocks and factors.
  auto i = merge(b);
  auto& bi = blocks[i];
  auto& fi = factors[i];
  assert(!fi.is_empty());
  Con c_bi = detail::convert(c, bi);
  auto res_fi = integral ? fi.integral_split(c_bi) : fi.split(c_bi, t);

  if (not fi.is_empty() && not res_fi.is_empty()) {
    // Both non-empty
    is_normalized = false;
    res = *this;
    res.factors[i] = std::move(res_fi);
    assert(check_inv() && res.check_inv());
    return res;
  }
  if (res_fi.is_empty()) {
    // Note: integral split can make both empty.
    if (integral && fi.is_empty())
      set_empty();
    return res;
  }
  assert(fi.is_empty());
  // Note: integral split can make both empty.
  if (integral && res_fi.is_empty())
    set_empty();
  else {
    m_swap(res);
    res.factors[i] = std::move(res_fi);
    res.is_normalized = false;
  }
  return res;
}

void
F_Poly::time_elapse_assign(const F_Poly& y) {
  auto& x = *this;
  assert(x.space_dim() == y.space_dim());
  if (x.is_empty())
    return;
  if (y.is_empty()) {
    x.set_empty();
    return;
  }
  if (x.space_dim() == 0)
    return;

  Gens y_rays;
  detail::add_as_rays(y.copy_gens(), y_rays);
  Index_Set y_vars;
  for (const auto& g : y_rays) {
    const auto& ex = g.linear_expr();
    for (auto i : bwd_dim_range(ex))
      if (!ex.get(i).is_zero())
        y_vars.set(i);
  }
  // Check if y was (just) the origin of the vector space.
  if (y_vars.empty())
    return;

  Block b(y_vars.begin(), y_vars.end());
  auto i = merge(b);
  auto& x_bi = x.blocks[i];
  auto& x_fi = x.factors[i];
  for (const auto& g : y_rays)
    x_fi.add_gen(detail::convert(g, x_bi));
  if (x_fi.space_dim() > 1)
    detail::factorize(i, x.blocks, x.factors, false);
  x.blocks_to_itvs();
  assert(x.check_inv());
}

void
F_Poly::normalize() const {
  if (is_normalized)
    return;
  // Only semantically const
  auto& x = const_cast<F_Poly&>(*this);
  // Compute maximal factorization
  x.factorize(true);
  // Sort blocks (and factors) lexicographically
  detail::synchronous_sort(x.blocks, x.factors);
  x.is_normalized = true;
}

size_t
F_Poly::hash() const {
  normalize();
  size_t res = hash_size(space_dim());
  if (is_empty())
    return res;
  for (const auto& itv : itvs)
    hash_combine(res, itv.hash());
  hash_combine(res, hash_size(num_rows(blocks)));
  for (const auto& b : blocks) {
    for (auto i : b)
      hash_combine(res, i);
  }
  for (const auto& f : factors)
    hash_combine(res, f.hash());
  return res;
}

void
F_Poly::print(std::ostream& os) const {
  // FIXME
#if 1
  Poly ph = to_poly();
  ph.minimize();
  ph.print(os);
#else
  if (is_empty())
    os << "false";
  minimize();
  bool comma = false;
  for (const auto& c : cons()) {
    if (comma) os << ", ";
    using namespace IO_Operators;
    os << c;
    comma = true;
  }
#endif
}

} // namespace pplite
