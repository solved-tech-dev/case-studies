// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Solved Tech
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>
namespace solved {
inline void require(bool ok, const char *msg) {
  if (!ok)
    throw std::runtime_error(msg);
}
struct Sparse {
  int64_t n = 0;
  std::vector<int64_t> rp, cp, permutation;
  std::vector<int32_t> ci, ri;
  std::vector<double> csr, csc;
  bool same_graph(int64_t size, const std::vector<int64_t> &p,
                  const std::vector<int32_t> &c) const {
    return n == size && rp == p && ci == c;
  }
  void reset(int64_t size, const std::vector<int64_t> &p, const std::vector<int32_t> &c) {
    require(size > 0 && size <= INT32_MAX && p.size() == size_t(size + 1),
            "Invalid matrix dimension");
    require(p.front() == 0 && p.back() == int64_t(c.size()), "Invalid CSR pointers");
    for (int64_t i = 0; i < size; ++i)
      require(p[i] >= 0 && p[i] <= p[i + 1] && p[i + 1] <= int64_t(c.size()),
              "Invalid CSR pointer range");
    for (int64_t i = 0; i < size; ++i) {
      int32_t last = -1;
      for (auto j = p[i]; j < p[i + 1]; ++j) {
        require(c[j] > last && c[j] < size, "CSR columns must be sorted, unique and in range");
        last = c[j];
      }
    }
    n = size;
    rp = p;
    ci = c;
    cp.assign(n + 1, 0);
    ri.resize(c.size());
    permutation.resize(c.size());
    csr.assign(c.size(), 0);
    csc.assign(c.size(), 0);
    for (auto j : c)
      ++cp[j + 1];
    for (int64_t j = 0; j < n; ++j)
      cp[j + 1] += cp[j];
    auto cursor = cp;
    for (int64_t i = 0; i < n; ++i)
      for (auto j = rp[i]; j < rp[i + 1]; ++j) {
        auto q = cursor[ci[j]]++;
        ri[q] = int32_t(i);
        permutation[j] = q;
      }
  }
  bool refresh(const std::vector<double> &a) {
    require(a.size() == ci.size(), "Wrong coefficient count");
    bool changed = false;
    for (size_t j = 0; j < a.size(); ++j) {
      require(std::isfinite(a[j]), "Nonfinite coefficient");
      changed |= (csr[j] != a[j]);
    }
    // Preserve explicit zeros and update every value, including newly zero entries.
    std::copy(a.begin(), a.end(), csr.begin());
    for (size_t j = 0; j < a.size(); ++j)
      csc[permutation[j]] = a[j];
    return changed;
  }
  double residual(const double *b, const double *x) const {
    long double rr = 0, bb = 0;
    for (int64_t i = 0; i < n; ++i) {
      require(std::isfinite(b[i]) && std::isfinite(x[i]), "Nonfinite RHS or solution");
      long double ax = 0;
      for (auto j = rp[i]; j < rp[i + 1]; ++j)
        ax += static_cast<long double>(csr[j]) * x[ci[j]];
      const long double d = ax - b[i];
      rr += d * d;
      bb += static_cast<long double>(b[i]) * b[i];
    }
    const auto value = static_cast<double>(std::sqrt(bb > 0 ? rr / bb : rr));
    require(std::isfinite(value), "Nonfinite residual");
    return value;
  }
};
} // namespace solved
