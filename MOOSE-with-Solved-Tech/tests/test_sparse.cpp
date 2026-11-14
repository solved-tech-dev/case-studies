// SPDX-License-Identifier: AGPL-3.0-only
#include "sparse.hpp"
#include <cassert>
#include <iostream>
int main() {
  solved::Sparse a;
  a.reset(3, {0, 2, 5, 7}, {0, 2, 0, 1, 2, 1, 2});
  a.refresh({4, 0, 2, 5, 1, 3, 6});
  assert(a.cp == std::vector<int64_t>({0, 2, 4, 7}));
  assert(a.csc == std::vector<double>({4, 2, 5, 3, 0, 1, 6}));
  assert(a.ci.size() == 7);
  double x[] = {1, 2, 3}, b[] = {4, 15, 24};
  assert(a.residual(b, x) == 0);
  assert(!a.refresh({4, 0, 2, 5, 1, 3, 6}));
  assert(a.refresh({4, 7, 2, 5, 0, 3, 6}));
  assert(a.csc[4] == 7 && a.csc[5] == 0); // zero-to-nonzero and nonzero-to-zero
  bool rejected = false;
  try {
    a.reset(2, {0, 2, 2}, {0, 0});
  } catch (...) {
    rejected = true;
  }
  assert(rejected);
  rejected = false;
  try {
    a.reset(2, {0, 100, 2}, {0, 1});
  } catch (...) {
    rejected = true;
  }
  assert(rejected); // reject an out-of-bounds row before reading its columns
  std::cout << "Sparse conversion, exact zeros, value updates and original CSR residual: passed\n";
}
