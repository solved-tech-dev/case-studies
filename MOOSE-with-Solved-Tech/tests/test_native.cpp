// SPDX-License-Identifier: AGPL-3.0-only
// Runs on a licensed CUDA GPU. Checks actual phase reuse and native block RHS.
#include "pangulu_solver.hpp"
#include <iostream>
int main() {
  try {
    solved::Solver s;
    std::vector<int64_t> p = {0, 3, 6, 9};
    std::vector<int32_t> c = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<double> a = {4, 1, 0, 2, 5, 1, 0, 3, 6};
    s.factor(3, p, c, a);
    double b[] = {5, 8, 9, 10, 16, 18}, x[6]{};
    s.solve(b, x, 2);
    s.factor(3, p, c, a);
    s.solve(b, x, 2);
    solved::require(s.analyses == 1 && s.numerics == 1, "Identical A did not reuse factors");
    a[2] = 2;
    a[5] = 0;
    double b2[] = {7, 7, 9};
    s.factor(3, p, c, a);
    s.solve(b2, x);
    solved::require(s.analyses == 1 && s.numerics == 2, "Changed values did not reuse analysis");
    s.factor(3, {0, 1, 2, 3}, {0, 1, 2}, {2, 3, 4});
    double b3[] = {2, 3, 4};
    s.solve(b3, x);
    solved::require(s.analyses == 2 && s.numerics == 3, "Changed graph did not trigger analysis");
    std::cout << "NATIVE_REUSE_TEST_PASSED\n";
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
