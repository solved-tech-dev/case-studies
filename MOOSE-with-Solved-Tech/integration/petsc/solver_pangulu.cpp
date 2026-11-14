// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Solved Tech
#include "pangulu_solver.hpp"
#include <petsc/private/matimpl.h>
#include <petscvec.h>
#if defined(PETSC_USE_COMPLEX) || !defined(PETSC_USE_REAL_DOUBLE)
#error This adapter requires real double-precision PETSc.
#endif
namespace {
constexpr const char *solver_name = "pangulu";
using Shared = std::shared_ptr<solved::Solver>;
struct State {
  Shared backend;
  PetscErrorCode (*base_destroy)(Mat) = nullptr;
  double tolerance = 1e-8;
};
PetscErrorCode DeleteCache(PetscCtxRt ptr) {
  auto **p = static_cast<Shared **>(ptr);
  delete static_cast<Shared *>(*p);
  *p = nullptr;
  return PETSC_SUCCESS;
}
PetscErrorCode Destroy(Mat f) {
  auto *s = static_cast<State *>(f->spptr);
  auto destroy = s ? s->base_destroy : nullptr;
  PetscFunctionBegin;
  delete s;
  f->spptr = nullptr;
  PetscCall(PetscObjectComposeFunction((PetscObject)f, "MatFactorGetSolverType_C", nullptr));
  if (destroy)
    PetscCall(destroy(f));
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode Solve(Mat f, Vec b, Vec x) {
  auto *s = static_cast<State *>(f->spptr);
  const PetscScalar *ba = nullptr;
  PetscScalar *xa = nullptr;
  std::string error;
  PetscInt n, m;
  PetscFunctionBegin;
  PetscCall(VecGetSize(b, &n));
  PetscCall(VecGetSize(x, &m));
  PetscCheck(n == s->backend->matrix.n && m == n, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
             "RHS size mismatch");
  PetscCall(VecGetArrayRead(b, &ba));
  PetscCall(VecGetArray(x, &xa));
  try {
    s->backend->solve(ba, xa, 1, s->tolerance);
  } catch (const std::exception &e) {
    error = e.what();
  }
  PetscCall(VecRestoreArray(x, &xa));
  PetscCall(VecRestoreArrayRead(b, &ba));
  PetscCheck(error.empty(), PETSC_COMM_SELF, PETSC_ERR_LIB, "%s", error.c_str());
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode MatSolveBlock(Mat f, Mat b, Mat x) {
  auto *s = static_cast<State *>(f->spptr);
  PetscInt n, k, m, l, ldb, ldx;
  const PetscScalar *ba = nullptr;
  PetscScalar *xa = nullptr;
  std::string error;
  PetscFunctionBegin;
  PetscCall(MatGetSize(b, &n, &k));
  PetscCall(MatGetSize(x, &m, &l));
  PetscCheck(n == s->backend->matrix.n && m == n && l == k, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
             "Dense RHS size mismatch");
  PetscCall(MatDenseGetLDA(b, &ldb));
  PetscCall(MatDenseGetLDA(x, &ldx));
  std::vector<double> rhs(size_t(n) * k), sol(rhs.size());
  PetscCall(MatDenseGetArrayRead(b, &ba));
  for (PetscInt j = 0; j < k; ++j)
    std::copy_n(ba + j * ldb, n, rhs.data() + j * n);
  PetscCall(MatDenseRestoreArrayRead(b, &ba));
  try {
    s->backend->solve(rhs.data(), sol.data(), k, s->tolerance);
  } catch (const std::exception &e) {
    error = e.what();
  }
  PetscCheck(error.empty(), PETSC_COMM_SELF, PETSC_ERR_LIB, "%s", error.c_str());
  PetscCall(MatDenseGetArray(x, &xa));
  for (PetscInt j = 0; j < k; ++j)
    std::copy_n(sol.data() + j * n, n, xa + j * ldx);
  PetscCall(MatDenseRestoreArray(x, &xa));
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode Numeric(Mat f, Mat a, const MatFactorInfo *) {
  const PetscInt *rp = nullptr, *ci = nullptr;
  const PetscScalar *v = nullptr;
  PetscBool done;
  PetscInt n;
  std::vector<int64_t> p;
  std::vector<int32_t> c;
  std::vector<double> values;
  std::string error;
  auto *s = static_cast<State *>(f->spptr);
  PetscFunctionBegin;
  PetscCall(MatGetRowIJ(a, 0, PETSC_FALSE, PETSC_FALSE, &n, &rp, &ci, &done));
  PetscCheck(done, PETSC_COMM_SELF, PETSC_ERR_SUP, "Sequential AIJ CSR access is required");
  PetscCall(MatSeqAIJGetArrayRead(a, &v));
  try {
    solved::require(n <= INT32_MAX, "Matrix row indices exceed 32 bits");
    p.assign(rp, rp + n + 1);
    c.resize(size_t(rp[n]));
    values.assign(v, v + rp[n]);
    for (PetscInt j = 0; j < rp[n]; ++j) {
      solved::require(ci[j] >= 0 && ci[j] < n, "Column out of range");
      c[j] = int32_t(ci[j]);
    }
  } catch (const std::exception &e) {
    error = e.what();
  }
  const PetscInt matrix_n = n; // MatRestoreRowIJ resets n to zero.
  PetscCall(MatSeqAIJRestoreArrayRead(a, &v));
  PetscCall(MatRestoreRowIJ(a, 0, PETSC_FALSE, PETSC_FALSE, &n, &rp, &ci, &done));
  PetscCheck(error.empty(), PETSC_COMM_SELF, PETSC_ERR_LIB, "%s", error.c_str());
  try {
    s->backend->factor(matrix_n, p, c, values);
  } catch (const std::exception &e) {
    error = e.what();
  }
  PetscCheck(error.empty(), PETSC_COMM_SELF, PETSC_ERR_LIB, "%s", error.c_str());
  f->ops->solve = Solve;
  f->ops->matsolve = MatSolveBlock;
  f->assembled = PETSC_TRUE;
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode Symbolic(Mat f, Mat, IS, IS, const MatFactorInfo *) {
  PetscFunctionBegin;
  f->ops->lufactornumeric = Numeric;
  PetscFunctionReturn(PETSC_SUCCESS);
}
PetscErrorCode Type(Mat, MatSolverType *type) {
  *type = solver_name;
  return PETSC_SUCCESS;
}
PetscErrorCode GetFactor(Mat a, MatFactorType type, Mat *f) {
  PetscInt n, m;
  PetscMPIInt size;
  PetscContainer cache = nullptr;
  Shared *backend = nullptr;
  PetscFunctionBegin;
  PetscCallMPI(MPI_Comm_size(PetscObjectComm((PetscObject)a), &size));
  PetscCheck(size == 1 && type == MAT_FACTOR_LU, PETSC_COMM_SELF, PETSC_ERR_SUP,
             "PanguLU requires one MPI rank and general LU");
  PetscCall(MatGetSize(a, &n, &m));
  PetscCheck(n == m, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Matrix must be square");
  // Cache belongs to the source matrix, so recreating a PETSc factor object does
  // not discard a valid analysis. Every numeric request still compares the graph
  // and all coefficients; PETSc's reuse hints alone never decide what to reuse.
  PetscCall(PetscObjectQuery((PetscObject)a, "solved_pangulu_cache", (PetscObject *)&cache));
  if (!cache) {
    backend = new Shared(std::make_shared<solved::Solver>());
    PetscCall(PetscContainerCreate(PETSC_COMM_SELF, &cache));
    PetscCall(PetscContainerSetPointer(cache, backend));
    PetscCall(PetscContainerSetCtxDestroy(cache, DeleteCache));
    PetscCall(PetscObjectCompose((PetscObject)a, "solved_pangulu_cache", (PetscObject)cache));
    PetscCall(PetscContainerDestroy(&cache));
  } else {
    void *ptr = nullptr;
    PetscCall(PetscContainerGetPointer(cache, &ptr));
    backend = static_cast<Shared *>(ptr);
  }
  // Never share a mutable numeric factor with another live PETSc factor.
  if (backend->use_count() > 1)
    *backend = std::make_shared<solved::Solver>();
  PetscCall(MatCreate(PETSC_COMM_SELF, f));
  PetscCall(MatSetSizes(*f, n, n, n, n));
  PetscCall(MatSetType(*f, MATSEQAIJ));
  PetscCall(MatSeqAIJSetPreallocation(*f, 0, nullptr));
  auto *s = new State();
  s->backend = *backend;
  s->base_destroy = (*f)->ops->destroy;
  PetscCall(PetscOptionsGetReal(nullptr, nullptr, "-pangulu_residual_tol", &s->tolerance, nullptr));
  PetscCheck(std::isfinite(s->tolerance) && s->tolerance > 0, PETSC_COMM_SELF,
             PETSC_ERR_ARG_OUTOFRANGE, "Invalid residual tolerance");
  (*f)->spptr = s;
  (*f)->ops->destroy = Destroy;
  (*f)->ops->lufactorsymbolic = Symbolic;
  (*f)->trivialsymbolic = PETSC_TRUE;
  (*f)->factortype = MAT_FACTOR_LU;
  (*f)->assembled = PETSC_TRUE;
  (*f)->preallocated = PETSC_TRUE;
  PetscCall(PetscFree((*f)->solvertype));
  PetscCall(PetscStrallocpy(solver_name, &(*f)->solvertype));
  PetscCall(PetscObjectComposeFunction((PetscObject)*f, "MatFactorGetSolverType_C", Type));
  PetscFunctionReturn(PETSC_SUCCESS);
}
} // namespace
PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_moose_pangulu(void) {
  PetscFunctionBegin;
  PetscCall(MatSolverTypeRegister(solver_name, MATSEQAIJ, MAT_FACTOR_LU, GetFactor));
  PetscFunctionReturn(PETSC_SUCCESS);
}
