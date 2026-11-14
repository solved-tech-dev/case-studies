// SPDX-License-Identifier: AGPL-3.0-only
#include <cmath>
#include <petscmat.h>
#include <petscvec.h>

static PetscErrorCode exercise() {
  Mat a, f, b, x;
  Vec rhs, sol;
  MatFactorInfo info;
  IS r, c;
  PetscFunctionBeginUser;
  PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF, 3, 3, 3, nullptr, &a));
  PetscCall(MatSetOption(a, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE));
  PetscInt indices[] = {0, 1, 2};
  PetscScalar values[] = {4, 1, 0, 2, 5, 1, 0, 3, 6};
  PetscCall(MatSetValues(a, 3, indices, 3, indices, values, INSERT_VALUES));
  PetscCall(MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY));
  PetscCall(MatFactorInfoInitialize(&info));
  PetscCall(ISCreateStride(PETSC_COMM_SELF, 3, 0, 1, &r));
  PetscCall(ISDuplicate(r, &c));
  PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, 3, 2, nullptr, &b));
  PetscCall(MatDuplicate(b, MAT_DO_NOT_COPY_VALUES, &x));
  PetscScalar *ba = nullptr;
  PetscCall(MatDenseGetArray(b, &ba));
  const double input[] = {5, 8, 9, 10, 16, 18};
  for (int j = 0; j < 6; ++j)
    ba[j] = input[j];
  PetscCall(MatDenseRestoreArray(b, &ba));
  for (int pass = 0; pass < 2; ++pass) {
    PetscCall(MatGetFactor(a, "pangulu", MAT_FACTOR_LU, &f));
    PetscCall(MatLUFactorSymbolic(f, a, r, c, &info));
    PetscCall(MatLUFactorNumeric(f, a, &info));
    PetscCall(MatMatSolve(f, b, x));
    const PetscScalar *xa = nullptr;
    PetscCall(MatDenseGetArrayRead(x, &xa));
    for (int j = 0; j < 6; ++j)
      PetscCheck(std::abs(xa[j] - (j < 3 ? 1.0 : 2.0)) < 1e-10, PETSC_COMM_SELF, PETSC_ERR_PLIB,
                 "Block solution mismatch");
    PetscCall(MatDenseRestoreArrayRead(x, &xa));
    PetscCall(MatDestroy(&f));
  }
  PetscCall(MatSetValue(a, 0, 2, 2, INSERT_VALUES));
  PetscCall(MatSetValue(a, 1, 2, 0, INSERT_VALUES));
  PetscCall(MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY));
  PetscCall(MatGetFactor(a, "pangulu", MAT_FACTOR_LU, &f));
  PetscCall(MatLUFactorSymbolic(f, a, r, c, &info));
  PetscCall(MatLUFactorNumeric(f, a, &info));
  PetscCall(VecCreateSeq(PETSC_COMM_SELF, 3, &rhs));
  PetscCall(VecDuplicate(rhs, &sol));
  PetscScalar bv[] = {7, 7, 9};
  PetscCall(VecSetValues(rhs, 3, indices, bv, INSERT_VALUES));
  PetscCall(VecAssemblyBegin(rhs));
  PetscCall(VecAssemblyEnd(rhs));
  PetscCall(MatSolve(f, rhs, sol));
  const PetscScalar *sv = nullptr;
  PetscCall(VecGetArrayRead(sol, &sv));
  for (int j = 0; j < 3; ++j)
    PetscCheck(std::abs(sv[j] - 1) < 1e-10, PETSC_COMM_SELF, PETSC_ERR_PLIB,
               "Updated solution mismatch");
  PetscCall(VecRestoreArrayRead(sol, &sv));
  PetscCall(VecDestroy(&rhs));
  PetscCall(VecDestroy(&sol));
  PetscCall(MatDestroy(&f));
  PetscCall(MatDestroy(&a));
  PetscCall(MatDestroy(&b));
  PetscCall(MatDestroy(&x));
  PetscCall(ISDestroy(&r));
  PetscCall(ISDestroy(&c));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "PETSC_CALLBACK_TEST_PASSED\n"));
  PetscFunctionReturn(PETSC_SUCCESS);
}
int main(int argc, char **argv) {
  PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));
  PetscCall(exercise());
  PetscCall(PetscFinalize());
  return 0;
}
