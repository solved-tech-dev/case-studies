"""Shared one-RHS PanguLU bridge and independent FEALPy solution checks."""
import argparse
import ctypes
from pathlib import Path

import numpy as np


def arguments(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--small", action="store_true", help="use the small verification mesh")
    parser.add_argument("--output", type=Path, help="save validated FE coefficients to a .npy file")
    args = parser.parse_args()
    if args.output and args.output.suffix != ".npy":
        parser.error("--output expects a .npy path")
    return args


def solve(matrix, rhs):
    a = matrix.to_scipy().tocsr()
    a.sum_duplicates()
    a.sort_indices()
    n = a.shape[0]
    if a.shape != (n, n) or not 0 < n <= np.iinfo(np.int32).max:
        raise ValueError("expected a square matrix with 32-bit column indices")
    if a.dtype not in (np.float64, np.complex128):
        raise TypeError(f"unsupported matrix type {a.dtype}")
    b = np.ascontiguousarray(rhs)
    if b.dtype != a.dtype or b.shape != (n,):
        raise ValueError("expected one matching double-precision RHS")
    if not (np.isfinite(a.data).all() and np.isfinite(b).all()):
        raise ValueError("non-finite coefficient or RHS")
    ptr = np.ascontiguousarray(a.indptr, dtype=np.int64)
    col = np.ascontiguousarray(a.indices, dtype=np.int32)
    values = np.ascontiguousarray(a.data)
    precision = "complex" if np.iscomplexobj(values) else "real"
    path = Path(__file__).resolve().parents[1] / "build" / f"libfealpy_pangulu_{precision}.so"
    lib = ctypes.CDLL(str(path))
    native = lib.solve
    native.argtypes = [ctypes.c_int] + [ctypes.c_void_p] * 5
    native.restype = ctypes.c_int
    x = np.empty_like(b)
    if native(n, *(arr.ctypes.data for arr in (ptr, col, values, b, x))):
        raise RuntimeError("PanguLU solve failed; see the native error above")
    residual = np.linalg.norm(a @ x - b) / np.linalg.norm(b)
    print(f"equations={n} csr_entries={a.nnz} relative_residual={residual:.8e}")
    if not np.isfinite(residual) or residual > 1e-8:
        raise ValueError("original CSR residual exceeds tolerance")
    return x


def validate_field(field, error, label, tolerance, output):
    print(f"{label}={error:.8e}")
    coefficients = np.asarray(field.array)
    if not np.isfinite(coefficients).all() or not np.isfinite(error) or error > tolerance:
        raise ValueError("finite-element field validation failed")
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        np.save(output, coefficients)
    print("validation=passed")
