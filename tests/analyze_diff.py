#!/usr/bin/env python3
# /// script
# dependencies = [
#   "numpy",
# ]
# ///

import numpy as np
import sys

def analyze_file_diff(orig_file, omp_file, name):
    print(f"\n{'='*60}")
    print(f"Analyzing: {name}")
    print(f"{'='*60}")

    orig = np.loadtxt(orig_file)
    omp = np.loadtxt(omp_file)

    diff = np.abs(orig - omp)
    max_diff = np.max(diff, axis=0)

    print(f"\nMaximum absolute differences per column:")
    for i, md in enumerate(max_diff):
        print(f"  Column {i+1}: {md:.10e}")

    # Find where maximum difference occurs
    idx = np.unravel_index(np.argmax(diff), diff.shape)
    print(f"\nLargest difference at row {idx[0]+1}, column {idx[1]+1}:")
    print(f"  Original: {orig[idx]:.15f}")
    print(f"  OMP:      {omp[idx]:.15f}")
    print(f"  Diff:     {diff[idx]:.15e}")

    # Relative error analysis (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_diff = np.abs((orig - omp) / np.where(np.abs(orig) > 1e-20, orig, 1e-20))
        rel_diff = np.where(np.abs(orig) > 1e-20, rel_diff, 0)

    max_rel = np.max(rel_diff, axis=0)
    print(f"\nMaximum relative differences per column:")
    for i, mr in enumerate(max_rel):
        print(f"  Column {i+1}: {mr:.10e}")

    # Count how many rows have any difference
    rows_with_diff = np.any(diff > 0, axis=1)
    num_diff_rows = np.sum(rows_with_diff)
    print(f"\nRows with differences: {num_diff_rows} / {len(orig)} ({100*num_diff_rows/len(orig):.2f}%)")

    # Show first few rows with differences
    if num_diff_rows > 0:
        diff_indices = np.where(rows_with_diff)[0]
        print(f"\nFirst 5 rows with differences:")
        for i in diff_indices[:5]:
            print(f"  Row {i+1}: max diff = {np.max(diff[i]):.10e}")

if __name__ == "__main__":
    analyze_file_diff("orig/snl_omp.dat", "omp/snl_omp.dat", "snl_omp.dat")
    analyze_file_diff("orig/sno2_omp.dat", "omp/sno2_omp.dat", "sno2_omp.dat")
    analyze_file_diff("orig/noise_omp.dat", "omp/noise_omp.dat", "noise_omp.dat")
    analyze_file_diff("orig/snc_omp.dat", "omp/snc_omp.dat", "snc_omp.dat")
