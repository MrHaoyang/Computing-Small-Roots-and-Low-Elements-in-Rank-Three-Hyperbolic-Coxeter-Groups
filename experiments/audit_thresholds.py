"""Targeted post-review diagnostic; does not modify recorded benchmark data.

Run: python audit_thresholds.py --output threshold_audit_rerun.json
Only the three level-one root sets are enumerated. No low-element search,
timing benchmark, or rendering is run. Decimal checks are high-precision
numerical evidence, not exact-arithmetic certificates.
"""
from collections import deque
from decimal import Decimal, localcontext
from pathlib import Path
import argparse
import json
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reference_algebra import CoxeterSystem3


def in_strict_band(value, tolerance):
    width = 100 * tolerance
    return 1 - width <= value < 1 or -1 < value <= -1 + width


def decimal_pairing(parameters, seed, path, reflection):
    """Reconstruct a root from parameters, not from rounded root coefficients."""
    with localcontext() as context:
        context.prec = 80
        one, two = Decimal(1), Decimal(2)

        def cosine(order):
            if order in ("inf", float("inf")):
                return one
            if order == 3:
                return one / two
            if order == 4:
                return two.sqrt() / two
            if order == 5:
                return (one + Decimal(5).sqrt()) / Decimal(4)
            raise ValueError("This targeted audit supports orders 3,4,5,inf only.")

        gram = [[Decimal(int(i == j)) for j in range(3)] for i in range(3)]
        for (i, j), order in zip(((0, 1), (0, 2), (1, 2)), parameters, strict=True):
            gram[i][j] = gram[j][i] = -cosine(order)
        beta = [Decimal(int(i == seed)) for i in range(3)]
        for index in path:
            pairing = sum((gram[index][j] * beta[j] for j in range(3)), Decimal(0))
            beta[index] -= two * pairing
        return +sum((gram[reflection][j] * beta[j] for j in range(3)), Decimal(0))


def run_audit():
    tolerance = 1e-9
    cases = []
    for parameters in ((3, 4, 5), (4, 4, 4), (4, 4, "inf")):
        system = CoxeterSystem3(*parameters, tolerance=tolerance)
        roots = system.enumerate_m_small_roots(1)
        paths = {}
        pending = deque()
        for seed, vector in enumerate(system.simple_roots):
            key = system.root_key(vector)
            paths[key] = (seed, ())
            pending.append(key)
        while pending:
            key = pending.popleft()
            beta = roots[key].coefficients
            seed, path = paths[key]
            for index, matrix in enumerate(system.generator_matrices.values()):
                reflected = matrix @ beta
                if np.all(reflected <= 10 * tolerance):
                    continue
                other = system.root_key(reflected)
                if other in roots and other not in paths:
                    paths[other] = (seed, path + (index,))
                    pending.append(other)
        if set(paths) != set(roots):
            raise AssertionError("Not all accepted roots have a recovered path.")
        positive, excluded, near, hits, mismatches = 0, 0, [], [], 0
        for key, record in roots.items():
            beta = record.coefficients
            seed, path = paths[key]
            for index, matrix in enumerate(system.generator_matrices.values()):
                reflected = matrix @ beta
                if np.all(reflected <= 10 * tolerance):
                    excluded += 1
                    continue
                positive += 1
                value = float(beta.T @ system.gram @ system.simple_roots[index])
                high = decimal_pairing(parameters, seed, path, index)
                if abs(high - Decimal(1)) < Decimal("1e-60"):
                    mathematical_step = -1
                elif abs(high + Decimal(1)) < Decimal("1e-60"):
                    mathematical_step = 1
                else:
                    mathematical_step = -1 if high >= 1 else (1 if high <= -1 else 0)
                numerical_step = -1 if value >= 1 - 100 * tolerance else (
                    1 if value <= -1 + 100 * tolerance else 0
                )
                mismatches += int(numerical_step != mathematical_step)
                item = {
                    "root_coefficients": beta.tolist(),
                    "seed_simple_root": seed + 1,
                    "reflection_path": [i + 1 for i in path],
                    "reflection": index + 1,
                    "binary64_pairing": value,
                    "decimal80_pairing": str(high),
                    "numerical_depth_change": numerical_step,
                    "high_precision_depth_change": mathematical_step,
                }
                if in_strict_band(value, tolerance):
                    hits.append(item)
                if min(abs(value - 1), abs(value + 1)) <= 100 * tolerance:
                    near.append(item)
        cases.append({
            "parameters": list(parameters),
            "level": 1,
            "root_count": len(roots),
            "all_candidates": 3 * len(roots),
            "excluded_negative_candidates": excluded,
            "positive_candidates": positive,
            "strict_band_hits": hits,
            "near_boundary": near,
            "high_precision_branch_mismatches": mismatches,
        })
    return {
        "kind": "New post-review targeted diagnostic, not historical benchmark output",
        "date": "2026-10-21",
        "tolerance": tolerance,
        "decimal_precision": 80,
        "boundary_classification_epsilon": "1e-60",
        "scope": "Positive reflections of the three accepted Sigma_1 root sets only",
        "limitations": "Numerical branch check, not an exact certificate or check of higher levels.",
        "cases": cases,
        "positive_candidates_total": sum(c["positive_candidates"] for c in cases),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    result = json.dumps(run_audit(), indent=2)
    if arguments.output:
        if arguments.output.name == "results.json":
            parser.error("Use another filename; the recorded benchmark must remain unchanged.")
        if arguments.output.exists():
            parser.error("The diagnostic output already exists; choose a new filename.")
        arguments.output.write_text(result + "\n", encoding="utf-8")
    else:
        print(result)


if __name__ == "__main__":
    main()
