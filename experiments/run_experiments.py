"""Reproduce the manuscript's bounded rank-three numerical experiments.

The reference algebra is an unchanged extraction from shi_arrangement.ipynb.
This runner never executes rendering cells or edits the supplied notebooks.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import contextlib
import cProfile
import ctypes
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import pstats
import statistics
import sys
import time

# Set before importing NumPy. These requests affect this process only.
for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from reference_algebra import CoxeterSystem3


HERE = Path(__file__).resolve().parent
PARAMETERS = ((3, 4, 5), (4, 4, 4), (4, 4, "inf"))
EXPECTED = {
    (3, 4, 5): ((9, 22), (16, 43)),
    (4, 4, 4): ((9, 22), (15, 40)),
    (4, 4, "inf"): ((7, 14), (13, 26)),
}
SUPPLIED_345_0 = (
    "e", "1", "2", "3", "12", "13", "21", "23", "31", "32",
    "121", "131", "232", "313", "323", "2323", "3121", "3131",
    "3232", "23131", "23232", "123232",
)
SUPPLIED_345_1 = SUPPLIED_345_0 + (
    "123", "132", "213", "231", "312", "321", "1213", "1232",
    "1323", "2131", "2313", "12313", "12323", "13121", "13232",
    "23121", "123131", "213232", "323131", "1213232", "3123232",
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def notebook_algebra(path):
    notebook = json.loads(Path(path).read_text(encoding="utf-8"))
    code = "\n".join(
        "".join(cell["source"]) for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
    start = code.index('"""Algebraic source of truth')
    stop = code.index('"""Analytic, paper-style rendering')
    return ast.dump(ast.parse(code[start:stop]), include_attributes=False)


def environment():
    information = {
        "python": platform.python_version(), "numpy": np.__version__,
        "os": platform.platform(), "machine": platform.machine(),
        "logical_processors": os.cpu_count(),
        "blas_thread_requests": {name: os.environ[name] for name in
                                 ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS")},
    }
    if sys.platform == "win32":
        import winreg
        with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                            r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as key:
            information["cpu"] = winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()

        class MemoryStatus(ctypes.Structure):
            _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong)] + [
                (name, ctypes.c_ulonglong) for name in (
                    "total_physical", "available_physical", "total_pagefile",
                    "available_pagefile", "total_virtual", "available_virtual",
                    "available_extended_virtual",
                )
            ]

        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
            information["physical_memory_gib"] = status.total_physical / 2**30
    with contextlib.redirect_stdout(io.StringIO()) as buffer:
        np.show_config()
    information["numpy_build"] = buffer.getvalue()
    try:
        from threadpoolctl import threadpool_info
        information["loaded_threadpools"] = threadpool_info()
    except ImportError:
        information["loaded_threadpools"] = "Not independently inspected; threadpoolctl unavailable."
    return information


def relative_residual(left, right):
    return float(np.max(np.abs(left - right)) /
                 (1 + max(float(np.max(np.abs(left))), float(np.max(np.abs(right))))))


def sets_match(left, right, bound=1e-10):
    """Compare full vectors/matrices, not instance-dependent numerical keys."""
    if len(left) != len(right):
        return False, None
    unused = set(range(len(right)))
    maximum = 0.0
    for item in left:
        matches = [(relative_residual(item, right[index]), index) for index in unused]
        error, index = min(matches)
        if error > bound:
            return False, error
        maximum = max(maximum, error)
        unused.remove(index)
    return not unused, maximum


def enumerate_case(parameters, level, tolerance):
    system = CoxeterSystem3(*parameters, tolerance=tolerance)
    begin = time.perf_counter()
    small = system.enumerate_m_small_roots(level)
    root_seconds = time.perf_counter() - begin
    begin = time.perf_counter()
    low = system.enumerate_m_low_elements(level)
    low_seconds = time.perf_counter() - begin
    assert (len(small), len(low)) == EXPECTED[parameters][level]
    return system, small, low, {
        "parameters": list(parameters), "level": level, "tolerance": tolerance,
        "small_roots": len(small), "low_elements": len(low),
        "root_seconds": root_seconds, "low_seconds": low_seconds,
        "maximum_word_length": max(record.length for record in low.values()),
        "length_histogram": dict(sorted(Counter(record.length for record in low.values()).items())),
        "max_matched_root_residual": system.max_root_residual,
        "max_matched_element_residual": system.max_element_residual,
        "completed_normally": True,
    }


def audit_case(system, small, low, parameters, level):
    """Invariants and bounded cross-checks; not an exact arithmetic certificate."""
    maximum_norm_error = 0.0
    maximum_negative_sign_error = 0.0
    inversions_checked = 0
    for root in small.values():
        beta = root.coefficients
        assert np.all(beta >= -1e-10) and np.max(beta) > 0
        error = abs(float(beta @ system.gram @ beta) - 1.0)
        maximum_norm_error = max(maximum_norm_error, error)
        assert error <= 1e-8
    for record in low.values():
        assert len(record.inversions) == record.length
        assert len(set(record.inversions)) == record.length
        inverse = np.linalg.inv(record.matrix)
        for key in record.inversions:
            beta = system._root_registry[key].coefficients
            assert np.all(beta >= -1e-10) and np.max(beta) > 0
            transformed = inverse @ beta
            scale = 1.0 + np.max(np.abs(transformed))
            error = max(0.0, float(np.max(transformed))) / scale
            maximum_negative_sign_error = max(maximum_negative_sign_error, error)
            assert error <= 1e-8
            inversions_checked += 1
        if record.length:
            suffix = record.word[1:] or "e"
            assert system.element_key(system.word_matrix(suffix)) in low
    for i, j, order in (("1", "2", parameters[0]), ("1", "3", parameters[1]),
                         ("2", "3", parameters[2])):
        if isinstance(order, int):
            product = system.generator_matrices[i] @ system.generator_matrices[j]
            assert np.allclose(np.linalg.matrix_power(product, order), np.eye(3), atol=1e-10)
            left = (i + j) * order
            right = (j + i) * order
            assert system.same_element(left[:order], right[:order])
    assert system.same_element("11", "e")
    if parameters[0] == 3:
        assert system.canonicalize_word("1212").length == 2
        assert system.same_element("1212", "21")
    # One ball larger than the maximum output length provides a bounded
    # alternative to the suffix frontier. It shares the reference predicates.
    bound = max(record.length for record in low.values()) + 1
    ball = system.enumerate_chambers(bound)
    passing_keys = {key for key, record in ball.items() if system.is_low_element(record, small)}
    assert passing_keys == set(low)
    supplied_match = None
    if parameters == (3, 4, 5):
        words = (SUPPLIED_345_0, SUPPLIED_345_1)[level]
        supplied_keys = {system.element_key(system.word_matrix(word)) for word in words}
        supplied_match = supplied_keys == set(low)
        assert supplied_match
    return {
        "inversions_checked": inversions_checked,
        "maximum_root_norm_error": maximum_norm_error,
        "maximum_negative_sign_error": maximum_negative_sign_error,
        "bounded_ball_length": bound, "bounded_ball_size": len(ball),
        "bounded_ball_set_equality": True, "suffix_closure": True,
        "braid_and_cancellation": True, "supplied_345_set_equality": supplied_match,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--pilot", action="store_true", help="Run one case without full audits.")
    parser.add_argument("--notebooks-dir", type=Path,
                        help="Optionally verify the extraction against both original notebooks.")
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be at least one")
    result = {
        "scope": "Numerical rank-three, n=0,1; no exact certificate or optimization speedup claim.",
        "environment": environment(), "reference_sha256": sha256(HERE / "reference_algebra.py"),
        "protocol": {"repeats": args.repeats, "warmup": "One untimed baseline enumeration per case.",
                     "tolerances": [1e-10, 1e-9, 1e-8], "comparison_bound": 1e-10,
                     "root_timer": "Root enumeration only; constructor excluded.",
                     "low_timer": "Low enumeration including its internal root enumeration; constructor excluded.",
                     "cache_policy": "Fresh system per repetition; root timer precedes low timer on that system.",
                     "excluded": "Audits, imports, rendering, source extraction and output I/O."},
        "cases": [], "sensitivity": [],
    }
    if args.notebooks_dir:
        before = {}
        extracted = ast.dump(ast.parse((HERE / "reference_algebra.py").read_text(encoding="utf-8")),
                             include_attributes=False)
        for name in ("shi_arrangement.ipynb", "shi_polyhedron.ipynb"):
            path = args.notebooks_dir / name
            before[name] = sha256(path)
            assert notebook_algebra(path) == extracted, name
        result["original_notebooks"] = before
        result["algebra_ast_matches_both_notebooks"] = True
    if args.pilot:
        _, _, _, entry = enumerate_case((3, 4, 5), 1, 1e-9)
        print(json.dumps(entry), flush=True)
        return
    snapshots = {}
    for parameters in PARAMETERS:
        for level in (0, 1):
            system, small, low, _ = enumerate_case(parameters, level, 1e-9)
            audit = audit_case(system, small, low, parameters, level)
            snapshots[(parameters, level)] = (
                [record.coefficients.copy() for record in small.values()],
                [record.matrix.copy() for record in low.values()],
            )
            if level == 1:
                previous_roots, previous_elements = snapshots[(parameters, 0)]
                assert {system.root_key(root) for root in previous_roots} <= set(small)
                assert {system.element_key(matrix) for matrix in previous_elements} <= set(low)
                audit["layer_nesting"] = True
            runs = []
            for repeat in range(args.repeats):
                _, _, _, run = enumerate_case(parameters, level, 1e-9)
                runs.append(run)
                print(f"baseline {parameters} n={level} repetition {repeat+1}: "
                      f"roots={run['small_roots']} low={run['low_elements']} "
                      f"root={run['root_seconds']:.6f}s low={run['low_seconds']:.6f}s", flush=True)
            entry = {"parameters": list(parameters), "level": level, "audit": audit, "runs": runs}
            for timer in ("root_seconds", "low_seconds"):
                values = [run[timer] for run in runs]
                entry[timer] = {"median": statistics.median(values), "minimum": min(values),
                                "maximum": max(values)}
            result["cases"].append(entry)
    for tolerance in (1e-10, 1e-8):
        for parameters in PARAMETERS:
            for level in (0, 1):
                system, small, low, entry = enumerate_case(parameters, level, tolerance)
                baseline_roots, baseline_elements = snapshots[(parameters, level)]
                root_match, root_error = sets_match(
                    baseline_roots, [record.coefficients for record in small.values()])
                element_match, element_error = sets_match(
                    baseline_elements, [record.matrix for record in low.values()])
                assert root_match and element_match
                entry.update({"root_set_matches_baseline": root_match,
                              "element_set_matches_baseline": element_match,
                              "maximum_comparison_root_residual": root_error,
                              "maximum_comparison_element_residual": element_error})
                result["sensitivity"].append(entry)
                print(f"sensitivity tau={tolerance:g} {parameters} n={level}: full sets matched", flush=True)
    profile = cProfile.Profile()
    profile.enable()
    profiling_system = CoxeterSystem3(3, 4, 5)
    profile_begin = time.perf_counter()
    profiling_low = profiling_system.enumerate_m_low_elements(1)
    profile_total = time.perf_counter() - profile_begin
    profile.disable()
    assert len(profiling_low) == 43
    stats = pstats.Stats(profile)
    result["profile"] = {"case": [3, 4, 5], "level": 1,
                         "profiled_elapsed_seconds": profile_total, "methods": {}}
    for (_, _, name), (primitive, calls, own, cumulative, _) in stats.stats.items():
        if name in {"canonicalize_word", "enumerate_chambers", "cone_contains",
                    "is_low_element", "enumerate_m_small_roots"}:
            result["profile"]["methods"][name] = {
                "primitive_calls": primitive, "total_calls": calls,
                "self_seconds": own, "cumulative_seconds": cumulative,
            }
    if args.notebooks_dir:
        assert all(sha256(args.notebooks_dir / name) == digest for name, digest in before.items())
        result["original_notebooks_unchanged"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Completed: {args.output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
