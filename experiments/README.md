# Reproducing the rank-three case studies

## Files and scope

- `reference_algebra.py`: English-language edition of the algebraic section extracted from the supplied visualization notebook; AST-equivalent to that section in both current notebooks. The original source is preserved in commit `fe2c187aafd91105a6f1330752d1220c1dc7d029`.
- `run_experiments.py`: standalone runner; it does not execute rendering cells, edit original notebooks, or implement an optimized enumerator.
- `results.json`: recorded results underlying this manuscript's tables and experimental statements.
- `requirements.txt`: the NumPy version used for the recorded run. The recorded Python version is 3.12.14.

The script uses only NumPy and the Python standard library. The current notebooks are optional for source-integrity checks, not required to reproduce enumeration. Their comments, messages and saved text outputs have been translated to English without changing mathematical algorithms or saved numeric results. Historical hashes in `results.json` refer to the original snapshot; current English-edition hashes are in the repository's `MANIFEST.json`.

## Run

Open a terminal in this directory, using a Python 3.12 environment:

```text
python -m pip install -r requirements.txt
python run_experiments.py --repeats 3 --output rerun_results.json
```

Keep `results.json` unchanged when checking the manuscript's recorded numbers. A rerun generates new timing measurements. To additionally verify the extraction against the two current English-language visualization notebooks, supply the directory containing them:

```text
python run_experiments.py --repeats 3 --output rerun_results.json --notebooks-dir "PATH_TO_NOTEBOOK_DIRECTORY"
```

Both `shi_arrangement.ipynb` and `shi_polyhedron.ipynb` must be present for that optional check. The script compares algebraic ASTs and checks SHA-256 hashes before and after execution. A quick extraction-and-enumeration check is available with `--pilot`; it prints one (3,4,5), n=1 case and does not write the full output.

## Protocol and interpretation

The full run has six untimed baseline enumeration/audit passes, eighteen timed repetitions, twelve tolerance-variant runs, and one separate profiling run. Each timed repetition uses a fresh system, then times root enumeration followed by low-element enumeration. The low routine internally enumerates roots again but starts with the root registry populated by the preceding call. These measurements are not a cold-start end-to-end pipeline benchmark. Constructors, imports, audits, rendering, and output I/O are excluded. Median times are reported; individual measurements and ranges are retained.

The tested tolerances are 1e-10, 1e-9 and 1e-8. Cross-instance comparisons use full vectors and matrices with a fixed normalized maximum-entry residual bound of 1e-10; numerical identifiers are not compared between instances. Baseline audits include root norms, positivity, inversion length and sign, suffix closure, braid/cancellation checks, layer nesting, and a full ball through maximum accepted length plus one using the same cone predicate. The supplied-word set-equality comparison applies only to the (3,4,5) lists.

The output records processor, usable memory, Python/NumPy/OS versions and NumPy build information. One-thread BLAS/OpenMP settings are requested before NumPy import, not asserted as verified backend thread counts. If threadpoolctl is already installed, its read-only backend information is also included; it is not required.

This is numerical validation of a bounded case study, not an exact-arithmetic certificate, a proof of completeness from the finite ball alone, a memory benchmark, or a measured optimization speedup. The profile has instrumentation overhead and its cumulative method times overlap. It must not be used as a fourth timing repetition or summed across nested functions. No general-rank prototype notebook was repaired in these experiments.

## Separate post-review threshold diagnostic (2026-10-21)

The historical results.json did not record individual reflection pairings.
The new audit_thresholds.py therefore performs a separate narrowly scoped
check; it does not revise that historical file or rerun low-element timings.

Run:

    python audit_thresholds.py --output threshold_audit_rerun.json

Choose a new output filename. The script refuses results.json and refuses to
overwrite an existing diagnostic. The supplied threshold_audit.json records
the post-review check, not the original benchmark.

Only the three level-one root sets are enumerated. They contain 16, 15 and
13 roots, giving 45, 42 and 36 positive reflection candidates after the
original sign cutoff, 123 in total. No binary64 pairing lies in either
strict discrepancy band [1-1e-7,1) or (-1,-1+1e-7].

For an additional numerical check, each accepted root is reconstructed
from a simple root and a recovered reflection path using the standard
library Decimal at 80-digit precision. Gram entries are rebuilt from
1/2, sqrt(2)/2, (1+sqrt(5))/4 and 1, rather than by converting binary64
root vectors. Values within 1e-60 of either exact boundary are classified
as numerical boundary matches. The six near-boundary candidates occur
in (4,4,infinity) at b=+1 or b=-1; all 123 branch classifications agree.

This remains high-precision numerical evidence, not interval or exact
certification. It does not cover roots outside the accepted level-one
sets or larger levels. Both scripts use the same unchanged reference
algebra and NumPy dependency.
