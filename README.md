# Computing Small Roots and Low Elements in Rank-Three Hyperbolic Coxeter Groups

Research code and recorded data accompanying the English manuscript of the same
title, by **Haoyang Tang**. The case studies use the parameter order
`(m12, m13, m23)` and cover `(3,4,5)`, `(4,4,4)`, and `(4,4,infinity)` at
levels zero and one.

## Start here

Download this repository using **Code > Download ZIP** and extract it, or clone it.
Open a terminal in the extracted repository root. Use Python 3.12; the recorded
benchmark used Python 3.12.14 and NumPy 2.3.5 on Windows 11.

Install the reference-runner dependency:

```console
python -m pip install -r experiments/requirements.txt
```

Run one small pilot and check that the extracted algebra matches both original
visualization notebooks:

```console
python experiments/run_experiments.py --pilot --notebooks-dir notebooks
```

The pilot prints one `(3,4,5)`, level-one result: **16 small roots and 43 low
elements**. It does not run the full audits or save the full benchmark output.
Its timings are new measurements, not the paper's historical measurements.

To reproduce all six cases, repeated timings, structural audits, tolerance
comparisons and the separate profile, explicitly choose a new output file:

```console
python experiments/run_experiments.py --repeats 3 --output rerun_results.json --notebooks-dir notebooks
```

Always specify `--output` for a full run: the original runner's default is
`experiments/results.json`. Preserve that file when comparing the paper's
recorded numbers. See [experiments/README.md](experiments/README.md) for the exact
measurement protocol and the distinction between root and element timings.

## Recorded counts

These values come from the preserved [results.json](experiments/results.json).
Level-one counts include the corresponding level-zero sets.

| Parameters | Level | Small roots | Low elements |
| --- | ---: | ---: | ---: |
| (3,4,5) | 0 | 9 | 22 |
| (3,4,5) | 1 | 16 | 43 |
| (4,4,4) | 0 | 9 | 22 |
| (4,4,4) | 1 | 15 | 40 |
| (4,4,infinity) | 0 | 7 | 14 |
| (4,4,infinity) | 1 | 13 | 26 |

## Files and their roles

| Location | Role |
| --- | --- |
| `experiments/reference_algebra.py` | Rank-three algebra used for the paper's tables; extracted without changing the algebraic section shared by the two visualization notebooks. |
| `experiments/run_experiments.py` | Standalone numerical experiment and audit runner. |
| `experiments/results.json` | Original recorded counts, timings, audits, environment and profiling data. |
| `experiments/audit_thresholds.py` | Separate post-review numerical threshold diagnostic. |
| `experiments/threshold_audit.json` | Recorded output of that separate diagnostic. |
| `notebooks/n small roots.ipynb` | Preserved interactive exploratory prototype. |
| `notebooks/n low elements.ipynb` | Preserved interactive exploratory prototype. |
| `notebooks/shi_arrangement.ipynb` | Original zero-Shi arrangement PNG renderer with static low-element input lists. |
| `notebooks/shi_polyhedron.ipynb` | Original zero-Shi polyhedron PNG renderer using inverse chambers. |
| `figures/` | Supplied manuscript SVG/PDF assets and their complete-label versions. |
| `MANIFEST.json` | File sizes, SHA-256 hashes and Git blob hashes for this publication snapshot. |

The two general-rank prototypes are archived for provenance. They use bounded
search and numerical heuristics and have known limitations discussed in the
manuscript appendix. Use the rank-three reference runner to reproduce the
paper's tables; the prototypes are not a certified arbitrary-rank implementation.
All four supplied notebooks, the reference algebra, runner and historical
measurements are preserved byte for byte.

## Open and use the visualization notebooks

Install the optional notebook dependencies, then open JupyterLab:

```console
python -m pip install -r requirements-notebooks.txt
python -m jupyterlab
```

The optional dependency file adds Pillow, SciPy and JupyterLab; their versions
are not recorded as part of the original numerical benchmark environment.

Open a notebook under `notebooks/` and run its code cell. The original
visualization notebooks are each provided as one complete code cell. Their
settings are at the top:

- Change `COXETER_PARAMETERS` using `(m12, m13, m23)`, for example `(4,4,"inf")`.
- Supply the corresponding complete `LOW_0_WORDS`; for the arrangement also
  supply the complete `LOW_1_WORDS`, including every level-zero element.
- `OUTPUT_FILE` selects the PNG filename; `None` uses an automatic name.
- `RENDER_DEPTH`, image-size and label settings control the rendering.

The arrangement uses dark purple for level zero, light purple for level one
outside level zero, and gray for the identity. The polyhedron is the chamber
union indexed by inverses of the supplied level-zero elements. Changing the
parameters also requires changing the input word lists; the plotting cells do
not automatically import a newly computed low-element set from the runner.

The notebooks generate PNG images. The SVG/PDF files in `figures/` are supplied
manuscript assets; executing the original notebooks does not regenerate the
paper's vector export and label placement exactly. See [figures/README.md](figures/README.md).

## Separate threshold diagnostic

```console
python experiments/audit_thresholds.py --output threshold_audit_rerun.json
```

Choose a filename that does not already exist. This checks the three accepted
level-one root sets (123 positive reflection candidates in total), including
an 80-digit Decimal comparison. It does not enumerate low elements or rerun
the timing benchmark. It supplies numerical evidence, not an exact certificate.

## Scope and citation

The study evaluates the stated rank-three systems, levels and tolerances.
Structural checks, finite-ball comparisons and tolerance agreement are not
independent exact-arithmetic proofs. No measured optimization speedup or peak
memory result is claimed. New timings vary by machine and environment.

The manuscript follows the theory in Humphreys, *Reflection Groups and Coxeter
Groups*, and the cited work of Dyer, Hohlweg, Fishel and Mark on small roots,
low elements and Shi arrangements. This repository does not redistribute those
books or papers. Use [CITATION.cff](CITATION.cff) for the repository citation;
the manuscript has no publication DOI recorded here.

No software license has been selected for this repository.
