# Computing Small Roots and Low Elements in Rank-Three Hyperbolic Coxeter Groups

Research code and recorded data accompanying the English manuscript of the same
title, by **Haoyang Tang**. The case studies use the parameter order
`(m12, m13, m23)` and cover `(3,4,5)`, `(4,4,4)`, and `(4,4,infinity)` at
levels zero and one.

The default branch provides an **English-language edition** of the four
notebooks and shared reference source. Chinese comments, messages and saved
text output have been translated; mathematical algorithms and recorded numeric
results are retained. The original publication snapshot is preserved in
[commit fe2c187](https://github.com/MrHaoyang/Computing-Small-Roots-and-Low-Elements-in-Rank-Three-Hyperbolic-Coxeter-Groups/tree/fe2c187aafd91105a6f1330752d1220c1dc7d029),
which is the source snapshot cited by the manuscript.

## Start here

Download this repository using **Code > Download ZIP** and extract it, or clone it.
Open a terminal in the extracted repository root. Use Python 3.12; the recorded
benchmark used Python 3.12.14 and NumPy 2.3.5 on Windows 11.

Install the reference-runner dependency:

```console
python -m pip install -r experiments/requirements.txt
```

Run one small pilot and check that the extracted algebra matches both current
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
| `experiments/reference_algebra.py` | English-language edition of the rank-three reference algebra; its algebraic AST matches the corresponding section in both current visualization notebooks. |
| `experiments/run_experiments.py` | Standalone numerical experiment and audit runner. |
| `experiments/results.json` | Original recorded counts, timings, audits, environment and profiling data. |
| `experiments/audit_thresholds.py` | Separate post-review numerical threshold diagnostic. |
| `experiments/threshold_audit.json` | Recorded output of that separate diagnostic. |
| `notebooks/n small roots.ipynb` | English-language edition of the interactive exploratory prototype. |
| `notebooks/n low elements.ipynb` | English-language edition of the interactive exploratory prototype. |
| `notebooks/shi_arrangement.ipynb` | English-language edition of the zero-Shi arrangement PNG renderer with static low-element input lists. |
| `notebooks/shi_polyhedron.ipynb` | English-language edition of the zero-Shi polyhedron PNG renderer using inverse chambers. |
| `figures/` | Supplied manuscript SVG/PDF assets and their complete-label versions. |
| `MANIFEST.json` | File sizes, SHA-256 hashes and Git blob hashes for this publication snapshot. |

The two general-rank prototypes are archived for provenance. They use bounded
search and numerical heuristics and have known limitations discussed in the
manuscript appendix. Use the rank-three reference runner to reproduce the
paper's tables; the prototypes are not a certified arbitrary-rank implementation.
The four current notebooks and reference source contain the authorized English
translations, so their file hashes differ from the original snapshot.
`experiments/results.json`, its historical measurements and source hashes, and
the supplied vector figures remain byte for byte unchanged. The historical
hashes describe the original snapshot, while `MANIFEST.json` describes this
English-language edition. The runner's numerical logic is unchanged; its
module description now identifies the translation accurately. The original
notebooks remain accessible in the linked historical commit.

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
  Infinity accepts `"inf"`, `"infinity"`, `"unbounded"`, or the infinity symbol.
  The former Chinese-word alias has been replaced by `"unbounded"` in this edition.
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
