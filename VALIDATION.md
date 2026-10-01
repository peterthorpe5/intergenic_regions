# Version 1.0.0 validation — 1 October 2026

Base: master c36293fc8f40117a8401a3e42aba255e172c914f.
Local runtime: Linux, Python 3.12.14.

- 195 pytest tests passed, including complete subprocess workflows.
- Combined statement/branch coverage: 98.92% (95% required).
- Ruff lint, 79-column format check and mypy passed.
- Source distribution and pure-Python wheel built successfully.
- Wheel installed into an isolated virtual environment containing only pip,
  packaging, pyfaidx and intergenic-regions. All 20 selected synthetic flanks
  were extracted successfully without numpy/scipy/matplotlib/scikit-learn.
- Complete 40-target, mixed-strand synthetic pipeline ran motif enrichment
  and automatic grouped ML, with 19 permutations and optional evidence.
- Headless Chromium browser checks passed: seven embedded plots, offline
  rendering, table filtering/reset, numerical ascending/descending sorting,
  mobile layout, no JavaScript errors and no external requests.
- The dashboard screenshot was inspected visually.

Synthetic demonstration: 2,081 hypotheses tested; six at q <= 0.05;
held-out ROC AUC 0.975, average precision 0.9735; composition baseline AUC
0.690. Permutation p-value 0.05 is the resolution limit with 19 replicates.
These planted results demonstrate software behaviour, not biological validity.

GitHub Actions is configured for Python 3.11–3.13 plus a minimal-dependency
extraction job; those remote jobs were not run in this environment. Real HOMER
execution/reference downloads use controlled stand-ins in tests. A live HOMER
installation and real biological datasets were not part of local validation.

The GitHub integration rejected writes (HTTP 403); command-line pushing lacked
credentials. The original remote master was not modified. The downloadable
delivery contains a self-contained Git bundle, full source, a patch against the
stated master, distributions, synthetic inputs/results and validation records.
