# Development

Install Python 3.11+ and an editable development environment:

```bash
python -m pip install --editable '.[dev]'
python -m pytest --cov=intergenic_regions --cov-report=term-missing
python -m ruff check .
python -m ruff format --check .
python -m mypy src
python -m build
```

Use PEP8 (79 columns), UK English, type annotations, Google docstrings, logging,
defensive input checks and keyword-only public APIs. Each added public function
needs a meaningful unit/focused workflow test, including malformed inputs.
Do not change the strand/boundary contract without new independent-oracle
regressions. Package-native outputs must not become comma-separated.

Keep optional scientific/ML dependencies out of extraction imports. Model
preprocessing must remain inside folds. Do not tune on reported held-out scores
or label sequence-only predictions as validated enhancers. Maintain the
distinction between missing evidence, no overlap and contradictory evidence.

Use a new output directory for test runs; integration tests run in pytest's
temporary directories. The coverage threshold is 95% including branches;
coverage is evidence of execution, not proof of correct biology. Reference
network calls and HOMER failures are unit-tested with controlled stand-ins.
A real external HOMER deployment is not required by the local/CI test suite.

The GitHub workflow checks Python 3.11–3.13 and builds distributions. A separate
minimal-dependency job verifies extraction without the optional analysis stack.
Do not claim every Python/platform combination was tested locally.

