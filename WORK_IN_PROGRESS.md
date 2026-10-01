# Recovery checkpoint — 1 October 2026

This branch preserves the recovered overhaul after an interrupted chat.
It is based on master commit c36293fc8f40117a8401a3e42aba255e172c914f.

Implemented: indexed strand-aware extraction; full annotation blockers;
GFF3/GTF/legacy coordinate parsing; native motif enrichment; held-out sequence
learning; optional accessibility/functional/reference evidence; inert JSON
models; HOMER hand-off; CLI; atomic result bundles and a unified offline HTML
dashboard. Motif workflows now attempt ML by default. Candidate priorities
retain the explicit `unvalidated_candidate` interpretation.

Verification before the latest additions: 153 tests passed, 98.12% combined
statement/branch coverage and mypy passed. The focused suite for automatic ML,
prioritisation, annotation and extraction subsequently passed 70 tests.

Latest verification: 195 pytest tests passed; combined statement/branch
coverage is 98.92%. Ruff lint/format, mypy and distribution builds passed.
A clean virtual environment containing only the wheel and its base dependencies
successfully extracted all 20 positive demo flanks, without the optional ML or
plotting stack. Full documentation, CI configuration and a reproducible
mixed-strand/blocker example are now present.

Still in progress: final rendered HTML/browser verification and complete
delivery archive. GitHub integration rejected writes with HTTP 403; source and
Git-history checkpoints have therefore been saved as downloadable artifacts.
