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

Still in progress: README/migration/method documentation and runnable examples;
PEP8 long-line cleanup; final complete pytest/coverage, packaging and rendered
HTML checks. This checkpoint is deliberately not a release.
