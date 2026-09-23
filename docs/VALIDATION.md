# Validation of Matthew's final analysis

Completed 2026-09-23 13:49 UTC using team-input revision `9ad1fb66039963b66a82680c3bf70a9f6c9e7e36`.

- Fresh Jupyter kernel: all four code cells completed on the real uploaded study.
- Real analysis: seven MSA conditions and seven supplied prediction structures across T1188, T1106s1 and T1122.
- Published output: 18 PNG figures and 28 CSV tables, plus report and captions.
- Software tests: **147 passed, 1 skipped**. Synthetic fixtures stay in temporary test folders.
- Optional skip: live equivalence against an independently installed official TMscore executable; this executable is not required for the accepted team-reported scores.
- Original-file preservation: all 57 files present in the team source revision retained their original bytes.
- Published source paths are repository-relative; figure families were visually checked.
- The committed notebook contains only code cells with transient outputs cleared. Run All regenerates the results.

## Environment

Python: `3.10.11 (tags/v3.10.11:7d4cc5a, Apr  5 2023, 00:38:17) [MSC v.1929 64 bit (AMD64)]`

- numpy: 2.2.6
- biopython: 1.88
- matplotlib: 3.10.9
- nbformat: 5.10.4
- nbclient: 0.10.4
- pytest: 8.3.4

## Test output

```text
........................................................................ [ 48%]
......................s................................................. [ 97%]
....                                                                     [100%]
=========================== short test summary info ===========================
SKIPPED [1] tests\test_metrics.py:219: Live official-tool equivalence unavailable: set PROJECT1_TMSCORE to an existing official TMscore executable; none installed by workflow
147 passed, 1 skipped in 15.02s
```

Reproduce with `python -m pytest -q`; see [workflow setup](../MATTHEW_WORKFLOW.md).
