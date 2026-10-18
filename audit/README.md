# Bounded exact-rational certificate check

This is a new, data-only checker for two finite claims in Section 2.5 of the
published manuscript. It is **not a complete proof verifier or peer review**.

## Run

From the packaged kit root, with Python 3 (standard library only):

```sh
python3 audit/check_certificate.py --code-dir Code --output-dir audit/results
```

The input `Code/` files are read without modification. No original solver,
certificate generator, replay wrapper, or mathematical implementation is imported
or executed. Output is written only to the specified output directory.

## What is checked

* Three archived JSON files are read, SHA-256 hashed, and their direct recorded
  dependency hash links compared.
* The 17-block grids, positive constants, triangular complete-memory factors and
  several shared reference/defect arrays are checked for consistency.
* All 17 rectangle-invariance inequalities are recomputed using exact fractions.
* A fresh three-step prefix-error majorant is computed without rounding, using
  the supplied unperturbed reference caps and every prior time block.
* All 17 strict third-iterate-entry comparisons are checked.
* All 51 archived prefix recurrence inequalities are separately checked, as are
  the archived entry comparisons and domination of the fresh caps by archived caps.

The final interval is `[3,infinity)`; it is included through the supplied bounds.
Every acceptance comparison uses Python `Fraction`. Decimal output is for display
only. Malformed or inconsistent inputs fail with a nonzero exit status; failed
target inequalities produce `FAIL_DATA_ONLY` and exit status 1.

## Meaning of PASS_DATA_ONLY

It means these finite arithmetic implications hold **assuming the supplied
analytic enclosures are valid**. The checker does not prove that the archived
reference caps, signed Fourier residual bounds, kernel factors or infinite-time
tails actually enclose their mathematical objects. Hash matches are consistency
checks, not a mathematical authenticity or correctness guarantee.

It also does not verify the PDE compactness/continuation argument, the complete
rate/stability certificate, Stage16 or N3, originality, publication suitability,
or any unrestricted Navier–Stokes claim. It is an AI-written separately
implemented finite audit, not independent human mathematical review.

## Result from this run

`PASS_DATA_ONLY`: 17/17 invariance comparisons, 17/17 fresh entry comparisons,
and 51/51 archived recurrence comparisons passed. The smallest invariance margin
was approximately `1.1507455879810338e-19`; the largest fresh `E3/rho` was
approximately `0.9127337584808773`. Exact fractions and source hashes are retained
in `W15_data_only_audit.json`, with per-block results in `W15_blocks.csv`.

The next substantive task is to audit construction of the input enclosures,
starting with the full signed third/fourth heat defects and kernel bounds.
