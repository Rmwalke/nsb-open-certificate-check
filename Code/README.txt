W15 / W16 / N3: selected proof and computational sources — version 1.1

Version 1.1 clarifies references, the small-data baseline, notation and the
interpretation of quantitative estimates. See CHANGES.txt. The published
version 1.0 is identified by https://doi.org/10.5281/zenodo.23428404.
All scientific sources, archived results, replay tools and replay receipts
are retained unchanged from that package; this revision is editorial.

This package contains the mathematical source dependencies for the manuscript
on global heat-seeded Picard convergence near W=(sin z,sin x,sin y), its
separate enlarged regularity classes, and the supplementary fixed-W N3 bound.
It is a selected source package, not a collection of all the author's research.

Contents
  proof_sources/W_research/  Original mathematical notes, W checker programs
                            and archived comparison certificates.
  checks/replay_W_chain.py   Portable orchestration of the fifteen W checks.
  checks/N3/                Original N3 programs, exact algebra input and runner.
  checks/N3/archived/        Three original N3 comparison outputs.
  SOURCE_MAP.json           Original/public hashes and relative source locations.
  WRAPPER_CHANGES.txt        Exact scope of the W wrapper's packaging changes.
  manuscript/              Consolidated preprint TeX and build instructions.
  checks/W/replay_receipt.json   Retained version 1.0 W replay record.
  checks/N3/replay_receipt.json  Retained version 1.0 N3 replay record.
  LICENSE_SCOPE.txt         Per-material CC BY 4.0 / MIT license mapping.
  CITATION.cff              Author and manuscript identification.
  CHANGES.txt               Version 1.1 changes and unchanged proof material.
  MANIFEST.json            File sizes and SHA-256 digests.
  verify_bundle.py         Integrity check, not a mathematical proof checker.

The scientific programs, their archived result files and all selected original
proof notes are byte-for-byte unchanged. The W wrapper's default paths and
receipt-location metadata were adapted for this package. The source map
identifies the one changed pre-existing program. Generated package documents
are listed separately from original mathematical material.

The original notes retain their historical headings and workflow descriptions,
including the word 'private' in some headers. These describe the source notes
at the time they were written; they are not a new instruction to the reader.
They are preserved because this is source documentation, and one note's bytes
are explicitly hashed into a computational certificate. The consolidated
manuscript is the statement of the present results. Earlier sufficient classes,
bounded failed tests and limits in these notes must be read in that context.

One historical documentation cross-reference is deliberately not included:
endpoint_critical_gate_audit.md names pair_geometry_adversarial_audit.md.
That secondary audit is not a runtime input. The needed complex symmetric
pair bound is proved in the manuscript's "Complex pair bound and heat
estimate" lemma; incompressible_pair_geometry.md is also supplied as an
original related derivation. The omitted audit is not required to run the
certificate chain or to read the consolidated proof.

Related W tests retained
The archived W14 output contains seven declared block tests, W15 contains
eleven radius trials, and W16 contains eleven candidate classes. The successful
certificates and the related unsuccessful tests are retained together so the
original result hashes and checker outputs remain reproducible. A failed
sufficient test does not establish divergence or singularity. These data are
part of the W computation; unrelated research branches are not included.

W reproduction
Use Python 3 with its standard library, with assertions enabled. The version 1.0
release checks used Python 3.12.14 on 16 October 2026 UTC (17 October IST).
All fifteen W outputs were byte-identical to their archived counterparts.
All N3 mathematical data agreed exactly after excluding the wall-clock field.
The accompanying receipts and manuscript record those actual version 1.0
release replays, retained unchanged in version 1.1. No new scientific replay
is claimed for this editorial revision. The documentation makes no verified
Python 3.13 reproduction claim.

From this code-package root:
  python3 checks/replay_W_chain.py

The source defaults to proof_sources/W_research relative to the wrapper's
location. The output defaults to replay_new_W under the caller's current
working directory. The output directory must be absent or empty. Explicit
relative arguments are also supported:
  python3 checks/replay_W_chain.py --source proof_sources/W_research --out replay_new_W

The wrapper creates a working_copy containing the selected programs and
generates result JSONs in dependency order. It compares each with its archived
counterpart. Do not run individual W programs directly in their separate
original-stage directories: their input dependencies are assembled by the
wrapper. It does not edit the source directory. Avoid python -O or -OO.

N3 reproduction
Make a separate copy of checks/N3 before running its writer. In that copy:
  python3 run_replay.py

This runs the original algebra audit, 224-bit/order16 finite interval
validator through time 8 and its analytic tail, then the rational margin
check. It writes fresh outputs and a receipt into that copy. The original
programs use standard-library dependencies and the included n3_interval_core.
The N3 mathematical outputs are compared as parsed data, excluding the
wall-clock field where specified; do not expect whole-file byte equality for
a result that includes a new runtime.

Output and evidence handling
Only the retained version 1.0 release replay receipts are included, not earlier workspace
receipts, execution logs, reviewer programs, raw chat records or review
verdicts. A later replay produces its own
receipt and logs. The W receipt deliberately uses logical public-relative
locations and a replayable command instead of absolute interpreter/workspace
paths. Some unchanged scientific programs print an absolute result-file path
to their raw stdout log; inspect generated raw logs before distributing them.
Those display paths do not enter the archived scientific result JSON bytes.

The original N3 proof notes are in proof_sources/W_research/01_N3_W; the
programs and exact algebra input are placed once in checks/N3. Their original
relative locations are recorded in SOURCE_MAP.json. This avoids duplicate
program/data copies while keeping the reproduction command intact.

Scope of verification
Reexecution of the same programs checks reproducibility of their encoded
certificates. It does not replace the written Fourier, time-history,
Galerkin-to-PDE and continuation arguments. File hashes identify content;
they do not establish correctness, priority or external peer review. The
W16 classes do not automatically inherit W15 Picard convergence, and the
fixed-W N3 bounds do not extend unchanged to the perturbation classes.

Licensing and citation
Read LICENSE_SCOPE.txt: software source (.py) is MIT; the non-software
manuscript, documentation and certificate data are CC BY 4.0, to the extent
applicable rights exist. These are different scopes, not alternative licenses
for each file. The accompanying license notices specify the terms.
CITATION.cff identifies the manuscript and author; scientific citation is
requested and is not an extra restriction on MIT software use.

Integrity
Before generating new outputs, run from this archive's root:
  python3 verify_bundle.py
This checks the file manifest and reports extra files. It does not execute
the mathematical checks. Run scientific replays into the separate copies
described above so the delivered files and receipts remain intact.
