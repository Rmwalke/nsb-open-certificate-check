"""Portable exact replay of existing N=3 certificate; standard-library only.

Run: python run_replay.py
This reruns the archived algebra audit and the original 224-bit/order-16
finite certificate with its analytic tail; it does not claim a separate
implementation or a new simulation design. Outputs remain in this directory.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parent
SOURCE_NAMES = ["n3_exact_algebra_intermediate.json", "n3_interval_core.py",
                "n3_validated_run.py", "n3_algebra_audit.py", "n3_check_margins.py"]

def hashes():
    return {name: hashlib.sha256((BASE/name).read_bytes()).hexdigest()
            for name in SOURCE_NAMES}

def load(name, archived=False):
    return json.loads((BASE/"archived"/name if archived else BASE/name).read_text())

def main():
    started = datetime.now(timezone.utc).isoformat()
    before = hashes()
    runs = []
    commands = [
        ["n3_algebra_audit.py"],
        ["n3_validated_run.py", "--precision", "224", "--denominator", "128",
         "--order", "16", "--horizon", "8", "--progress", "32"],
        ["n3_check_margins.py"],
    ]
    for args in commands:
        start = time.monotonic()
        log = BASE/(Path(args[0]).stem + ".log")
        with log.open("w") as out:
            result = subprocess.run([sys.executable] + args, cwd=BASE,
                                    stdout=out, stderr=subprocess.STDOUT)
        runs.append({"arguments": args, "exit_code": result.returncode,
                     "elapsed_seconds": time.monotonic()-start, "log": log.name})
        if result.returncode:
            raise RuntimeError(f"Replay failed: {args[0]}; see {log.name}")
    algebra = load("n3_algebra_audit_result.json")
    result = load("n3_validated_result.json")
    archived = load("n3_validated_result.json", True)
    # Wall-clock duration is the only intentionally variable result field.
    comparable = lambda obj: {k:v for k,v in obj.items() if k != "elapsed_seconds"}
    data = load("n3_exact_algebra_intermediate.json")
    support = {tuple(row["k"]) for field in data["basis"] for row in field}
    checks = {
        "source_bytes_unchanged": before == hashes(),
        "algebra_pass": algebra["status"] == "PASS_EXACT_ALGEBRA",
        "algebra_exact_archived_agreement": algebra == load("n3_algebra_audit_result.json", True),
        "all_finite_certificate_fields_exact_archived_agreement_except_elapsed_seconds": comparable(result) == comparable(archived),
        "rational_check_exact_archived_agreement": load("n3_rational_margin_check.json") == load("n3_rational_margin_check.json", True),
        "original_gate_pass": result["status"] == "PASS_ORIGINAL_GATE",
        "completed_horizon_8": result["completed_horizon"] == "8",
        "retained_modes_108": len(support) == 108,
        "support_squared_radius_sum_660": sum(sum(k*k for k in mode) for mode in support) == 660,
    }
    receipt = {"started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "python_version": sys.version, "method": "Reexecution of archived implementation in an isolated copy; exact JSON comparison except wall time",
               "all_checks_pass": all(checks.values()), "checks": checks,
               "source_sha256": before, "runs": runs,
               "precision_bits": result["precision_bits"], "order": result["order"],
               "step": result["h"], "steps": 1024, "saved_records": len(result["records"]),
               "kappa_original_upper": result["kappa_original_upper"],
               "scope": "Direct finite N=3 ODE interval enclosure plus analytic tail. PDE conclusions additionally require the written residual-error estimate, Sobolev bound, all-cutoff comparison, and continuation argument. H1 perturbation ball has distinct initial-energy scope."}
    (BASE/"replay_receipt.json").write_text(json.dumps(receipt, indent=2))
    print(json.dumps({"all_checks_pass": receipt["all_checks_pass"], "checks": checks}, indent=2))
    if not receipt["all_checks_pass"]:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
