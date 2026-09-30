#!/usr/bin/env python3
"""Recompute the Appendix E diagnostics in a fresh directory and check its table.

Uses only Python's standard library to orchestrate the separately supplied
mpmath calculation. Numerical agreement is not an interval certificate.
"""
from __future__ import annotations
import argparse
import concurrent.futures
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_job(work: Path, dps: int, timeout: int) -> dict:
    name = f"screw_kernel_{dps}.json"
    command = [sys.executable, "-B", "screw_kernel_diagnostics.py", "--dps", str(dps), "--output", name]
    start = time.monotonic()
    process = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=timeout)
    (work / f"run_{dps}.log").write_text(process.stdout + process.stderr, encoding="utf-8")
    if process.returncode != 0:
        raise RuntimeError(f"Diagnostic at {dps} digits failed: {process.stderr}")
    data = json.loads((work / name).read_text(encoding="utf-8"))
    if data.get("status") != "PASS" or not data["exact_checks"]["all_pass"]:
        raise AssertionError(f"Failed numerical or exact check at {dps} digits")
    return {"dps": dps, "command": command, "exit_code": process.returncode,
            "elapsed_seconds": round(time.monotonic() - start, 3),
            "output": name, "sha256": sha(work / name), "data": data}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="New, nonexistent directory")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    root = Path(__file__).resolve().parent
    output = args.output_dir.resolve()
    if output.exists():
        parser.error("output directory already exists; use a fresh path")
    work = output / "work"
    work.mkdir(parents=True)
    for name in ("screw_kernel_diagnostics.py", "requirements.txt"):
        shutil.copy2(root / name, work / name)
    # No expected JSON, source PDFs, zero tables, or manuscript are copied into work.
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        jobs = list(pool.map(lambda d: run_job(work, d, args.timeout), (40, 60)))
    source = (root / "main.tex").read_text(encoding="utf-8")
    start = source.index(r"\label{tab:screw_diagnostics}")
    table = source[start:source.index(r"\end{table}", start)]
    printed = re.findall(r"\$(1/2|1|2)\$\s*&\s*\$(\d+)\$\s*&\s*\$(0\.\d+)\$", table)
    if len(printed) != 3:
        raise AssertionError("Could not extract all three printed Appendix E rows")
    comparisons = []
    for job in jobs:
        data = job.pop("data")
        reference = root / job["output"]
        job["byte_identical_to_packaged_output"] = sha(reference) == job["sha256"]
        expected = json.loads(reference.read_text(encoding="utf-8"))
        # Python version is environmental, not mathematical output.
        comparable = {k: v for k, v in data.items() if k != "python"}
        expected_comparable = {k: v for k, v in expected.items() if k != "python"}
        job["json_equal_except_python_version"] = comparable == expected_comparable
        if not job["json_equal_except_python_version"]:
            raise AssertionError("Regenerated output differs beyond Python version metadata")
        job["exact_checks"] = data["exact_checks"]["count"]
        for (a, cap, value), window in zip(printed, data["arithmetic_windows"]):
            observed = Decimal(window["eigenvalues"][0])
            target = Decimal(value)
            tolerance = Decimal(5) * Decimal(10) ** (target.as_tuple().exponent - 1)
            expected_a = Decimal("0.5") if a == "1/2" else Decimal(a)
            passed = (abs(observed - target) <= tolerance
                      and int(cap) == window["complete_integer_cutoff"]
                      and expected_a == Decimal(window["a"]))
            comparisons.append({"dps": job["dps"], "a": a, "cutoff": int(cap),
                                "printed": value, "recomputed": str(observed),
                                "rounding_tolerance": str(tolerance), "pass": passed})
            if not passed:
                raise AssertionError(f"Printed row failed at a={a}, dps={job['dps']}")
    report = {"status": "PASS", "manuscript_sha256": sha(root / "main.tex"),
              "driver_sha256": sha(Path(__file__)), "jobs": jobs,
              "printed_comparisons": comparisons, "printed_comparison_count": len(comparisons),
              "distinct_exact_checks": jobs[0]["exact_checks"],
              "scope": "Fresh-directory regression and exact rational checks; no interval certification, global positivity, RH proof, or source-derived arithmetic operator."}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "jobs": len(jobs), "printed_comparisons": len(comparisons),
                      "distinct_exact_checks": jobs[0]["exact_checks"]}))


if __name__ == "__main__":
    main()
