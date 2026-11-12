#!/usr/bin/env python3
"""Small reproducible comparison run for the README example."""
from __future__ import annotations
import argparse, csv, os, statistics, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.aks import AKSTimeout, aks
from src.baselines import fixed_miller_rabin, gmpy2_probable_prime
from src.miller import Cutoff, strong_miller

CASES = {
    "64-bit balanced semiprime": 13373461803509102023,
    "64-bit certified prime": 17548706721817263871,
}
ALGORITHMS = ("miller_grh_reference", "aks_textbook", "fixed_miller_rabin", "gmpy2_probable_prime")


def timed(name, n, repeats, timeout):
    values, verdict, proof = [], None, ""
    for _ in range(repeats):
        start = time.perf_counter()
        if name == "miller_grh_reference":
            result = strong_miller(n, Cutoff("miller_grh_reference"), timeout=timeout)
            verdict, proof = result.prime, "conditional_ERH"
        elif name == "aks_textbook":
            try:
                result = aks(n, timeout=timeout)
            except AKSTimeout:
                return None, None, "timeout", "unconditional_exact_algorithm"
            verdict, proof = result.verdict, "unconditional_exact_algorithm"
        elif name == "fixed_miller_rabin":
            result = fixed_miller_rabin(n)
            verdict, proof = result["prime"], "probable_fixed_bases"
        elif name == "gmpy2_probable_prime":
            result = gmpy2_probable_prime(n)
            verdict, proof = result["prime"], "probable_GMP"
        values.append(time.perf_counter() - start)
    return statistics.median(values), max(values), "completed", proof


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=2.0)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "quick_benchmark.csv")
    args = parser.parse_args()
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    rows = []
    for case, n in CASES.items():
        for algorithm in ALGORITHMS:
            median, maximum, status, proof = timed(algorithm, n, args.repeats, args.timeout)
            rows.append({"case": case, "bits": n.bit_length(), "algorithm": algorithm,
                         "median_seconds": median, "max_seconds": maximum,
                         "status": status, "proof_status": proof})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    print(f"wrote {args.output}")


if __name__ == "__main__": main()
