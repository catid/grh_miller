"""Focused regression tests for the explicit QRH Miller API."""

import math
import unittest

try:
    import gmpy2  # The implementation intentionally uses one GMP backend.
except ModuleNotFoundError:  # pragma: no cover - exercised on bare systems
    gmpy2 = None

if gmpy2 is not None:
    from src.qrh_miller import (
        assumed_qrh_bound,
        assumed_qrh_metadata,
        exhaustive_miller,
        qrh_exact_with_aks_fallback,
    )


@unittest.skipUnless(gmpy2 is not None, "gmpy2 is required by the benchmark implementation")
class QRHMillerTests(unittest.TestCase):
    def test_assumed_bound_requires_explicit_positive_constant(self):
        for bad in (None, 0, -1, float("nan"), float("inf")):
            with self.subTest(C=bad):
                with self.assertRaises(ValueError):
                    assumed_qrh_bound(101, bad)

    def test_bound_uses_requested_log_convention_and_caps(self):
        n = 1009
        natural = assumed_qrh_bound(n, 1, exponent=2, log_base="natural")
        self.assertEqual(natural, math.ceil(math.log(n) ** 2))
        bits = assumed_qrh_bound(n, 1, exponent=2, log_base="bits")
        self.assertEqual(bits, n.bit_length() ** 2)
        # The base range cannot extend to n or beyond.
        self.assertEqual(assumed_qrh_bound(7, 10**9, exponent=1), 5)

    def test_metadata_makes_conditional_status_explicit(self):
        metadata = assumed_qrh_metadata(3, exponent=8, log_base="bits")
        self.assertEqual(metadata["C"], "3")
        self.assertEqual(metadata["exponent"], 8)
        self.assertEqual(metadata["log_base"], "bits")
        self.assertIn("conditional", metadata["proof_status"])
        self.assertIn("not official", metadata["source"])

    def test_known_composite_can_exit_on_prefix_witness(self):
        # 91 is immediately rejected by base 2.  C=1, k=1 gives a prefix
        # containing base 2 and therefore does not need the AKS fallback.
        result = qrh_exact_with_aks_fallback(91, C=1, exponent=1)
        self.assertFalse(result.prime)
        self.assertEqual(result.status, "prefix_witness")
        self.assertEqual(result.proof_status, "unconditional_witness")
        self.assertTrue(result.prefix_used)
        self.assertEqual(result.prefix_result["first_witness"], 2)

    def test_prefix_prechecks_keep_provenance_distinct_from_witnesses(self):
        result = qrh_exact_with_aks_fallback(4, C=1, exponent=1)
        self.assertFalse(result.prime)
        self.assertEqual(result.status, "prefix_precheck")
        self.assertEqual(result.proof_status, "unconditional_precheck")
        self.assertIsNone(result.prefix_result["first_witness"])

    def test_prime_uses_exact_aks_fallback(self):
        result = qrh_exact_with_aks_fallback(97, C=1, exponent=1)
        self.assertTrue(result.prime)
        self.assertEqual(result.status, "aks_fallback")
        self.assertEqual(result.proof_status, "unconditional_aks")
        self.assertIsNotNone(result.aks_result)

    def test_exact_wrapper_agrees_with_trial_division(self):
        cases = [2, 3, 4, 5, 9, 15, 25, 31, 49, 91, 97, 101]
        for n in cases:
            expected = n >= 2 and all(n % d for d in range(2, math.isqrt(n) + 1))
            with self.subTest(n=n):
                result = qrh_exact_with_aks_fallback(n)
                self.assertEqual(result.prime, expected)
                self.assertEqual(result.proof_status, "unconditional_aks")

    def test_exhaustive_diagnostic_is_exact_on_small_inputs(self):
        composite = exhaustive_miller(91)
        self.assertFalse(composite.prime)
        self.assertEqual(composite.cutoff, 89)
        self.assertEqual(composite.first_witness, 2)
        prime = exhaustive_miller(97)
        self.assertTrue(prime.prime)
        self.assertEqual(prime.cutoff, 95)
        self.assertEqual(prime.status, "all_bases_pass")


if __name__ == "__main__":
    unittest.main()
