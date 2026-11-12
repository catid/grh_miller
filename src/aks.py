"""A small, transparent reference implementation of textbook AKS.

This module implements the algorithm of Agrawal--Kayal--Saxena (2004), not
the substantially faster Lenstra--Pomerance family.  Polynomial arithmetic is
over ``Z/nZ[X]/(X**r - 1)`` and uses GMP Kronecker multiplication. It is deliberately
an audit/reference implementation; it is not intended to compete with a
production primality prover at thousand-bit sizes.

The implementation uses ``ceil(log2(n))`` in the order threshold and in the
bound on ``a``.  This is a conservative integer rounding of the textbook
real-valued quantities and can test a few more cases than the minimum
specification.  The conservative rounding does not weaken correctness.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import time
from typing import Optional

try:  # gmpy2 is used when present, with Python integers as a fallback.
    import gmpy2
except ImportError:  # pragma: no cover - exercised only outside the benchmark env
    gmpy2 = None


class AKSTimeout(TimeoutError):
    """The per-case AKS timeout expired."""

    def __init__(self, message: str, stats: Optional["AKSStats"] = None):
        super().__init__(message)
        self.stats = stats


@dataclass
class AKSStats:
    verdict: Optional[bool] = None
    bits: int = 0
    r: int = 0
    order: int = 0
    order_target: int = 0
    a_max: int = 0
    bases_tested: int = 0
    gcd_bases_tested: int = 0
    gcd_shortcuts: int = 0
    perfect_power: bool = False
    polynomial_tests: int = 0
    polynomial_multiplications: int = 0
    coefficient_multiplications: int = 0
    packed_integer_multiplications: int = 0
    setup_seconds: float = 0.0
    parameter_seconds: float = 0.0
    gcd_seconds: float = 0.0
    polynomial_seconds: float = 0.0
    total_seconds: float = 0.0
    stage: str = ""
    interrupted_stage: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _gcd(a: int, b: int) -> int:
    if gmpy2 is not None:
        return int(gmpy2.gcd(a, b))
    return math.gcd(a, b)


def _isqrt(n: int) -> int:
    if gmpy2 is not None:
        return int(gmpy2.isqrt(n))
    return math.isqrt(n)


def _iroot(n: int, k: int) -> tuple[int, bool]:
    if gmpy2 is not None:
        root, exact = gmpy2.iroot(n, k)
        return int(root), bool(exact)
    # Binary search avoids floating-point errors for large n.
    lo, hi = 1, 1 << ((n.bit_length() + k - 1) // k)
    while lo <= hi:
        mid = (lo + hi) // 2
        v = mid**k
        if v == n:
            return mid, True
        if v < n:
            lo = mid + 1
        else:
            hi = mid - 1
    return hi, False


def is_perfect_power(n: int) -> bool:
    """Return whether n is a nontrivial integer perfect power."""
    if n < 4:
        return False
    if gmpy2 is not None:
        return bool(gmpy2.is_power(n))
    # It suffices to try prime exponents: if n=x**k, a prime divisor p of k
    # gives n=(x**(k/p))**p.  This matters for the 1K--8K-bit corpus.
    limit = n.bit_length()
    prime_exponents = []
    for k in range(2, limit + 1):
        if all(k % p for p in prime_exponents if p * p <= k):
            prime_exponents.append(k)
    for k in prime_exponents:
        _, exact = _iroot(n, k)
        if exact:
            return True
    return False


def _phi(r: int) -> int:
    result = r
    x = r
    p = 2
    while p * p <= x:
        if x % p == 0:
            while x % p == 0:
                x //= p
            result -= result // p
        p += 1
    if x > 1:
        result -= result // x
    return result


def _check_timeout(deadline: Optional[float]) -> None:
    if deadline is not None and time.perf_counter() >= deadline:
        raise AKSTimeout("AKS per-case timeout expired")


def _order_mod(n: int, r: int, target: int, deadline: Optional[float]) -> int:
    """Return ord_r(n), or target+1 when it is known to exceed target."""
    if _gcd(n, r) != 1:
        return 0
    x = n % r
    for k in range(1, target + 1):
        _check_timeout(deadline)
        if x == 1:
            return k
        x = (x * n) % r
    return target + 1


def _find_r(n: int, target: int, deadline: Optional[float], stats: AKSStats) -> tuple[int, int]:
    # AKS chooses the first r with ord_r(n) > log2(n)^2.
    r = 2
    while True:
        _check_timeout(deadline)
        order = _order_mod(n, r, target, deadline)
        if order > target:
            return r, order
        r += 1


def _poly_mul_schoolbook(a: list[int], b: list[int], n: int) -> list[int]:
    """Independent, deliberately simple cyclic multiplication for validation."""
    r = len(a)
    out = [0] * r
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[(i + j) % r] += int(ai) * int(bj)
    return [x % n for x in out]


def _poly_mul(a: list[int], b: list[int], n: int, stats: AKSStats,
              deadline: Optional[float]) -> list[int]:
    """Exact cyclic product using non-overlapping Kronecker digits.

    For canonical nonnegative input coefficients, every ordinary product
    coefficient is at most r*(n-1)**2. A base exceeding that number ensures
    no carry crosses coefficients. Byte-aligned slots permit O(r) packing
    through Python byte strings; the large multiplication itself uses GMP.
    The ordinary coefficients are then folded cyclically and reduced mod n.
    """
    _check_timeout(deadline)
    r = len(a)
    if not r or len(b) != r:
        raise ValueError("polynomials must have the same positive length")
    width = max(1, (r * (n - 1) ** 2).bit_length())
    size = (width + 7) // 8
    packed_a = int.from_bytes(b"".join(int(x % n).to_bytes(size, "little") for x in a), "little")
    packed_b = int.from_bytes(b"".join(int(x % n).to_bytes(size, "little") for x in b), "little")
    if gmpy2 is not None:
        product = gmpy2.mpz(packed_a) * gmpy2.mpz(packed_b)
        stats.packed_integer_multiplications += 1
    else:  # outside the pinned benchmark environment only
        product = packed_a * packed_b
    stats.polynomial_multiplications += 1
    # Do not claim r*r scalar modular multiplications: this backend performs
    # one packed product, and GMP's internal multiplication count is opaque.
    data = int(product).to_bytes((2 * r - 1) * size, "little")
    out = [0] * r
    for j in range(2 * r - 1):
        coefficient = int.from_bytes(data[j * size:(j + 1) * size], "little")
        out[j % r] += coefficient
    out = [x % n for x in out]
    _check_timeout(deadline)
    return out


def _poly_pow_x_plus_a(n: int, r: int, a: int, exponent: int,
                       stats: AKSStats, deadline: Optional[float]) -> list[int]:
    base = [0] * r
    base[0] = a % n
    if r == 1:
        base[0] = (base[0] + 1) % n
    else:
        base[1] = 1
    out = [0] * r
    out[0] = 1 % n
    e = exponent
    while e:
        _check_timeout(deadline)
        if e & 1:
            out = _poly_mul(out, base, n, stats, deadline)
        e >>= 1
        if e:
            base = _poly_mul(base, base, n, stats, deadline)
    return out


def _poly_identity_rhs(n: int, r: int, a: int) -> list[int]:
    rhs = [0] * r
    rhs[0] = a % n
    rhs[n % r] = (rhs[n % r] + 1) % n
    return rhs


def aks(n: int, timeout: Optional[float] = None) -> AKSStats:
    """Run textbook AKS and return a detailed result.

    ``timeout`` is a wall-clock budget in seconds for this input.  A timeout
    raises :class:`AKSTimeout`; callers should record it as censored data.
    """
    n = int(n)
    start = time.perf_counter()
    deadline = None if timeout is None else start + float(timeout)
    stats = AKSStats(bits=max(0, n.bit_length()))
    active_phase = "setup"
    phase_start = start
    try:
        stats.stage = "boundary"
        if n in (2, 3):
            stats.verdict = True
            return stats
        if n < 2:
            stats.verdict = False
            return stats
        if is_perfect_power(n):
            stats.perfect_power = True
            stats.verdict = False
            return stats
        if n % 2 == 0:
            stats.verdict = False
            return stats

        setup_end = time.perf_counter()
        stats.setup_seconds = setup_end - start
        stats.stage = "parameter"
        active_phase = "parameter"
        phase_start = time.perf_counter()
        # Exact ceil(log2(n)); a floating conversion can round down close to
        # a power of two or overflow on the largest input sizes.
        L = max(1, (n - 1).bit_length())
        target = L * L
        stats.order_target = target
        pstart = time.perf_counter()
        r, order = _find_r(n, target, deadline, stats)
        stats.r, stats.order = r, order
        phi = _phi(r)
        # ceil(sqrt(phi(r))) * ceil(log2(n)) is a conservative integer bound.
        a_max = _isqrt(phi)
        if a_max * a_max < phi:
            a_max += 1
        a_max *= L
        stats.a_max = a_max
        stats.parameter_seconds = time.perf_counter() - pstart

        # Textbook AKS checks gcd(a,n) for 2 <= a <= r before the polynomial
        # congruences.  This is required for composites with a small factor
        # that need not lie in the polynomial base range.
        stats.stage = "gcd"
        active_phase = "gcd"
        phase_start = time.perf_counter()
        for a in range(2, r + 1):
            _check_timeout(deadline)
            stats.gcd_bases_tested += 1
            g = _gcd(a, n)
            if 1 < g < n:
                stats.gcd_shortcuts += 1
                stats.verdict = False
                return stats

        if n <= r:
            stats.stage = "n_le_r"
            stats.verdict = True
            return stats

        stats.gcd_seconds = time.perf_counter() - phase_start

        stats.stage = "polynomial"
        active_phase = "polynomial"
        phase_start = time.perf_counter()
        pstart = time.perf_counter()
        for a in range(1, a_max + 1):
            _check_timeout(deadline)
            g = _gcd(a, n)
            if 1 < g < n:
                stats.gcd_shortcuts += 1
                stats.verdict = False
                return stats
            stats.bases_tested += 1
            lhs = _poly_pow_x_plus_a(n, r, a, n, stats, deadline)
            rhs = _poly_identity_rhs(n, r, a)
            stats.polynomial_tests += 1
            if lhs != rhs:
                stats.verdict = False
                return stats
        stats.polynomial_seconds = time.perf_counter() - pstart
        stats.stage = "final"
        # Once all congruences have passed, textbook AKS accepts.
        stats.verdict = True
        return stats
    except AKSTimeout as exc:
        stats.interrupted_stage = stats.stage
        stats.stage = "timeout"
        exc.stats = stats
        raise
    finally:
        finish = time.perf_counter()
        setattr(stats, active_phase + "_seconds", finish - phase_start)
        stats.total_seconds = finish - start


def is_prime_aks(n: int, timeout: Optional[float] = None) -> bool:
    """Boolean convenience wrapper around :func:`aks`."""
    return bool(aks(n, timeout=timeout).verdict)


__all__ = ["AKSStats", "AKSTimeout", "aks", "is_prime_aks", "is_perfect_power"]
