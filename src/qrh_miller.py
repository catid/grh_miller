"""QRH-style Miller scans with explicit proof-status boundaries.

The OpenAI quasi-RH material reviewed by this project does not publish a
numeric constant for an executable strong-Miller witness bound.  This module
therefore exposes two different operations:

``assumed_qrh_bound``
    Evaluate a caller-supplied ``C * log(n)**k`` model.  The caller is making
    an assumption; the result is not an unconditional theorem-backed cutoff.

``qrh_exact_with_aks_fallback``
    Optionally scan that prefix, then run exact textbook AKS whenever the
    prefix did not produce a compositeness witness.  The complete decision is
    unconditional because AKS supplies the final decision.  The QRH model
    only offers a possible early exit.

``exhaustive_miller``
    Scan every admissible base.  This is an exact diagnostic for small inputs,
    but has no polynomial-time claim.

No guessed value of ``C`` is supplied here.  In particular, ``C=1`` is an
exploratory model and must not be described as an unconditional consequence
of the QRH papers.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import time
from typing import Any

from .aks import AKSStats, AKSTimeout, aks
from .miller import Cutoff, MillerResult, cutoff_for, strong_miller


def _validate_model(C: Any, exponent: int, log_base: str) -> None:
    if C is None:
        raise ValueError("an explicit positive C is required")
    try:
        c = float(C)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("C must be a finite positive number") from exc
    if not math.isfinite(c) or c <= 0:
        raise ValueError("C must be a finite positive number")
    if isinstance(exponent, bool) or not isinstance(exponent, int) or exponent < 1:
        raise ValueError("exponent must be a positive integer")
    if log_base not in {"natural", "bits"}:
        raise ValueError("log_base must be 'natural' or 'bits'")


def assumed_qrh_bound(
    n: int,
    C: Any,
    exponent: int = 32,
    *,
    log_base: str = "natural",
) -> int:
    """Return ``ceil(C * log(n)**exponent)``, capped at ``n - 2``.

    This is a parameterized assumption model.  It is deliberately not named
    or documented as an unconditional QRH theorem.  ``C`` is mandatory so a
    missing constant cannot silently become ``1``.
    """

    _validate_model(C, exponent, log_base)
    n = int(n)
    if n <= 2:
        return 0
    spec = Cutoff(
        mode="miller_qrh_parametric",
        C=C,
        k=exponent,
        log_base=log_base,
        source="caller-supplied QRH assumption model",
        proof_status="conditional/exploratory; not an explicit theorem cutoff",
    )
    value = cutoff_for(n, spec)
    assert value is not None
    return int(value)


def assumed_qrh_metadata(
    C: Any,
    exponent: int = 32,
    *,
    log_base: str = "natural",
) -> dict[str, Any]:
    """Describe the assumption model in a machine-readable form."""

    _validate_model(C, exponent, log_base)
    argument = "ln(n)" if log_base == "natural" else "bit_length(n)"
    return {
        "formula": f"ceil({C} * {argument}^{exponent})",
        "C": str(C),
        "exponent": exponent,
        "log_base": log_base,
        "proof_status": "conditional/exploratory",
        "source": "caller-supplied model; not official (no numerical QRH cutoff)",
    }


@dataclass
class QRHExactResult:
    """Result and provenance for the exact AKS-fallback wrapper."""

    prime: bool
    status: str
    proof_status: str
    cutoff: int | None = None
    prefix_used: bool = False
    prefix_result: dict[str, Any] | None = None
    aks_result: dict[str, Any] | None = None
    total_seconds: float = 0.0

    @property
    def verdict(self) -> bool:
        """Alias matching :class:`src.aks.AKSStats`."""

        return self.prime

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def asdict(self) -> dict[str, Any]:
        """Compatibility alias for :class:`MillerResult`."""

        return self.as_dict()


def qrh_exact_with_aks_fallback(
    n: int,
    C: Any | None = None,
    *,
    exponent: int = 32,
    log_base: str = "natural",
    timeout: float | None = None,
    gcd_shortcut: bool = True,
) -> QRHExactResult:
    """Return an unconditional exact decision using an optional QRH prefix.

    When ``C`` is omitted this is simply exact textbook AKS.  When supplied,
    a caller-selected prefix is scanned first.  A witness proves compositeness
    immediately; every other case is delegated to AKS.  Thus the wrapper is
    exact without assuming that the QRH theorem supplies a usable constant.
    ``timeout`` is a total wall-clock budget; an :class:`AKSTimeout` is raised
    if the final AKS call exhausts it.
    """

    n = int(n)
    started = time.perf_counter()
    prefix: MillerResult | None = None
    cutoff: int | None = None

    if C is not None:
        _validate_model(C, exponent, log_base)
        cutoff = assumed_qrh_bound(n, C, exponent, log_base=log_base)
        spec = Cutoff(
            mode="miller_qrh_parametric",
            C=C,
            k=exponent,
            log_base=log_base,
            source="caller-supplied QRH assumption model",
            proof_status="conditional/exploratory; prefix only",
        )
        prefix_timeout = None
        if timeout is not None:
            prefix_timeout = max(0.0, float(timeout))
        prefix = strong_miller(
            n,
            spec,
            gcd_shortcut=gcd_shortcut,
            timeout=prefix_timeout,
        )
        if prefix.prime is False:
            # A strong Miller witness is a certificate by itself.  Boundary,
            # parity, and perfect-power exits are also exact, but they are
            # prechecks rather than witness bases; preserve that provenance
            # in the wrapper's status fields.
            has_witness = prefix.first_witness is not None
            return QRHExactResult(
                prime=False,
                status="prefix_witness" if has_witness else "prefix_precheck",
                proof_status=(
                    "unconditional_witness" if has_witness
                    else "unconditional_precheck"
                ),
                cutoff=cutoff,
                prefix_used=True,
                prefix_result=prefix.asdict(),
                total_seconds=time.perf_counter() - started,
            )

    elapsed = time.perf_counter() - started
    remaining = None if timeout is None else max(0.0, float(timeout) - elapsed)
    # n < 2 and small boundary values are handled exactly by AKS as well.
    try:
        aks_stats: AKSStats = aks(n, timeout=remaining)
    except AKSTimeout:
        raise
    return QRHExactResult(
        prime=bool(aks_stats.verdict),
        status="aks_fallback" if prefix is not None else "aks_only",
        proof_status="unconditional_aks",
        cutoff=cutoff,
        prefix_used=prefix is not None,
        prefix_result=None if prefix is None else prefix.asdict(),
        aks_result=aks_stats.as_dict(),
        total_seconds=time.perf_counter() - started,
    )


def exhaustive_miller(
    n: int,
    *,
    timeout: float | None = None,
    gcd_shortcut: bool = True,
) -> MillerResult:
    """Run the exact all-base diagnostic scan through ``n - 2``.

    The returned object is the regular :class:`MillerResult`; its ``status``
    and ``cutoff`` identify the exhaustive mode.  This is intended for small
    inputs and witness-rank diagnostics, not production-sized primality tests.
    """

    spec = Cutoff(
        mode="miller_exact_cutoff_diagnostic",
        source="exhaustive base scan",
        proof_status="unconditional exact diagnostic",
    )
    return strong_miller(n, spec, timeout=timeout, gcd_shortcut=gcd_shortcut)


__all__ = [
    "QRHExactResult",
    "assumed_qrh_bound",
    "assumed_qrh_metadata",
    "exhaustive_miller",
    "qrh_exact_with_aks_fallback",
]
