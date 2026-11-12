"""Strong Miller scans: conditional, exploratory and exhaustive diagnostic modes.

B denotes the largest tested base, inclusive: 2..min(B,n-2). Tiny
exploratory bounds can produce an empty scan; never silently add base 2.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from fractions import Fraction
import time
import gmpy2

@dataclass(frozen=True)
class Cutoff:
    mode: str
    C: float | None = None
    k: int | None = None
    log_base: str = 'natural'
    source: str = ''
    proof_status: str = ''

    def describe(self):
        if self.mode == 'miller_grh_reference': return 'ceil(2*ln(n)^2); inclusive bases 2..min(B,n-2)'
        if self.mode == 'miller_qrh_official': return 'not runnable as theorem-backed code; no explicit cutoff'
        if self.mode == 'miller_exact_cutoff_diagnostic': return 'full scan 2..n-2; small-input diagnostic'
        arg='ln(n)' if self.log_base=='natural' else 'bit_length(n)'
        return f'ceil({self.C}*{arg}^{self.k}); inclusive bases 2..min(B,n-2)'


def cutoff_for(n: int, spec: Cutoff) -> int | None:
    if spec.mode=='miller_qrh_official': return None
    if n<=3: return max(0,n-2)
    if spec.mode=='miller_exact_cutoff_diagnostic': return n-2
    C,k = (2,2) if spec.mode=='miller_grh_reference' else (spec.C,spec.k)
    if C is None or k is None or C<=0 or k<1: raise ValueError('positive C,k required')
    if spec.log_base=='bits' and spec.mode!='miller_grh_reference':
        b=Fraction(str(C))*n.bit_length()**k
        return min(n-2,-(-b.numerator//b.denominator))
    # An interval computation prevents floating-point underestimation of a
    # theorem cutoff. Increase precision until both endpoint ceilings agree.
    for precision in (128,256,512,1024):
        ends=[]
        for rounding in (gmpy2.RoundDown,gmpy2.RoundUp):
            with gmpy2.context(precision=precision,round=rounding):
                v=gmpy2.mpfr(str(C))*gmpy2.log(gmpy2.mpfr(n))**k
                ends.append(min(n-2,int(gmpy2.ceil(v))))
        if ends[0]==ends[1]: return ends[0]
    raise ArithmeticError('unable to resolve cutoff ceiling')

@dataclass
class MillerResult:
    prime: bool | None = False
    status: str = ''
    cutoff: int | None = None
    bases_tested: int = 0
    first_witness: int | None = None
    modular_exponentiations: int = 0
    explicit_squarings: int = 0
    # GMP's internal powmod multiplication count is unavailable.
    modular_multiplications: int | None = None
    gcd_checks: int = 0
    gcd_shortcut_used: bool = False
    setup_seconds: float = 0.0
    cutoff_seconds: float = 0.0
    arithmetic_seconds: float = 0.0
    total_seconds: float = 0.0
    def asdict(self): return asdict(self)


def _strong_pass(n, a, d, s, result):
    result.modular_exponentiations+=1
    x=gmpy2.powmod(a,d,n)
    if x==1 or x==n-1: return True
    for _ in range(1,s):
        x=x*x%n
        result.explicit_squarings+=1
        if x==n-1: return True
        if x==1: return False
    return False


def strong_miller(n: int, spec: Cutoff, *, gcd_shortcut=True, timeout=None):
    t0=time.perf_counter(); r=MillerResult(); ar0=None
    try:
        if spec.mode=='miller_qrh_official':
            r.prime=None; r.status='not_runnable'; return r
        if n<2: r.status='below_2'; return r
        if n in (2,3): r.prime=True; r.status='small_prime'; return r
        if n%2==0: r.status='even_composite'; return r
        n=gmpy2.mpz(n)
        if gmpy2.is_power(n): r.status='perfect_power_composite'; return r
        s=int(gmpy2.bit_scan1(n-1)); d=(n-1)>>s
        r.setup_seconds=time.perf_counter()-t0
        ct=time.perf_counter(); r.cutoff=cutoff_for(int(n),spec)
        r.cutoff_seconds=time.perf_counter()-ct
        ar0=time.perf_counter()
        for a in range(2,r.cutoff+1):
            if timeout is not None and time.perf_counter()-t0>=timeout:
                r.prime=None; r.status='timeout'; return r
            r.bases_tested+=1
            if gcd_shortcut:
                r.gcd_checks+=1
                if gmpy2.gcd(a,n)>1:
                    r.first_witness=a; r.gcd_shortcut_used=True; r.status='witness_gcd'; return r
            if not _strong_pass(n,a,d,s,r):
                r.first_witness=a; r.status='witness_strong'; return r
        r.prime=True; r.status='all_bases_pass'; return r
    finally:
        r.total_seconds=time.perf_counter()-t0
        if ar0 is not None: r.arithmetic_seconds=time.perf_counter()-ar0
        elif not r.setup_seconds: r.setup_seconds=r.total_seconds


def exact_least_witness(n):
    """Exhaustive small-input oracle; no probable-prime shortcut."""
    if n<4: return None
    if n%2==0: return 2
    m=gmpy2.mpz(n); s=int(gmpy2.bit_scan1(m-1)); d=(m-1)>>s
    for a in range(2,n-1):
        if not _strong_pass(m,a,d,s,MillerResult()): return a
    return None


def spec_dict(spec):
    return {**asdict(spec),'formula':spec.describe()}
