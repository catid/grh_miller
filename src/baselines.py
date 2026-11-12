"""Practical probable-prime baselines, explicitly separate from AKS."""
import random,time
import gmpy2
from sympy.external.gmpy import is_strong_bpsw_prp
from .miller import _strong_pass, MillerResult
FIXED_BASES=(2,3,5,7,11,13,17,19,23,29,31,37)

def mr(n,bases):
    r=MillerResult()
    if n<2:return r
    if n in (2,3):r.prime=True;return r
    if n%2==0:return r
    m=gmpy2.mpz(n);s=int(gmpy2.bit_scan1(m-1));d=(m-1)>>s
    for a in bases:
        if a>=n:continue
        r.bases_tested+=1
        if not _strong_pass(m,a,d,s,r):r.first_witness=a;return r
    r.prime=True;return r

def fixed_miller_rabin(n):
    t=time.perf_counter();r=mr(n,FIXED_BASES);r.total_seconds=time.perf_counter()-t
    return r.asdict()

def randomized_miller_rabin(n,rounds=16,seed=20261111):
    t=time.perf_counter(); rng=random.Random(seed ^ n)
    bases=[rng.randrange(2,n-1) for _ in range(rounds)] if n>3 else []
    r=mr(n,bases);r.total_seconds=time.perf_counter()-t
    return r.asdict()

def gmpy2_probable_prime(n):
    t=time.perf_counter(); p=bool(gmpy2.is_prime(n,25))
    return {'prime':p,'total_seconds':time.perf_counter()-t}

def bpsw(n):
    t=time.perf_counter();p=bool(is_strong_bpsw_prp(n)) if n>=2 else False
    return {'prime':p,'total_seconds':time.perf_counter()-t}
