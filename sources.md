# Source notes

The implementation and README use these primary references:

- Gary L. Miller, “Riemann’s hypothesis and tests for primality”: <https://www.cs.cmu.edu/~glmiller/Publications/Papers/Mi75.pdf>
- Eric Bach, “Explicit Bounds for Primality Testing and Related Problems”: <https://doi.org/10.1090/S0025-5718-1990-1023756-8>
- Agrawal, Kayal, Saxena, “PRIMES is in P”: <https://annals.math.princeton.edu/2004/160-2/p12>
- Lenstra and Pomerance, “Primality testing with Gaussian periods”: <https://math.dartmouth.edu/~carlp/aks111216.pdf>
- OpenAI, “The Quasi-Riemann Hypothesis,” September 30, 2026: <https://github.com/openai/math/tree/main/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026>
- OpenAI, “The Quasi-Riemann Hypothesis,” October 5, 2026: <https://github.com/openai/math/blob/main/preprints/The-Quasi-Riemann-Hypothesis-October-5-2026/build/paper2.tex>

Bach’s ERH conditional consequence gives the reference cutoff `B(n) = ceil(2 ln(n)^2)` using the natural logarithm. The OpenAI QRH material gives no explicit constant for a runnable arbitrary-composite strong-witness cutoff. Its September corollary concerns a least quadratic nonresidue modulo a prime and leaves `C` unspecified; its October introduction states a deterministic polynomial-time Miller application without a numeric cutoff. Those claims remain separate from the measured exploratory `C,k` sweeps.
