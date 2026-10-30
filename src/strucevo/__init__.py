"""strucevo: pairwise sequence / 3Di / structure divergence measures and curve models.

Modules
-------
distances : p / Poisson / Kimura distances, batched GPU ML distances under any
            reversible 20-state model (+G, +I), GTR estimation (for 3Di), TM/LDDT
            transforms, batched Kabsch RMSD.
models    : y-vs-x model fits (linear through origin, linear, saturating, power)
            with AIC/BIC and out-of-range prediction error.
"""
__version__ = "0.1.0"
GLOBAL_SEED = 20261030
