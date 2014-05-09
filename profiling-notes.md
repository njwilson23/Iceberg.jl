# Profiling notes

## Results from profiling 2d run

Ran `test_hill_onephase_2d` to `tend=0.1` with `n=32`

- 9600 samples
- about one-third spent solving the heat equation (4100 samples)
    - only about a small amount of this time doing assembly (500 samples)
    - mostly doing linalg (3600 samples)
        - are there better sparse solvers - problem is tridiagonal
        - majority of time doing the LU factorization in umfpack
- ten-percent of time computing the front velocity (1032 samples)
- nearly half of time reinitalizing (3900 samples)
    - again, three-quarters of which is computing the LU factorization




