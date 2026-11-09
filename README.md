# Degrees of Prolongation Varieties of Monomial Planar Systems

Deep Bhattacharjee

For x1' = x2^p, x2' = x1^q with p, q >= 2, the paper gives the degree of every prolongation variety in closed form. This answers Questions 6 and 7 of G. Pogudin's list of open problems in differential algebra for this family, and proves Puddu's conjecture.

- `paper/`: LaTeX source (amsart, pdflatex, no BibTeX step), figures and PDF
- `code/check_formula.py`: compares the data with the formula and tests the gcd theorem (needs sympy)
- `code/run_degrees.sh p q M`: recomputes deg X_m for m <= M (needs Singular)
- `code/make_figures.py`: draws the figures (needs numpy, sympy, matplotlib)
- `data/fPQ.txt`: degrees computed with Singular, m <= 30

The proofs in the paper do not use the scripts.

Archive: Zenodo, concept DOI [10.5281/zenodo.24094789](https://doi.org/10.5281/zenodo.24094789).

License: MIT.
