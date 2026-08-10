# Run sNMF (LEA) admixture on a genotype matrix

Writes a `.geno` file (missing coded as 9) and runs
[`LEA::snmf()`](https://rdrr.io/pkg/LEA/man/main_sNMF.html) over a range
of `K`, returning the sNMF project. Cross-entropy is computed so
[`snmf_best_k()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/snmf_best_k.md)
can pick `K`.

## Usage

``` r
run_snmf(
  geno,
  K = 1:10,
  rep = 10,
  alpha = 10,
  seed = 42,
  cpu = 1,
  cache = TRUE,
  cache_dir = NULL,
  verbose = FALSE,
  log_file = NULL
)
```

## Arguments

- geno:

  A genotype matrix (samples x SNPs, 0/1/2, `NA` allowed) or the list
  from
  [`load_genotypes()`](https://nickjhathaway.github.io/plasgenomicsutilsR/reference/load_genotypes.md).

- K:

  Integer vector of ancestral-population counts to fit (default `1:10`).

- rep:

  Repetitions per `K`.

- alpha:

  Regularisation.

- seed:

  Random seed.

- cpu:

  CPU cores.

- cache:

  Reuse a previously computed project when the genotypes and all
  parameters match (sNMF is slow); set `FALSE` to always recompute.

- cache_dir:

  Directory for the cached `.geno`/project (default a per-session temp
  dir). Point it at a persistent path to reuse across sessions.

- verbose:

  Show LEA's (voluminous) console output. Default `FALSE` runs it
  quietly; set `log_file` to capture it instead of discarding it.

- log_file:

  Optional file to write LEA's output to when `verbose = FALSE`.

## Value

An `snmf_fit` object: a list with the LEA `project`, the fitted `K`
range, and `samples`.
