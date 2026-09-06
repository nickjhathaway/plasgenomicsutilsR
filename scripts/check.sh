#!/usr/bin/env bash
# Run what CI runs, before pushing.
#
# CI for this repo is two jobs (.github/workflows/R-CMD-check.yaml): `R CMD check` across
# an R matrix, and `pkgdown::check_pkgdown()`, which fails on a topic in _pkgdown.yml that
# no longer exists or an export missing from the reference index. The second takes seconds
# -- `--fast` runs it plus the test suite and skips the full check, for the loop between
# edits.
#
#   scripts/check.sh          # everything CI runs
#   scripts/check.sh --fast   # pkgdown config + testthat only (under a minute)
set -uo pipefail
cd "$(dirname "$0")/.."

FAST=0
[[ "${1:-}" == "--fast" ]] && FAST=1
FAILED=()
R=${RSCRIPT:-Rscript}

run() {                       # run <label> <r-expression>
  local label=$1; shift
  printf '\n\033[1m== %s ==\033[0m\n' "$label"
  if "$R" -e "$1"; then
    printf '\033[32mOK\033[0m  %s\n' "$label"
  else
    printf '\033[31mFAILED\033[0m  %s\n' "$label"
    FAILED+=("$label")
  fi
}

# roxygen output committed and current: CI checks a tree where man/ and NAMESPACE are
# whatever was pushed, so a forgotten document() shows up there as a check failure
run "roxygen is up to date" '
  before <- tools::md5sum(list.files(c("man", "."), pattern = "[.](Rd|NAMESPACE)$",
                                     full.names = TRUE, recursive = FALSE))
  suppressMessages(devtools::document(quiet = TRUE))
  after <- tools::md5sum(names(before))
  changed <- names(before)[before != after | is.na(after)]
  if (length(changed)) stop("run devtools::document() and commit: ",
                            paste(basename(changed), collapse = ", "))
  cat("man/ and NAMESPACE match the roxygen comments\n")'

run "pkgdown config" 'pkgdown::check_pkgdown()'

run "testthat" 'suppressMessages(devtools::load_all(quiet = TRUE));
  r <- as.data.frame(testthat::test_local(reporter = "silent"));
  bad <- sum(r$failed) + sum(r$error);
  if (bad) stop(bad, " test failure(s)") else cat(sum(r$passed), "tests passed\n")'

if [[ $FAST -eq 0 ]]; then
  run "R CMD check" 'devtools::check(document = FALSE,
    args = c("--no-manual"), error_on = "warning")'
fi

printf '\n'
if [[ ${#FAILED[@]} -eq 0 ]]; then
  printf '\033[32mall checks passed\033[0m\n'
else
  printf '\033[31m%d check(s) failed:\033[0m %s\n' "${#FAILED[@]}" "${FAILED[*]}"
  exit 1
fi
