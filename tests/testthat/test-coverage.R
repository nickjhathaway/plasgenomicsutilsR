# Coverage QC: the verdict table and the three plots. The two floors are independent --
# selective whole-genome amplification can give a respectable mean while leaving much of the
# genome at zero, and only the breadth column shows that -- so the tests check they can fail
# separately and that `fail_reason` says which.

.cov_table <- function(n = 10) {
  set.seed(7)
  ch <- sprintf("Pf3D7_%02d_v3", 1:5)
  s <- paste0("s", seq_len(n))
  per_chr <- expand.grid(sample = s, chrom = ch, stringsAsFactors = FALSE)
  per_chr$mean <- stats::runif(nrow(per_chr), 5, 80)
  # a genome-wide row per sample: the one coverage_qc() reduces to
  genome <- data.frame(sample = s, chrom = "genome",
                       mean = c(60, 40, 30, 3, 55, 45, 50, 35, 2, 70),
                       stringsAsFactors = FALSE)
  df <- rbind(per_chr, genome)
  df$median <- df$mean
  df$sd <- df$mean / 4
  # s3 is deep but narrow (amplification), s4/s9 are simply shallow
  df$pct_ge_10x <- pmin(100, df$mean * 2.4)
  df$pct_ge_10x[df$sample == "s3" & df$chrom == "genome"] <- 41
  df$pct_zero <- pmax(0, 100 - df$pct_ge_10x)
  tibble::as_tibble(df)
}

test_that("coverage_qc reduces to one row per sample and applies both floors", {
  qc <- coverage_qc(.cov_table())
  expect_equal(nrow(qc), 10)
  expect_setequal(qc$sample, paste0("s", 1:10))
  expect_true(all(c("mean", "pct_ge_10x", "pass", "fail_reason") %in% names(qc)))

  # shallow samples fail on the mean; the deep-but-narrow one fails on breadth alone
  by_s <- as.data.frame(qc)[match(paste0("s", 1:10), qc$sample), ]
  expect_false(by_s$pass[by_s$sample == "s4"])     # mean 3
  expect_false(by_s$pass[by_s$sample == "s9"])     # mean 2
  expect_false(by_s$pass[by_s$sample == "s3"])     # mean 30 but only 41% at 10x
  expect_true(by_s$pass[by_s$sample == "s1"])

  # the reason distinguishes the two failure modes, and names both when both apply
  expect_equal(by_s$fail_reason[by_s$sample == "s3"], "low breadth")
  expect_equal(by_s$fail_reason[by_s$sample == "s4"], "low depth and breadth")
  expect_true(is.na(by_s$fail_reason[by_s$sample == "s1"]))
  # sorted worst-first, so a surprising drop is the first thing read
  expect_false(is.unsorted(qc$mean))
})

test_that("the floors move independently", {
  cov <- .cov_table()
  expect_true(all(coverage_qc(cov, min_mean = 0, min_breadth = 0)$pass))
  # raising only the breadth floor fails more samples, and only on breadth
  strict <- coverage_qc(cov, min_mean = 0, min_breadth = 99)
  expect_true(sum(!strict$pass) > 0)
  expect_setequal(strict$fail_reason[!strict$pass], "low breadth")   # never "low depth"
})

test_that("read_coverage reads a written table back with its text columns intact", {
  skip_if_not_installed("readr")
  f <- tempfile(fileext = ".tsv")
  cov <- .cov_table()
  utils::write.table(cov, f, sep = "\t", quote = FALSE, row.names = FALSE)
  back <- read_coverage(f)
  expect_equal(nrow(back), nrow(cov))
  # `sample` and `chrom` must stay character: a cohort of numeric ids would otherwise
  # come back as numbers and stop matching the metadata
  expect_type(back$sample, "character")
  expect_type(back$chrom, "character")
  expect_equal(sort(unique(back$chrom)), sort(unique(cov$chrom)))
})

test_that("the three coverage plots build", {
  skip_if_not_installed("ggplot2")
  cov <- .cov_table()
  expect_s3_class(plot_coverage_summary(cov), "ggplot")
  expect_s3_class(plot_coverage_summary(cov, label_failures = FALSE), "ggplot")
  expect_s3_class(plot_coverage_by_chrom(cov), "ggplot")
  expect_s3_class(plot_coverage_by_chrom(cov, relative = FALSE), "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot_coverage_summary(cov)))
  expect_silent(ggplot2::ggplotGrob(plot_coverage_by_chrom(cov)))

  # a metric the table does not carry is named rather than failing inside ggplot
  expect_error(plot_coverage_by_chrom(cov, metric = "nope"), "no `nope` column")
  expect_error(plot_coverage_by_chrom(cov[cov$chrom == "genome", ]), "only genome-wide")
})

test_that("plot_coverage_dropout takes per-window depths or merged regions", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  win <- expand.grid(sample = paste0("s", 1:6),
                     chrom = c("Pf3D7_01_v3", "Pf3D7_02_v3"),
                     start = seq(0, 4e5, by = 1e5), stringsAsFactors = FALSE)
  win$end <- win$start + 1e5
  win$mean_depth <- stats::runif(nrow(win), 0, 30)
  expect_s3_class(plot_coverage_dropout(win), "ggplot")

  # the already-merged form goes straight through
  agg <- data.frame(chrom = "Pf3D7_01_v3", start = c(0, 1e5), end = c(1e5, 2e5),
                    frac_samples_uncovered = c(0.1, 0.95))
  expect_s3_class(plot_coverage_dropout(agg), "ggplot")
  expect_error(plot_coverage_dropout(win[, c("sample", "chrom")]), "needs sample, chrom")
})
