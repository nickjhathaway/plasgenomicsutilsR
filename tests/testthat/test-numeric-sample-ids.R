# A sample id made only of digits is still a label. Read from a file it comes back numeric,
# and R's `named_vector[ids]` then indexes by POSITION rather than by name -- so a lookup
# meant to find a sample's group silently returns whoever sits in that row. Ids are coerced
# to character wherever one enters the package; these hold that line.

test_that("metadata sample ids are made character on the way in", {
  meta <- data.frame(sample = c(4064928010, 4074600944), region = c("East", "North"))
  out <- .normalise_meta(meta)
  expect_type(out$sample, "character")
  expect_identical(out$sample, c("4064928010", "4074600944"))
})

test_that("a large id does not come back in scientific notation", {
  expect_identical(.as_id_chr(4064928010), "4064928010")
  expect_identical(.as_id_chr(c(1e15, 2)), c("1000000000000000", "2"))
})

test_that("NA stays NA rather than becoming the string 'NA'", {
  expect_identical(.as_id_chr(c(1, NA)), c("1", NA_character_))
})

test_that("a case-variant sample column is coerced too", {
  meta <- data.frame(Sample = c(1, 2), region = c("East", "North"))
  expect_message(out <- .normalise_meta(meta), "reading metadata column")
  expect_type(out$sample, "character")
})

test_that("character ids are left exactly as they are", {
  meta <- data.frame(sample = c("RCN13010", "RCN13048"), region = c("East", "North"))
  expect_identical(.normalise_meta(meta)$sample, c("RCN13010", "RCN13048"))
})

# ---- the failure this file exists for -------------------------------------------------

.numeric_id_fixture <- function(ids, dir) {
  blocks <- do.call(rbind, lapply(1:3, function(i) do.call(rbind, lapply((i + 1):4, function(j)
    data.frame(sample1 = ids[i], sample2 = ids[j], chr = "Pf3D7_07_v3",
               start = 403000, end = 407000, different = 0, Nsnp = 50)))))
  meta <- data.frame(sample = ids, region = c("East", "East", "North", "North"))
  ibd_results(blocks = blocks, meta = meta, group_col_in_meta = "region",
              genes = PF_EXAMPLE_DRUG_GENES, min_block_snp = 0, min_block_kb = 0)
}

test_that("a numeric-id cohort groups exactly like the same cohort named with strings", {
  num <- .numeric_id_fixture(c(4064928010, 4074600944, 4074608135, 4089106120))
  chr <- .numeric_id_fixture(c("a", "b", "c", "d"))
  expect_type(num$get_analyzed_samples(), "character")

  as_num <- gene_ibd_overlap(num, genes = "pfcrt")
  as_chr <- gene_ibd_overlap(chr, genes = "pfcrt")
  expect_equal(as_num$n_pairs_total, as_chr$n_pairs_total)
  expect_equal(as_num$group_a, as_chr$group_a)
  expect_equal(as_num$group_b, as_chr$group_b)
})

test_that("small numeric ids are looked up by name, not by row position", {
  # ids whose order differs from their row order: positional indexing gets this wrong
  # without ever saying so, which is what makes it worth a test of its own
  ibd <- .numeric_id_fixture(c(3, 1, 4, 2))
  meta <- ibd$get_meta()
  s2g <- stats::setNames(as.character(meta$region), as.character(meta$sample))
  by_name <- s2g[as.character(c(3, 1, 4, 2))]
  expect_equal(unname(s2g[ibd$get_analyzed_samples()]), unname(by_name))
})

test_that("a bare numeric vector is still refused, and the error says which thing to do", {
  # deliberately NOT coerced: `subset_genotypes(G, 1:3)` almost always means "the first
  # three", and silently reading it as ids would answer a question nobody asked. Ids read
  # from a file are a different matter -- those are coerced where they are read.
  expect_error(.as_sample_ids(c(4064928010, 2)), "selects by id, not by position")
  expect_error(.as_sample_ids(c(4064928010, 2)), "as.character")
  expect_identical(.as_sample_ids(c("x", "y")), c("x", "y"))
})
