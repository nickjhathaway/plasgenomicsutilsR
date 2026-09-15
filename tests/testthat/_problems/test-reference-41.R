# Extracted from test-reference.R:41

# test -------------------------------------------------------------------------
skip_if_not_installed("yaml")
yml <- file.path("..", "..", "_pkgdown.yml")
skip_if_not(file.exists(yml))
cfg <- yaml::read_yaml(yml)
listed <- unlist(lapply(cfg$articles, `[[`, "contents"), use.names = FALSE)
present <- tools::file_path_sans_ext(
    list.files(file.path("..", "..", "vignettes"), pattern = "[.]Rmd$"))
expect_setequal(present, listed)
ref <- unlist(lapply(cfg$reference, `[[`, "contents"), use.names = FALSE)
ref <- sub("^starts_with\\(\"(.*)\"\\)$", "\\1", trimws(ref))
exports <- getNamespaceExports("plasgenomicsutilsR")
missing <- exports[!vapply(exports, function(f)
    any(f == ref) || any(startsWith(f, ref[!ref %in% exports])), logical(1))]
expect_equal(missing, character(0))
