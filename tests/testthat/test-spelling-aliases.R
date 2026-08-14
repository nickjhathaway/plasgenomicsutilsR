# Colour arguments accept either spelling. The package's own naming was mixed -- some
# functions took `colours`, others `colors`, two took both -- so the guard below is the
# point of these tests: a new colour argument that ships in one spelling only fails here.

.SPELLING_PAIRS <- list(c("colours", "colors"), c("colour", "color"),
                        c("color_group", "colour_group"), c("na_colour", "na_color"),
                        c("edge_colour", "edge_color"), c("border_colour", "border_color"),
                        c("group_colours", "group_colors"),
                        c("annotate_colours", "annotate_colors"))

test_that("every exported colour argument exists in both spellings and is resolved", {
  ns <- asNamespace("plasgenomicsutilsR")
  lonely <- character()
  unresolved <- character()
  for (nm in sort(getNamespaceExports(ns))) {
    f <- get(nm, envir = ns)
    if (!is.function(f)) next
    args <- names(formals(f))
    body_txt <- paste(deparse(body(f)), collapse = " ")
    for (p in .SPELLING_PAIRS) {
      if (any(p %in% args) && !all(p %in% args))
        lonely <- c(lonely, sprintf("%s(%s)", nm, paste(intersect(p, args), collapse = ", ")))
      if (all(p %in% args) &&
          !any(vapply(p, function(a) grepl(sprintf('alias_arg("%s"', a), body_txt,
                                           fixed = TRUE), logical(1))))
        unresolved <- c(unresolved, sprintf("%s(%s)", nm, p[1]))
    }
  }
  expect_equal(lonely, character())          # one spelling only
  expect_equal(unresolved, character())      # both formals, but neither read
})

test_that("the two spellings produce the same plot", {
  skip_if_not_installed("ggplot2")
  pal <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
  ps <- example_pop_structure(umap = FALSE)
  ps$run_snmf(K = 1:3, rep = 2, cache = FALSE)
  q <- ps$q(3)
  built <- function(p) ggplot2::ggplot_build(p)$data

  expect_equal(built(plot_admixture(q, colours = pal)),
               built(plot_admixture(q, colors = pal)))
  expect_equal(built(plot_admixture(q, border_colour = "grey30")),
               built(plot_admixture(q, border_color = "grey30")))
  expect_equal(built(ps$plot_pca(colour = "country")),
               built(ps$plot_pca(color = "country")))

  # naming both at once is a mistake, not a silent win for one of them
  expect_error(plot_admixture(q, colours = pal, colors = pal),
               "either `colours` or `colors`")
  expect_error(ps$plot_pca(colour = "country", color = "country"),
               "either `colour` or `color`")
})

test_that("an argument whose default is a real colour still defaults", {
  skip_if_not_installed("ggplot2")
  ps <- example_pop_structure(umap = FALSE)
  ps$run_snmf(K = 1:3, rep = 2, cache = FALSE)
  q <- ps$q(3)
  bar_colour <- function(p) {
    d <- ggplot2::ggplot_build(p)$data
    unique(d[[which(vapply(d, function(x) "ymin" %in% names(x), logical(1)))[1]]]$colour)
  }
  # `border_colour = "black"` is the default; passing neither spelling must not read the
  # alias's NULL, and either spelling must override it
  expect_equal(bar_colour(plot_admixture(q)), "black")
  expect_equal(bar_colour(plot_admixture(q, border_color = "grey30")), "grey30")
  expect_equal(bar_colour(plot_admixture(q, border_colour = "grey30")), "grey30")
})

test_that("PopStructure takes either spelling for its colour map", {
  ps <- example_pop_structure(umap = FALSE)
  cols <- list(country = c(Ghana = "#111111", Cambodia = "#222222"))
  ps$set_colors(colours = cols)
  expect_equal(ps$get_colors()$country, cols$country)
  expect_identical(ps$get_colours(), ps$get_colors())
})
