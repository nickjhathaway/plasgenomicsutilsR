test_that("the falciparum URLs are the ones the annotation actually lives at", {
  expect_equal(ensembl_gff_url("falciparum"),
    paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/gff3/",
           "plasmodium_falciparum/Plasmodium_falciparum.GCA000002765v3.63.gff3.gz"))
  expect_equal(ensembl_genome_url("falciparum"),
    paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/fasta/",
           "plasmodium_falciparum/dna/Plasmodium_falciparum.GCA000002765v3.dna.toplevel.fa.gz"))
})

test_that("a species can be named several ways", {
  want <- ensembl_gff_url("plasmodium_falciparum")
  for (nm in c("pf", "PF", "falciparum", "Falciparum", "P. falciparum", "P falciparum",
               "Plasmodium falciparum", " plasmodium_falciparum "))
    expect_equal(ensembl_gff_url(nm), want, info = nm)
})

test_that("the other species are reachable, by short name and in full", {
  expect_match(ensembl_gff_url("vivax"), "plasmodium_vivax_gca900093555/", fixed = TRUE)
  expect_match(ensembl_gff_url("pv"), "Plasmodium_vivax_gca900093555.GCA900093555v2.63.gff3.gz",
               fixed = TRUE)
  expect_match(ensembl_gff_url("knowlesi"), "GCA000006355v3", fixed = TRUE)
  expect_match(ensembl_gff_url("berghei"), "GCA900002375v2", fixed = TRUE)
  expect_match(ensembl_gff_url("chabaudi"), "Plasmodium_chabaudi.PCHAS01.63.gff3.gz", fixed = TRUE)
  expect_match(ensembl_gff_url("yoelii"), "GCA900002385v2", fixed = TRUE)
})

test_that("a genome filed under a collection gets the collection directory", {
  expect_equal(ensembl_gff_url("plasmodium_malariae_gca_900090045"),
    paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/gff3/",
           "protists_alveolata1_collection/plasmodium_malariae_gca_900090045/",
           "Plasmodium_malariae_gca_900090045.PmUG01.63.gff3.gz"))
  expect_match(ensembl_genome_url("plasmodium_malariae_gca_900090045"),
               "fasta/protists_alveolata1_collection/plasmodium_malariae_gca_900090045/dna/",
               fixed = TRUE)
  # a top-level genome must not gain one
  expect_false(grepl("collection", ensembl_gff_url("falciparum")))
})

test_that("masking picks the fasta flavour", {
  expect_match(ensembl_genome_url("pf", masking = "none"), ".dna.toplevel", fixed = TRUE)
  expect_match(ensembl_genome_url("pf", masking = "soft"), ".dna_sm.toplevel", fixed = TRUE)
  expect_match(ensembl_genome_url("pf", masking = "hard"), ".dna_rm.toplevel", fixed = TRUE)
  expect_error(ensembl_genome_url("pf", masking = "lower"), "arg")
})

test_that("the release sets the directory and the release number in the gff name", {
  u <- ensembl_gff_url("pf", release = 62)
  expect_match(u, "/release-62/", fixed = TRUE)
  expect_match(u, ".62.gff3.gz", fixed = TRUE)
  # the genome fasta is not stamped with the release, only its directory is
  g <- ensembl_genome_url("pf", release = 62)
  expect_match(g, "/release-62/", fixed = TRUE)
  expect_false(grepl("[.]62[.]", g))
})

test_that("`assembly` builds a URL without any lookup", {
  # release 60 carried a different assembly token, and is not the table we ship
  expect_equal(ensembl_gff_url("plasmodium_falciparum", release = 60, assembly = "ASM276v2"),
    paste0("https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-60/gff3/",
           "plasmodium_falciparum/Plasmodium_falciparum.ASM276v2.60.gff3.gz"))
  expect_match(ensembl_gff_url("some_other_genome", assembly = "ASM1v1",
                               collection = "protists_alveolata1_collection"),
               "gff3/protists_alveolata1_collection/some_other_genome/Some_other_genome.ASM1v1.63.gff3.gz",
               fixed = TRUE)
})

test_that("an unknown species errors, and names the near misses", {
  # P. malariae has two assemblies and no Ensembl reference, so it has to be named in full
  expect_error(ensembl_gff_url("malariae"), "plasmodium_malariae_gca_900088575")
  expect_error(ensembl_gff_url("malariae"), "plasmodium_malariae_gca_900090045")
  expect_error(ensembl_gff_url("nonesuch"), "see ensembl_species\\(\\)")
  # a bare epithet picks the genome Ensembl gives a top-level directory, over its
  # collection-filed relatives
  expect_equal(ensembl_species(pattern = "^plasmodium_vivax")$collection,
               c("protists_alveolata1_collection", "", rep("protists_alveolata1_collection", 3)))
  expect_match(ensembl_gff_url("vivax"), "plasmodium_vivax_gca900093555", fixed = TRUE)
  expect_error(ensembl_gff_url(c("pf", "pv")), "one species name")
  expect_error(ensembl_gff_url("pf", release = "sixty"), "release")
  expect_error(ensembl_gff_url("pf", release = -1), "release")
})

test_that("ensembl_species lists what the URL builders accept", {
  sp <- ensembl_species()
  expect_s3_class(sp, "tbl_df")
  expect_equal(names(sp), c("species", "assembly", "collection"))
  expect_true(all(grepl("^plasmodium", sp$species)))
  expect_true("plasmodium_falciparum" %in% sp$species)
  # every listed genome round-trips through the URL builder
  expect_true(all(vapply(sp$species, function(s) grepl(s, ensembl_gff_url(s), fixed = TRUE),
                         logical(1))))
  expect_setequal(ensembl_species(pattern = "malariae")$species,
                  c("plasmodium_malariae_gca_900088575", "plasmodium_malariae_gca_900090045"))
})
