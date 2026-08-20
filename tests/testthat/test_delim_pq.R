library(MiscMetabar)
data(data_fungi_mini)

test_that("is_delim_installed returns a single logical", {
  expect_length(is_delim_installed("asap"), 1)
  expect_type(is_delim_installed("abgd"), "logical")
  expect_false(is_delim_installed("asap", path = tempfile()))
  expect_error(is_delim_installed("gmyc"), "should be one of")
})

test_that("delim_pq aborts when delimtools is absent", {
  skip_if(requireNamespace("delimtools", quietly = TRUE))
  expect_error(delim_pq(data_fungi_mini), "delimtools")
})

test_that("delim_pq aborts without a refseq slot", {
  skip_if_not_installed("delimtools")
  pq <- data_fungi_mini
  pq@refseq <- NULL
  expect_error(delim_pq(pq), "no .*refseq")
})

test_that("delim_pq aborts when the executable is missing", {
  skip_if_not_installed("delimtools")
  skip_if(is_delim_installed("asap"))
  expect_error(delim_pq(data_fungi_mini, method = "asap"), "was not found")
})

test_that("delim_pq validates its method", {
  expect_error(delim_pq(data_fungi_mini, method = "gmyc"), "should be one of")
})

test_that("delim_partitions maps partitions onto taxa names", {
  taxa <- phyloseq::taxa_names(data_fungi_mini)
  res <- data.frame(
    labels = rev(taxa),
    asap = rev(seq_along(taxa)),
    stringsAsFactors = FALSE
  )
  part <- delim_partitions(res, data_fungi_mini, "asap")
  expect_named(part, c("taxa", "partition"))
  expect_equal(part$taxa, taxa)
  expect_equal(part$partition, seq_along(taxa))
})

test_that("delim_partitions returns NA for undelimited taxa", {
  taxa <- phyloseq::taxa_names(data_fungi_mini)
  res <- data.frame(
    labels = taxa[1:5],
    asap = 1:5,
    stringsAsFactors = FALSE
  )
  part <- delim_partitions(res, data_fungi_mini, "asap")
  expect_equal(nrow(part), length(taxa))
  expect_equal(sum(is.na(part$partition)), length(taxa) - 5)
})

test_that("delim_partitions rejects a table without labels", {
  res <- data.frame(x = 1, asap = 1)
  expect_error(
    delim_partitions(res, data_fungi_mini, "asap"),
    "no .*labels"
  )
})

test_that("resolve_delim_exe errors on a missing executable", {
  expect_error(resolve_delim_exe("asap", "/no/such/binary"), "was not found")
})

test_that("align_refseq leaves an aligned set untouched", {
  skip_if_not_installed("Biostrings")
  dna <- Biostrings::DNAStringSet(c(a = "ACGT", b = "ACGA"))
  expect_identical(align_refseq(dna, align = FALSE), dna)
})

test_that("align_refseq refuses unaligned sequences when align is FALSE", {
  skip_if_not_installed("Biostrings")
  dna <- Biostrings::DNAStringSet(c(a = "ACGT", b = "ACG"))
  expect_error(align_refseq(dna, align = FALSE), "different lengths")
})

test_that("align_refseq aligns unequal sequences", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("DECIPHER")
  dna <- Biostrings::DNAStringSet(c(
    a = "ACGTACGT",
    b = "ACGTCGT",
    c = "ACGTACG"
  ))
  aln <- align_refseq(dna, align = TRUE)
  expect_length(unique(Biostrings::width(aln)), 1)
  expect_equal(names(aln), names(dna))
})

test_that("delim_pq delimits and merges taxa with asap", {
  skip_if_not(is_delim_installed("asap"))
  res <- delim_pq(data_fungi_mini, method = "asap")
  expect_s4_class(res, "phyloseq")
  expect_lt(phyloseq::ntaxa(res), phyloseq::ntaxa(data_fungi_mini))
  expect_gt(phyloseq::ntaxa(res), 1)
})

test_that("delim_pq returns a partition table when merge_taxa is FALSE", {
  skip_if_not(is_delim_installed("asap"))
  res <- delim_pq(data_fungi_mini, method = "asap", merge_taxa = FALSE)
  expect_named(res, c("taxa", "partition"))
  expect_identical(res$taxa, phyloseq::taxa_names(data_fungi_mini))
  expect_gt(length(unique(res$partition)), 1)
})

test_that("delim_pq is reproducible across consecutive runs", {
  skip_if_not(is_delim_installed("asap"))
  first <- delim_pq(data_fungi_mini, method = "asap", merge_taxa = FALSE)
  second <- delim_pq(data_fungi_mini, method = "asap", merge_taxa = FALSE)
  expect_identical(first, second)
})

test_that("delim_pq reports a single-partition result instead of merging", {
  skip_if_not(is_delim_installed("abgd"))
  expect_error(
    delim_pq(data_fungi_mini, method = "abgd", slope = 1.5),
    "single partition"
  )
})

test_that("delim_pq delimits with abgd given a workable slope", {
  skip_if_not(is_delim_installed("abgd"))
  res <- delim_pq(data_fungi_mini, method = "abgd", slope = 0.5)
  expect_s4_class(res, "phyloseq")
  expect_lt(phyloseq::ntaxa(res), phyloseq::ntaxa(data_fungi_mini))
})
