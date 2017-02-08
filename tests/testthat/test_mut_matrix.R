context("Function 'mut_matrix'")

test_that("transforms correctly", {
    # To test mut_matrix, we need to load the reference genome first.
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
    library(BSgenome)
    library(ref_genome, character.only = TRUE)

    # We re-use the data that is shipped with the package.
    sample_names <- c ( "colon1", "colon2", "colon3",
                        "intestine1", "intestine2", "intestine3",
                        "liver1", "liver2", "liver3" )

    vcfs <- list.files (system.file("extdata", package="MutationalPatterns"),
                        pattern = ".vcf", full.names = TRUE)

    input <- read_vcfs_as_granges(vcfs, sample_names,
                                    "BSgenome.Hsapiens.UCSC.hg19")

    matrix <- mut_matrix (input, ref_genome)

    actual <- colSums(matrix)
    expected <- c(491, 488, 487, 488, 484, 486, 488, 489, 490)
    names(expected) <- sample_names

    expect_that(actual, equals(expected))
})
