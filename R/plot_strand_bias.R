#' Plot strand bias per base substitution type
#'
#' For each base substitution type and transcriptional strand the total
#' number of mutations and the relative contribution within a group is
#' returned.
#'
#' @param strand_bias data.frame, result from strand_bias function
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_alpha_discrete
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 position_dodge
#'
#' @examples
#' ## See the 'mut_matrix_stranded()' example for how we obtained the
#' ## following mutation matrix.
#' mut_mat_s <- readRDS(system.file("states/mut_mat_s_data.R",
#'                                     package="MutationalPatterns"))
#'
#' ## Load a reference genome.
#' ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#' library(ref_genome, character.only = TRUE)
#'
#' tissue <- c("colon", "colon", "colon",
#'             "intestine", "intestine", "intestine",
#'             "liver", "liver", "liver")
#'
#' ## Perform the strand bias test.
#' strand_counts = strand_occurences(mut_mat_s, by=tissue)
#' strand_bias = strand_bias_test(strand_counts)
#'
#' ## Plot the strand bias.
#' plot_strand_bias(strand_bias)
#'
#' @seealso
#' \code{\link{mut_matrix_stranded}},
#' \code{\link{strand_occurences}},
#' \code{\link{strand_bias_test}}
#' \code{\link{plot_strand}}
#'
#' @export

plot_strand_bias = function(strand_bias, colors)
{
    # if colors parameter not provided, set to default colors
    if (missing(colors))
        colors=COLORS6

    # determine max y value for plotting
    # = log2 ratio with pseudo counts
    max = round( max(abs(log2(strand_bias$ratio))), digits = 1) + 0.1

    type = NULL
    ratio = NULL
    significant = NULL

    # plot strand bias with poisson test results
    plot = ggplot(strand_bias, aes(x = type, y = log2(ratio), fill = type)) +
        scale_fill_manual(values = COLORS6) +
        geom_bar(colour = "black", stat ="identity", position = "identity") +
        scale_y_continuous(limits = c(-max, max), breaks = seq(-1, 1, 0.2)) +
        geom_text(
            aes(x = type,
                y = log2(ratio),
                ymax = log2(ratio), 
                label = significant,
                vjust = ifelse(sign(log2(ratio)) > 0, 0.5, 1)),
            size = 8,
            position = ggplot2::position_dodge(width = 1)) +
        facet_grid(. ~ group) +
        theme_bw()  +
        theme(axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                legend.title = element_blank()) +
        xlab("") + 
        ylab("log2(transcribed/untranscribed)") +
        scale_x_discrete(breaks=NULL)

    return(plot)
}
