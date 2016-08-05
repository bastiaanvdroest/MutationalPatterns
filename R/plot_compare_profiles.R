#' Compare two 96 mutation profiles
#' 
#' Plots two 96 mutation profiles and their difference, reports the residual sum of squares (RSS)
#' @param profile1 First 96 mutation profile
#' @param profile2 Second 96 mutation profile
#' @param profile_names Character vector with names of the mutations profiles used for plotting, default = c("profile 1", "profile 2")
#' @param profile_ymax Maximum value of y-axis (relative contribution) for profile plotting, default = 0.15
#' @param diff_ylim Y-axis limits for profile difference plot, default = c(-0.02, 0.02)
#' @param colors 6 value color vector
#' @return 96 spectrum plot of profile 1, profile 2 and their difference
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 element_blank
#' @importFrom BiocGenerics cbind
#' @export

plot_compare_profiles = function(profile1, profile2, profile_names = c("profile 1", "profile 2"), profile_ymax = 0.15, diff_ylim = c(-0.02, 0.02), colors)
{
  # if colors parameter not provided, set to default colors
  if(missing(colors)){colors = COLORS6}
  s1_relative = profile1 / sum(profile1)
  s2_relative = profile2 / sum(profile2)
  diff = s1_relative - s2_relative
  # residual sum of squares
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = T, digits = 3)
  
  x = cbind(s1_relative, s2_relative, diff)
  colnames(x) = c(profile_names, "Difference")
  
  substitutions = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
  index = c(rep(1,1,16), rep(2,1,16), rep(3,1,16), rep(4,1,16), rep(5,1,16), rep(6,1,16))
  # Context
  context = rownames(x)
  # Replace mutated base with dot
  substring(context,2,2) = "."
  # Construct dataframe for plotting
  df = data.frame(substitution = substitutions[index], context = context)
  rownames(x) = NULL
  df2 = cbind(df, as.data.frame(x))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  # Add dummy non_visible data points to force y axis limits per facet
  df4 = data.frame(substitution = rep("C>A", 4), context = rep("A.A",4), variable = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2]))
  
  plot = ggplot(data=df3, aes(x=context, y=value, fill=substitution, width=0.6)) +  
    geom_bar(stat="identity", position = "identity", colour="black", size=.2) + 
    geom_point(data = df4, aes(x = context, y = value), alpha = 0) +
    scale_fill_manual(values=colors) + 
    facet_grid(variable ~ substitution, scales = "free_y") + 
    ylab("Relative contribution") + 
    # ylim(-yrange, yrange) +
    # no legend
    guides(fill=FALSE) + 
    # white background
    theme_bw() +
    ggtitle(paste("RSS =", RSS)) + 
    # format text
    theme(axis.title.y=element_text(size=12,vjust=1),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=12),
          axis.text.x=element_text(size=5,angle=90,vjust=0.4),
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14),
          panel.grid.major.x = element_blank())
  return(plot)
}
