#' Plot point mutation spectrum
#'    
#' @param type_occurences Type occurences matrix
#' @param CT Distinction between C>T at CpG and C>T at other sites, default = F
#' @param by Optional grouping variable
#' @param colors Optional color vector with 7 values
#' @param legend Plot legend, default = T
#' @return Spectrum plot
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom BiocGenerics cbind
#' @importFrom plyr ddply
#' @importFrom plyr summarise
#' @export

plot_spectrum = function(type_occurences, CT = F, by, colors, legend = T)
{
  # if colors parameter not provided, set to default colors
  if(missing(colors)){colors = COLORS7}
  # check color vector length
  if(length(colors) != 7){stop("Colors parameter: supply color vector with length 7")}
  # distinction between C>T at CpG or not
  if(CT == FALSE){type_occurences = type_occurences[,1:6] }
  if(CT == TRUE){type_occurences = type_occurences[,c(1:2,8,7,4:6)]}
  # relative contribution per sample
  df2 = type_occurences / rowSums(type_occurences)
  # if grouping variable not provided, set to "all"
  if(missing(by)){by="all"}
  # add by info to df
  df2$by = by
  # reshape
  df3 = melt(df2, id.vars = "by")
  # count number of mutations per mutation type
  counts = melt(type_occurences, measure.vars = colnames(type_occurences))
  df4 = cbind(df3, counts$value)
  colnames(df4)[4] = "nmuts" 
  # calculate the mean and the stdev on the value for each group broken down by type + variable
  x = ddply(df4, c("by", "variable"), summarise, mean = mean(value), stdev = sd(value))
  info_x = ddply(df4, c("by"), summarise, total_individuals = sum(value), total_mutations = sum(nmuts))
  x = merge(x, info_x)
  info_type = data.frame(sub_type = c("C>A", "C>G", "C>T", "C>T", "C>T", "T>A", "T>C", "T>G"), variable = c("C>A", "C>G", "C>T", "C>T at CpG", "C>T other", "T>A", "T>C", "T>G"))
  x = merge(x,info_type)
  x$total_mutations = prettyNum(x$total_mutations, big.mark = ",")
  x$total_mutations = paste("N =", as.character(x$total_mutations))
  # define colors for plotting
  if(CT == F){colors = colors[c(1,2,3,5:7)]}
  # define positioning of error bars
  x$error_pos = x$mean
  # if C>T stacked bar (distinction between CpG sites and other)
  if(CT == T){
    # adjust positioning of error bars for stacked bars
    # mean of C>T at CpG should be plus the mean of C>T at other
    x = x[order(x$by),]
    CpG = which(x$variable == "C>T at CpG")
    other = which(x$variable == "C>T other")
    x$error_pos[CpG] = x$error_pos[other] + x$error_pos[CpG]
    # define order of bars
    order = order(factor(x$variable, levels = c("C>A", "C>G",  "C>T other", "C>T at CpG","T>A", "T>C", "T>G")))
    x = x[order,]
  }
  # make barplot
  plot = ggplot(data=x, aes(x=sub_type, y=mean, fill=variable, group=sub_type)) +
          geom_bar(stat="identity") +
          geom_errorbar(aes(ymin=error_pos-stdev, ymax=error_pos+stdev), width=0.2) + 
          scale_fill_manual(values=colors, name="Point mutation type") +
          theme_bw() +
          xlab("") +
          ylab("Relative contribution") + 
          theme(axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank()) 
  # facetting
  if(length(by) == 1){plot = plot + facet_wrap( ~ total_mutations)}
  else{plot = plot + facet_wrap(by ~ total_mutations)}
  # legend
  if(legend == FALSE){plot = plot + theme(legend.position="none")}
  return(plot)
} 
