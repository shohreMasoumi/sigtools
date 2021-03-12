#' Plots a heatmap of pair-wise correlation of sets of genomic signals
#'
#' @param path_mulColBedG_1 The path to the file of multi-column bedgraph to be considered as the first set of signals
#' @param header_1 Does your file include headers? [default:FALSE]
#' @param prefix_1 Your preferred prefix for naming the signals in the first set [default:'set1_']
#' @param path_mulColBedG_2 The path to the file of multi-column bedgraph to be considered as the second set of signals
#' @param header_2 Does your file include headers? [default:FALSE]
#' @param prefix_2 Your preferred prefix for naming the signals in the second set [default:'set2_']
#' @param outdir Name of the output directory [default:correlation]
#' @param sep Your field separator character [default:'\t']
#' @param enrichment Do you like to plot the enrichment? [default:FALSE]
#' @param x_label Label for the X axis
#' @param y_label label for the Y axis
#' @param img_title Title of the graph
#' @param img_width Width of the saved image (cm) [default:10]
#' @param img_height Height of the saved image (cm) [default:10]
#' @param fontSize The font size of the saved image [default:20]
#'
#' @return
#' @export
#'
#' @importFrom ggplot2 geom_tile labs xlab ylab
#' @importFrom viridis scale_fill_viridis
#' @importFrom reshape2 melt
#' @importFrom ggcorrplot ggcorrplot
#'
#' @examples
#' sigtools_correlation(path_mulColBedG_1 = "./test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph", header_1 = TRUE, path_mulColBedG_2 = "./test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph", header_2 = TRUE, outdir = "correlation", img_width = 15, img_height = 15)
sigtools_correlation <- function(path_mulColBedG_1, header_1 = FALSE, prefix_1 = 'set1',
                        path_mulColBedG_2, header_2 = FALSE, prefix_2 = 'set2',
                        outdir = "correlation",
                        sep = '\t',
                        enrichment = FALSE,
                        x_label = NA, y_label = NA, img_title = NA,
                        img_width = 10, img_height = 10, fontSize = 20){

  df_mcbg_1 <- read.csv(path_mulColBedG_1, header = header_1, sep = sep, quote = "")
  num_signals_1 <- ncol(df_mcbg_1) - 3
  if(! header_1){
    colnames(df_mcbg_1) = c('chr', 'start', 'end', paste(prefix_1, 1:num_signals_1, sep=''))
  }
  df_mcbg_1 <- df_mcbg_1[ ,4:(3+num_signals_1)]


  df_mcbg_2 <- read.csv(path_mulColBedG_2, header = header_2, sep = sep, quote = "")
  num_signals_2 <- ncol(df_mcbg_2) - 3
  if(! header_2){
    colnames(df_mcbg_2) = c('chr', 'start', 'end', paste(prefix_2, 1:num_signals_2, sep=''))
  }
  df_mcbg_2 <- df_mcbg_2[ ,4:(3+num_signals_2)]

  if(enrichment){
    for(i in 1:num_signals_1){
        mean <- mean(df_mcbg_1[ ,i])
        df_mcbg_1[ ,i] <-  df_mcbg_1[ ,i]/mean
    }
    for(i in 1:num_signals_2){
      mean <- mean(df_mcbg_2[ ,i])
      df_mcbg_2[ ,i] <-  df_mcbg_2[ ,i]/mean
    }
    rm(mean)
  }

  # chop the extra length
  l <- min(nrow(df_mcbg_1), nrow(df_mcbg_2))
  df_mcbg_1 <- df_mcbg_1[1:l, ]
  df_mcbg_2 <- df_mcbg_2[1:l, ]

  #head(df_mergedBedGraph_1_wide)
  mat_correlation <- matrix(0, nrow = num_signals_1, ncol = num_signals_2, byrow = TRUE)
  rownames(mat_correlation) <- names(df_mcbg_1)
  colnames(mat_correlation) <- names(df_mcbg_2)
  for (i in 1:num_signals_1) {
    for (j in 1:num_signals_2) {
      mat_correlation[i,j] <- cor(df_mcbg_1[i], df_mcbg_2[j], method = c("pearson"))
    }
  }

  # JUL 9 2020
  # TODO: insead of pairwise correlation calculation: df <- rbind(df_1, df_2), cor(df)
  #mat_correlation
#  ggcorrplot(mat_correlation)
#  ggcorrplot(mat_correlation, method = "circle")
#  ggcorrplot(mat_correlation, method = "circle", hc.order = TRUE)
#  ggcorrplot(mat_correlation, method = "circle", hc.order = TRUE, outline.col = "white")
#  ggcorrplot(mat_correlation, type = "lower",    hc.order = TRUE, outline.col = "white")
#  ggcorrplot(mat_correlation, type = "lower",    hc.order = TRUE)
#  ggcorrplot(mat_correlation, type = "lower")
#  ggcorrplot(mat_correlation)
#  ggcorrplot(mat_correlation, lab = TRUE)
#  ggcorrplot(mat_correlation, lab = TRUE, type = "lower")


  # commented on Jul 10 *********************************
#  mat_correlation_melted <- melt(mat_correlation)
#  head(mat_correlation_melted)
#  colnames(mat_correlation_melted) <- c('set1', 'set2', 'value')

  #mat_correlation_flattened <- wideToLong(mat_correlation, percentile = 100, enrichment = enrichment, nozero = FALSE)
  #mat_correlation_flattened
  #head(mat_correlation_flattened)

#  cat('Generating correlation plot...\n')
#  p <- ggplot(mat_correlation_melted, aes(set1, set2, fill= value)) +
#    geom_tile() +
#    labs(fill = "")+
#    #scale_fill_gradient(low="white", high="#CC0633") +
#    #theme(axis.title.x=) +
#    xlab(element_blank()) +
#    ylab(element_blank()) +
#    scale_fill_viridis(discrete=FALSE) +
#    theme_classic(base_size = fontSize)

#  if(! is.na(x_label)){p <- p + xlab(x_label)}
#  if(! is.na(y_label)){p <- p + ylab(y_label)}

  # commented on Jul 10 *********************************

  cat('Generating correlation plot...\n')
  p <- ggcorrplot(mat_correlation, method = "circle")

  dir.create(outdir)

  ggsave(paste0('./', outdir, '/','heatmap.png'),
         plot = p,
         #scale = 0.7,
         #width = 9.0, height = 4.5 , units = 'in',
         width = img_width, height = img_height, units = "cm",
         #width = 30, height = 7.5, units = "cm",
         dpi = 300, limitsize = TRUE)
}



