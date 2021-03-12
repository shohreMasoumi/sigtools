get_num_signals <- function(df){
  return(ncol(df) - 3)
}

import_signals <- function(path, header, sep, quote, prefix){
  df <- read.csv(path, header = header, sep = sep, quote = quote)
  num_signals <- get_num_signals(df[1, ])
  if(! header){
    colnames(df) <- c('chr', 'start', 'end', paste(prefix, 1:num_signals, sep=''))
  }
  return(df)
}

select_intervals <- function(df_bedg, path_intervals, chr){
  df_regions <- get_regions(path_intervals, chr)
  for(r in 1:nrow(df_regions)){
    print(r)
    print(df_regions[r, ])
    start_ind <- which(df_bedg$start <= df_regions$start[r] & df_bedg$end   > df_regions$start[r])
    end_ind   <- which(df_bedg$start < df_regions$end[r]   & df_bedg$end   >= df_regions$end[r])
    df_select_temp <- fill_region(df_bedg[(start_ind - neighborhood):(end_ind + neighborhood),  (4:ncol(df_bedg))], neighborhood)
    df_select <- df_select + df_select_temp

  }
  return(df_select)
}

#' Plots the autocorrelation of several signals
#'
#' @param path_mulColBedG The path to the input file of multi-column bedgraph
#' @param outdir Name of the output directory [default:autocorrelation]
#' @param header Does your file include headers? [default:FALSE]
#' @param prefix Your preferred prefix for naming the signals [default:'s']
#' @param sep Your field separator character [default:'\t']
#' @param resolution The bin size of the input file [default:200]
#' @param lag The number of demanded shifts [default:50]
#' @param mode If set to 'regions', the autocorrelation is calculated over adjacent elements mentioned in the additional bed file [default:bins] [options:'bins' or 'regions']
#' @param path_intervals If given an interval file (.bed), the analysis would be carried out over the these intervals [default:NA]
#' @param chr The intended chromosome [default:NA]
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
#'@importFrom ggplot2 ggplot ggsave geom_line aes theme_classic ylab xlab theme labs element_blank
#'@importFrom parallel detectCores mclapply
#'
#' @examples
#' sigtools_autocorrelation("./test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph", header = TRUE, resolution = 200L, lag = 20L)
sigtools_autocorrelation <- function(path_mulColBedG,
                            outdir = "autocorrelation",
                            header = FALSE, prefix = 's', sep = '\t',
                            resolution = 200L, lag = 50L,
                            mode = 'bins', path_intervals = NA, chr=NA,
                            x_label = NA, y_label = NA, img_title = NA,
                            img_width = 10, img_height = 10, fontSize = 20){

  num_cores <- detectCores()
  df_mergedBedGraph_wide <- import_signals(path = path_mulColBedG, header = header, sep = sep, quote = "", prefix)
  if(! is.na(path_intervals)){
    df_mergedBedGraph_wide <- select_intervals(df_bedg, path_intervals, chr)
  }
  dir.create(outdir)
  df_autoCor <- read.csv('', col.names = c('lag', 'feature', 'autoCorrelation'),header = FALSE, sep = ' ')

  for(i in 1:get_num_signals(df_mergedBedGraph_wide[1, ])){
    #signal_name <- paste(prefix, 1:num_signals, sep='')[i]
    signal_name <- colnames(df_mergedBedGraph_wide)[i+3]
    print(signal_name)
    df_temp <- data.frame('lag' = (1:lag)*resolution,
                          'feature' = rep(signal_name, lag),
                          'autoCorrelation' = rep(0, lag))
    #head(df_temp)
    for(l in 1:lag){
      #print(l)
      #l <- 1
      list_splited <- split(df_mergedBedGraph_wide, f = df_mergedBedGraph_wide$chr)
      x <- c()
      y <- c()
      for(df in list_splited){
        x <- c(x, df[(l:nrow(df)), i+3])
        y <- c(y, df[1:(nrow(df)-l+1), i+3])
      }
      df_temp$autoCorrelation[l] <- cor(x, y, method = c("pearson"))
    }
    df_autoCor <- rbind(df_autoCor, df_temp)
    #break
  }

  cat('Generating autocorrelation plot...\n')
  p <- ggplot(df_autoCor, aes(lag, autoCorrelation, group = feature, color= feature )) +
    geom_line() +
    theme(legend.title=element_blank()) +
    labs(fill = "", x = "distance (bp)", y ="autocorrelation") +
    #scale_color_manual(values=myColors) +
    #scale_color_discrete(name="") +
    theme_classic(base_size = fontSize) +
    theme(legend.title=element_blank())

  if(! is.na(x_label)){p <- p + xlab(x_label)}
  if(! is.na(y_label)){p <- p + ylab(y_label)}


  ggsave(paste0('./', outdir, '/', 'autoCorrelation.png'),
         plot = p,
         #scale = 0.7,
         #width = 5.8, height = 4.13 , units = 'in',
         width = img_width, height = img_height, units = "cm",
         dpi = 300, limitsize = TRUE)
}
