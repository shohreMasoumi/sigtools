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

library(ggplot2)
library(parallel)

sigtools_autocorrelation(path_mulColBedG = "./test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph",
                        outdir = "autocorrelation",
                        header = TRUE,
                        prefix = 's', sep = '\t',
                        resolution = 200L,
                        lag = 20L,
                        mode = 'bins', path_intervals = NA, chr=NA,
                        x_label = NA, y_label = NA, img_title = NA,
                        img_width = 10, img_height = 10, fontSize = 20)
