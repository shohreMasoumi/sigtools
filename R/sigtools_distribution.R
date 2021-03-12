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

get_regions <- function(path, chr){
  print("Loading Regions")
  df <- read.csv(path, header = FALSE, sep = '\t', quote = "")[ ,1:5]

  format <- strsplit(basename(path), split="\\.")[[1]][-1]
  if( format == "bed"){
    num <- ncol(df)
    colnames(df) <- c('chr', 'start', 'end', 'name', 'value')
  } else if(format == "gene_info"){
    colnames(df) <- c('name', 'chr', 'start', 'end', 'direction')
  } else{
    print("Sigtools current region formats are .bed and .gene_info")
    print("Unidentified file format")
    print("exiting")
    quit()
  }
  df <- df[df$chr == chr, ]
  df$center <- as.integer((df$start + df$end  )/2)
  df$width  <- as.integer( df$end   - df$start)

  #  for(r in 1:nrow(df)){
  #    if(df$direction[r] == '-1'){
  #      temp <- df$start[r]
  #      df$start[r] <- df$end[r]
  #      df$end[r] <- temp
  #    }
  #  }
  return(df)
}

fill_region <- function(df_sig, nei){
  df_agg <- get_aggregation_frame(nei, df_sig[1, ])
  df_agg[1:nei, ] <- df_sig[1:nei, ]
  #print(nrow(df_sig))
  #print((nrow(df_sig) - nei + 1):nrow(df_sig))
  df_agg[(nrow(df_agg) - nei + 1):nrow(df_agg), ] <-  df_sig[(nrow(df_sig) - nei + 1):nrow(df_sig), ]
  for(i in 1:ncol(df_sig)){
    df_agg[(nei + 1):(nrow(df_agg) - nei), i] <- fill_body(df_sig[ ,i], nei)
  }
  return(df_agg)
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

exclude_blacklist <- function(df_bedg, path_intervals, chr){
  df_regions <- get_regions(path_intervals, chr)
  #myData[-c(2, 4, 6), ]
  for(r in 1:nrow(df_regions)){
    print(r)
    print(df_regions[r, ])
    start_ind <- which(df_bedg$start <= df_regions$start[r] & df_bedg$end   > df_regions$start[r])
    end_ind   <- which(df_bedg$start < df_regions$end[r]   & df_bedg$end   >= df_regions$end[r])
    df_bedg <- df_bedg[-c(start_ind:end_ind), ]
  }
  return(df_bedg)
}

flatten_mulColBedG <- function(df, percentile, enrichment, nozero, num_cores){
  num_signals <- get_num_signals(df[1, ])
  list_signals <- list()

  if(num_signals == 1){
    df  <- data.frame('signal' = rep(names(df)[4], length(df)),
                      'value' =  df[ ,4])
  }
  else{
    df <-  df[ ,(4:(3+num_signals))]
    list_signals <- mclapply(df, FUN= function(col, per, en, nozero){
      col <- col[!is.na(col)]
      print(head(col))
      if(per != 100){ col <- col[col < quantile(col, c(per*0.01))]}
      if(en)        { col <- col/mean(col, na.rm=TRUE)}
      if(nozero)    {col <- col[which(col!= 0)]}
      return(col)
    }, percentile, enrichment, nozero, mc.cores = num_cores)

    df  <- read.csv('', col.names = c('signal', 'value'),header = FALSE, sep = ' ')
    for(i in 1:length(list_signals)){
      df <- rbind(df, data.frame('signal' = rep(names(list_signals)[i], length(list_signals[[i]])),
                                 'value' = list_signals[[i]]))
    }
  }
  return(df)
}

save_plot <- function(p, outdir, name, img_width, img_height, scale){
  ggsave(paste0('./', outdir, '/', name),
         plot = p,
         scale = scale,
         width = img_width, height = img_height, units = "cm",
         dpi = 300, limitsize = TRUE)
}

generate_density_curve <- function(df, fontSize, x_label, y_label){
  cat("Generating distribution curve...\n")
  p <- ggplot(df, aes(x=value, y=signal, fill=signal)) +
       geom_density_ridges2(scale=0.9) + # rel_min_height = 0.01
       #labs(colour = "") +
       #theme(legend.title=element_blank()) +
       #labs(fill = "", x= 'mean enrichment' , y= 'feature' )+
       theme_classic(base_size = fontSize) +
      theme(legend.title=element_blank()) +
      theme(legend.position = "none")

  if(! is.na(x_label)){p <- p + xlab(x_label)}
  if(! is.na(y_label)){p <- p + ylab(y_label)}
  return(p)
}

generate_ecdf <- function(df, fontSize, x_label, y_label){
  cat("Generating cumulative distribution...\n")
  p <- ggplot(df, aes(x=value, colour=signal)) +
       stat_ecdf(geom = "step", size = 0.8) +
       #xlab("value") +
       #ylab("empirical cumulative of \nmean enrichment") +
       ylab("empirical cumulative \n distribution") +
       labs(fill = "") +
       theme_classic(base_size = fontSize) +
       theme(legend.title=element_blank()) #legend.position = c(.9, .25)

  if(! is.na(x_label)){p <- p + xlab(x_label)}
  if(! is.na(y_label)){p <- p + ylab(y_label)}
  return(p)
}

generate_boxplot <- function(df, fontSize, x_label, y_label){
  cat('Generating boxplot...\n')
  p <- ggplot(df, aes(x=signal, y=value, fill=signal)) +
    geom_boxplot(alpha=1) + #alpha=0.2
    #labs(fill = "", x= 'Signals' ) + #  , y= ''
    #stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="red", fill="red") +
    theme_classic(base_size = fontSize)

  if(! is.na(x_label)){p <- p + xlab(x_label)}
  if(! is.na(y_label)){p <- p + ylab(y_label)}
  return(p)
}

generate_boxplot <- function(df, fontSize, x_label, y_label){
  cat('Generating histogram...\n')
  p <- ggplot(df, aes(x=value, y=feature, fill=feature)) +
       geom_density_ridges(alpha= 0.8, stat="binline", scale = 0.9, rel_min_height = 0.01) +
       theme_ridges() +
       theme(
         legend.position="none",
         panel.spacing = unit(0.1, "lines"),
         #strip.text.x = element_text(size = 8)
       ) +
       theme_classic(base_size = fontSize)

  if(! is.na(x_label)){p <- p + xlab(x_label)}
  if(! is.na(y_label)){p <- p + ylab(y_label)}
}

#' Plots the distribution of several signals
#'
#' @param path_mulColBedG The path to the input file of multi-column bedgraph to depict distribution for
#' @param outdir Name of the output directory [default:distribution]
#' @param plots Options of distribution plots [default:plots = c("curve", "ecdf", "box")]
#' @param header Does your file include headers? [default:FALSE]
#' @param prefix Your preferred prefix for naming the signals [default:'s']
#' @param sep Your field separator character [default:'\t']
#' @param percentile What percentile of your data you want to work with [default:'100']
#' @param enrichment Do you like to plot the enrichment? [default:FALSE]
#' @param nozero Do you want to remove the zero entries? [default:FALSE]
#' @param path_intervals If given an interval file (.bed), the analysis would be carried out over the these intervals [default:NA]
#' @param chr The intended chromosome [default:NA]
#' @param path_blacklist If given, the intervals in this file would be excluded from the analysis [default:NA]
#' @param x_label Label for the X axis
#' @param y_label label for the Y axis
#' @param img_title Title of the graph
#' @param img_width Width of the saved image (cm) [default:10]
#' @param img_height Height of the saved image (cm) [default:10]
#' @param fontSize The font size of the saved image [default:20]
#' @param scale [default:0.7]
#'
#' @return
#' @export
#'
#' @importFrom ggplot2 ggplot ggsave stat_ecdf geom_boxplot aes theme_classic last_plot ylab xlab theme labs
#' @importFrom ggplot2 element_blank
#' @importFrom ggridges geom_density_ridges2 geom_density_ridges
#' @importFrom parallel detectCores mclapply
#'
#' @examples
#' sigtools_distribution("./test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph", plots = c("curve"), header = TRUE, percentile = 99, enrichment = FALSE, nozero = TRUE, img_height = 10, img_width = 20, y_label = "")
sigtools_distribution <- function(path_mulColBedG,
                         outdir = "distribution", plots = c("curve", "ecdf", "box"),
                         header = FALSE, prefix = 's', sep = '\t',
                         percentile = 100, enrichment = FALSE, nozero = FALSE,
                         path_intervals = NA, chr=NA, path_blacklist = NA,
                         x_label = NA, y_label = NA, img_title = NA,
                         img_width = 10, img_height = 10, fontSize = 20, scale = 0.7){

  num_cores <- detectCores()
  df_signals <- import_signals(path = path_mulColBedG, header = header, sep = sep, quote = "", prefix)
  if(! is.na(path_intervals)){
    df_signals <- select_intervals(df_bedg, path_intervals, chr)
  }
  df_signals <-  flatten_mulColBedG(df_signals, percentile, enrichment, nozero, num_cores)
  dir.create(outdir)

  if("curve" %in% plots){
    p <- generate_density_curve(df_signals, fontSize, x_label, y_label)
    save_plot(p=p, outdir=outdir, name='distributionCurve.png', img_width=img_width, img_height=img_height, scale=scale)

    # if(! is.na(legend_title)){ p <- p + theme(legend.title = legend_title)}
    #
    # ggplot(df_signals, aes(x=signal, y=value, fill=signal)) + # fill=name allow to automatically dedicate a color for each group
    #   geom_violin() +
    #   theme_classic(base_size = 21)
    #
    # ggplot(df_signals, aes(x=feature, y=value, fill=feature, color=feature)) +
    #   geom_violin(width=2.1, size=0.2) +
    #   scale_fill_viridis(discrete=TRUE) +
    #   scale_color_viridis(discrete=TRUE)
  }

  if("ecdf" %in% plots){
    p <- generate_ecdf(df_signals, fontSize, x_label, y_label)
    save_plot(p=p, outdir=outdir, name= 'ecdf.png', img_width=img_width, img_height=img_height, scale=scale)
  }

  if('box' %in% plots){
    p <- generate_boxplot(df_signals, fontSize, x_label, y_label)
    save_plot(p=p, outdir=outdir, name= 'boxplot.png', img_width=img_width, img_height=img_height, scale=scale)
  }

  if('histogram' %in% plots){
    p <- generate_boxplot <- function(df_signals, fontSize, x_label, y_label)
    save_plot(p=p, outdir=outdir, name= 'histogram.png', img_width=img_width, img_height=img_height, scale=scale)
  }
}
