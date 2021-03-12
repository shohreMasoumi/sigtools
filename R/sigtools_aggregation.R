get_signal_num <- function(df_bedg){
  return(num_signals <- ncol(df_bed) - 3)
}

get_signals <- function(path, header, prefix, sep){
  print("Loading signals")
  df <- read.csv(path, header = header, sep = sep, quote = "")
  num_signals <- ncol(df) - 3
  if(! header){
    colnames(df) <- c('chr', 'start', 'end', paste(prefix, 1:num_signals, sep=''))
  }
  return(df)
}

get_enriched <- function(df_bedg, means){
  print("Converting to enrichment")
  num_signals <- ncol(df_bedg) - 3
  for(i in 1:num_signals){
    df_bedg[ ,(3+i)] <- df_bedg[ ,(3+i)]/means[i]
  }
  return(df_bedg)
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

get_aggregation_frame <- function(nei, df_sig){
  df_agg <- data.frame('position' =  rep((-(nei):(nei)), 2))
  for(i in 1:ncol(df_sig)){df_agg <- cbind(df_agg, rep(0, (nei*2+1)*2 ))}
  colnames(df_agg) = c('position', colnames(df_sig))
  df_agg <- df_agg[ , 2:(ncol(df_agg))]
  return(df_agg)
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

fill_body <- function(values, nei){
  coords <- data.frame('x' = 1:length(values), 'y' = values)
  # https://stackoverflow.com/questions/11356997/interpolating-a-path-curve-within-r
  df_new_coords <- as.data.frame(with(coords, aspline(x = x, y = y, n = (2*nei + 2))))
  return(df_new_coords$y)
}

flip_aggregation_frame <- function(df_agg){
  print("Flipping")
  for(j in 1:ncol(df_agg)){
    df_agg[ ,j] <- df_agg[nrow(df_agg):1, j]
  }
  return(df_agg)
}

save_plot <- function(p, outdir, name){
  ggsave(paste0(outdir, '/aggregation_', name,'.png'),
         plot = p,
         #scale = 0.7,
         width = 6, height = 5 , units = 'in',
         dpi = 300, limitsize = TRUE)
}

create_aggregation_plots <- function(df_agg, nei, res, outdir){
  position <- (((-2*nei-1):(2*nei))*res)
  mid_line <- 1
  for(i in 1:ncol(df_agg)){
    print(paste0("Plotting ", colnames(df_agg)[i]))
    signal <- df_agg[ ,i]
    plot <- ggplot(, aes(x = position, y = signal)) +
      #ggplot(agg, aes(x=position,y=signal)) +
      geom_ribbon(aes(ymin = pmin(signal, mid_line), ymax = mid_line), fill="#00BFC4", alpha=0.9) +
      geom_ribbon(aes(ymin = mid_line, ymax=pmax(signal, mid_line)),   fill="#F8766D", alpha=0.9) +
      geom_line(aes(y = mid_line)) +
      ylab(colnames(df_agg)[i]) +
      xlab("") +
      geom_vline(xintercept = c(-nei*res, nei*res), linetype="dashed") +
      #geom_vline(xintercept = 32, linetype="dashed") +
      theme_classic(base_size = 21.5)

    #plot
    save_plot(plot, outdir, colnames(df_agg)[i])
  }
}

#' Generates aggregation plots for all the signals in the given multi-column bedGraph file, over the specified elements in the given interval file
#'
#' @param path_mulColBedG The path to the input multi-column bedgraph file
#' @param path_regThe The path to the inpput interval file
#' @param chr The intended chromosome (example:21L)
#' @param mode Controls the length unification process [options: 'point' or 'stretch']
#' @param header Does your file include headers? [default:FALSE]
#' @param prefix Your preferred prefix for naming the signals [default:'s']
#' @param sep Your field separator character [default:'\t']
#' @param outdir Name of the output directory [default:aggregation]
#' @param resolution The bin size of the input file [default:200]
#' @param neighborhood The number of covered up and down stream bins [default:50]
#' @param percentile What percentile of your data you want to work with [default:'100']
#' @param enrichment Do you like to plot the enrichment? [default:FALSE]
#' @param means If given, istead of calculating mean values of the signals, the mean will be obtained from this list [default:c(NA)]
#'
#' @return
#' @export
#'
#' @importFrom ggplot2 ggplot geom_ribbon geom_line aes ylab xlab theme_classic geom_vline ggsave
#' @importFrom reshape2 melt
#' @importFrom parallel detectCores mclapply
#' @importFrom akima aspline
#'
#' @examples
#' sigtools_aggregation(path_mulColBedG = "./test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph", path_reg = "./test_data/Ensembl_chr21.gene_info", chr = 21, mode = 'stretch', header = TRUE , outdir = 'aggregation_E003-chr21_enriched_2', resolution = 200L , neighborhood = 10L, enrichment = TRUE, means = c(0.640, 0.601, 0.633, 0.575, 0.604, 0.597) )
sigtools_aggregation <- function(path_mulColBedG, path_reg, chr, mode,
                        header = FALSE, prefix = 's', sep = '\t', outdir = "aggregation",
                        resolution = 200L, neighborhood = 50L,
                        percentile = 100, enrichment = FALSE, means = c(NA)){

  # Loading signals
  df_bedg <- get_signals(path_mulColBedG, header, prefix, sep)
  if(enrichment){df_bedg <- get_enriched(df_bedg, means)}

  # Loading regions
  df_regions <- get_regions(path_reg, chr)

  if(mode == 'stretch'){
    df_agg <- get_aggregation_frame(neighborhood, df_bedg[1, (4:ncol(df_bedg))])
    df_agg
    r <- 1
    for(r in 1:nrow(df_regions)){
      print(r)
      print(df_regions[r, ])
      start_ind <- which(df_bedg$start <= df_regions$start[r] & df_bedg$end   > df_regions$start[r])
      end_ind   <- which(df_bedg$start < df_regions$end[r]   & df_bedg$end   >= df_regions$end[r])
      df_agg_temp <- fill_region(df_bedg[(start_ind - neighborhood):(end_ind + neighborhood),  (4:ncol(df_bedg))], neighborhood)
      if(df_regions$direction[r] == '-1'){
        df_agg_temp <- flip_aggregation_frame(df_agg_temp)
      }
      df_agg <- df_agg + df_agg_temp
      df_agg
    }
    #saving <- df_agg
    df_agg
    for(j in 1:ncol(df_agg)){df_agg[ ,j] <- df_agg[ ,j]/nrow(df_regions)}
  } # END if(mode == 'stretch')

  dir.create(outdir)
  create_aggregation_plots(df_agg, neighborhood, resolution, outdir)
}
