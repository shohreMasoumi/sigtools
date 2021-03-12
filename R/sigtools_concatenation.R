#' Title
#'
#' @param path The path to the directory containing your cell types
#' @param pattern The pattern in which you named your cell types
#' @param num_signals The number of signals presented in the dataset
#' @param header Does this file include headers? [default: FALSE]
#' @param sep The field separator character [default: '\t']
#' @param outdir The name of the output file [default: 'preprocessing']
#' @param sampling Wheather you want the input file(s) to be sampled? [default: TRUE]
#' @param population_per The size of the sample relative to the size of the dataset  [default: 0.01]
#' @param population_ratio The number of sampled intervals relative to their length [default: c(1,2)]
#'
#' @return
#' @export
#'
#' @examples
#' sigtools_concatenate('../compCan/ssmEdited', "exp73-7-E*_edited", 3)
sigtools_concatenate <- function(path, pattern, num_signals,
                                 header = FALSE, sep = '\t', outdir = "./preprocessing",
                                 sampling = TRUE, population_per = 0.1, population_ratio = c(1,2)){
  # Creating the output directory
  dir.create(outdir)

  # Setting the path to input files
  path_cellTypes <- Sys.glob(paste(path, pattern, sep='/'))

  # Setting the name of the new file
  if(sampling){
    output_name <- paste0(gsub('*', "STAR", pattern), '_sampled_concatenated.mulColBedg')
  } else{
    output_name <- paste0(gsub('*', "STAR", pattern), '_concatenated.mulColBedg')
  }


  # Creating a empty data frame
  df_cellTypeAgnostic <- read.csv('', header = header, sep = sep,
                                  col.names = c('chr', 'start', 'end', paste('s', 1:num_signals, sep='')))


  # stating the timer
  start_time <- Sys.time()
  for(i in 1:length(path_cellTypes)){
    print(paste0("sampling: ", sampling,
                 "  population_per: ", population_per,
                 "  population_ratio: ", population_ratio))

    temp <- con_chr(path_cellTypes[i], num_signals, sep = sep,
                   sampling = sampling , population_per = population_per, population_ratio = population_ratio)

    write.table(temp, paste(outdir, output_name, sep = "/"), append = TRUE,
                quote = FALSE, sep = '\t', col.names = FALSE, row.names = FALSE)

  }
  end_time <- Sys.time()
  print(end_time - start_time)
}


con_chr <- function(path, num_signals, sep = '\t',
                    sampling , population_per, population_ratio){

  files_chr <- list.files(path)
  df_concatenated <- read.csv('', header = header, sep = sep,
                              col.names = c('chr', 'start', 'end', paste('s', 1:num_signals, sep='')))

  for(f in files_chr){
    print(f)
    temp <- read.csv(paste(path, f, sep = '/'), header = FALSE, sep = sep)
    if(sampling){
      temp <- sampling(temp, pop_per = population_per, pop_ratio = population_ratio)
    }
    df_concatenated <- rbind(df_concatenated, temp)
  }

  print(object.size(df_concatenated), units="Mb")
  return(df_concatenated)
}


sampling <- function(df_signal, pop_per = 0.01, pop_ratio = NA){

  num_signals <- ncol(df_signal) - 3
  population_line <- ceiling(pop_per * nrow(df_signal))

  # num_head * num_neigh = pop_lines
  # num_head = (ratio[2]/ratio[1]) num_neigh
  # --> (ratio[2]/ratio[1]) num_neigh * num_neigh = pop_lines
  # num_neigh^2 = pop_lines * ratio[1]/ratio[2]

  #num_neigh <- sqrt(population_line*(ratio[1]/ratio[2]))
  #num_head = (ratio[2]/ratio[1]) * num_neigh

  if(is.na(pop_ratio)){
    num_neighbors <- round(sqrt(population_line*(pop_ratio[1]/pop_ratio[2])))
    num_heads <-  round((pop_ratio[2]/pop_ratio[1]) * num_neighbors)
  } else {
    num_heads <- 1
    num_neighbors <- population_line
  }

  heads <- sample(nrow(df_signal), num_heads, replace = FALSE)
  df_sample <- read.csv('', header = FALSE, sep = "\t", col.names = c('chr', 'start', 'end', paste('s', 1:num_signals, sep='')))
  for(i in 1:num_heads){
    df_sample <- rbind(df_sample, df_signal[heads[i]:(heads[i]+num_neighbors), ])
  }
  return(df_sample)
}
