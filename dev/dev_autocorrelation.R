# Sigtools code tutorial -- autocorrelation

# first thing first, let's develop some methods to divide our tasks
# complete the following code
get_num_signals <- function(df){
  # calculates and returns the number of signals in a given multi-column bedgraph file
  # to complete this code section, note how many columns are there in a multi-column bedgraph file
  # and what are those columns
  # there is an example of a multi-column bedgraph file in './test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph'
}

import_signals <- function(path, header, sep, quote, prefix){
  df <- read.csv(path, header = header, sep = sep, quote = quote)
  num_signals <- get_num_signals(df[1, ])
  if(! header){
    colnames(df) <- c('chr', 'start', 'end', paste(prefix, 1:num_signals, sep=''))
  }
  return(df)
}
