# Sigtools code tutorial -- autocorrelation

# first thing first, let's develop some methods to divide our tasks
# complete the following code
get_num_signals <- function(df){
  # calculates and returns the number of signals in a given multi-column bedgraph file
  # to complete this code section, note how many columns are there in a multi-column bedgraph file
  # and what are those columns
  # there is an example of a multi-column bedgraph file in './test_data/E003-assays_chr21_bin200_bwtool.mcBedGraph'

}

# we already did this in the previous script, but for the sake of practice, let's do it again
import_signals <- function(path, header, sep, quote, prefix){
  # read a dataframe from the `path` with other input parameters
  # if the file has a header, leave it be, if not,  rename the coloumns to
  # c('chr', 'start', 'end', 'prefix1', 'prefix2', ..., 'prefixN')
  # return the dataframe
  return(df)
}

# An interval file is a file that contains specific stretches of the genome or chromosome (see an example in './test_data/E003_activeGenes.gene_info')
# Imagine we imported a dataframe with the `import_signals` method
# and we also have an interval file (path_intervals),
# and we want to separate the DNA stretches mentioned in the interval file from the dataframe
select_intervals <- function(df_bedg, path_intervals, chr){
  # takes a dataframe, and a path to an interval file, and a chromosome number
  # returns a dataframe with the rows from the original dataframe that sit inside the mentioned intervals in the interval file

  # hints :
  #df_regions <- #import the data from path_intervals
  #for(r in 1:nrow(df_regions)){
  #}

  return(df_select)
}

# The main method,
sigtools_autocorrelation <- function(path_mulColBedG,
                                     outdir = "autocorrelation",
                                     header = FALSE, prefix = 's', sep = '\t',
                                     resolution = 200L, lag = 50L,
                                     mode = 'bins', path_intervals = NA, chr=NA,
                                     x_label = NA, y_label = NA, img_title = NA,
                                     img_width = 10, img_height = 10, fontSize = 20){
  # Takes several parameters such as
  # the path to a multiColumnBedGraph file

  # First thing first, read a bit on the autocorrelation,
  # comment here what the `lag` parameter is
  # contact me for more info

  # import the data
  #df_mergedBedGraph_wide <- import_signals(path = ?, header = ?, sep = ?, quote = "", prefix=?)

  # the path-interval is an optional parameter for this method,
  # if it is NA continue,
  # if in NOT is use the select_interval method to obtain a new df
  #if(! is.na(path_intervals)){
  #  df_mergedBedGraph_wide <- select_intervals(?, ?, ?)
  #}

  # create the output directory using outdir as its name

  # create an empty dataframe with columns: c('lag', 'feature', 'autoCorrelation')
  # this is actually going to be our final df that we use to generate the image from
  #df_autoCor <-

  # for every column in df_mergedBedGraph_wide (for every signal)
  #for(i in 1:?){
    #signal_name <-

    # Here, I am just creating a temp dataframe
    # to store the data for this iteration
    # Our main challenge in Sigtools is to process and configure data while being mindful that we have to use specific formats for `ggplot`
    # we will end up appending this df to our main df (df_autoCor) at the end of each round
    # Ho many rows should this df have?
    # could you figure out why I am doing this this way?
    #df_temp <- data.frame('lag' = (1:lag)*resolution,
    #                      'feature' = rep(signal_name, lag),
    #                      'autoCorrelation' = rep(0, lag))

    # R has a method for calculating correlation,
    # so, for each lag, we are going to calculate the correlation of the signal with itself but shifted (lagged) `l` times
    # first we need to format our data in a way to be accepted by the `cor` function
    #for(l in 1:lag){
      # in your df_mergedBedGraph_wide, different chromosomes might be included.
      # here we will devide our original df base on the chr column
      # if we oly have one chromosome, we would have a list with only on element, a df, our original df
      #list_splited <- split(df_mergedBedGraph_wide, f = df_mergedBedGraph_wide$chr)

      # the `corr` function takes three parameters:
      # x and y, together they form a set of points on a two axis coordinate plane
      # example: x=c(1, 2, 4, 6) and y=c(1, -1, 3, 12) for the set {(1,1), (2,-1), (4, 3), (6, 12)}
      # here, we will just create two empty lists for x and y variables
      # keep in mind that for `corr` to accept x and y, they need to be of the same length
      #x <- c()
      #y <- c()

      # for every separated df we obtained in list_splited

      #for(df in list_splited){
        # here is a bit tricky,
        # This whole for was to eliminate rare error in the case that
        # so we we are only working with one chromosome, this loop just runs 1
        # again, keep in mind that for `corr` to accept x and y, they need to be of the same length
        # I leave some section of the code to be,
        # for the sake of practice
        # if you are confident you are familiar with the code, just remove them and try to rewrite them yourself
        # lest take a look at:
        #x <- c(x, df[(l:nrow(df)), i+3])
        # here x is already an empty list,
        # the purpose of x <- c(x, ...) is to making sure x stays the same format (a list)
        # I know, it might not be the cleanest way to do this
        # in df[(l:nrow(df)), i+3]:
        # for the column choice `i` is the index of a signal, and because we have three other columns ('chr', 'start', 'end') before the first signal, we add +3
        # for rows: (l:nrow(df)): shifts the signal l times forward: instead of index 1 to N, pick l to N

        # Now, can you explain what the following line does?
        #y <- c(y, df[1:(nrow(df)-l+1), i+3])
      #}

      # calculating the correlation of lag l
      # this part is also a little tricky:
      #df_temp$autoCorrelation[?] <- cor(?, ?, method = c("pearson"))
      # so, we would have a df such as:
      # lag                 feature         autoCorrelation
      # 1*resolution        signal_name     cor_value_sigi_lag1
      # 2*resolution        signal_name     cor_value_sigi_lag1
      # 3*resolution        signal_name     cor_value_sigi_lag1
      # 4*resolution        signal_name     cor_value_sigi_lag1
      # .
      # .
      # .
    #}

    # binding df_temps together after the end of each iteration or each l (lag)
    #df_autoCor <- rbind(df_autoCor, df_temp)
  #}

  # after this loop df_autoCor would be:
  # lag                 feature         autoCorrelation
  # 1*resolution        signal_1        cor_value_sig1_lag1
  # 2*resolution        signal_1        cor_value_sig1_lag3
  # 3*resolution        signal_1        cor_value_sig1_lag4
  # .
  # .
  # .
  # 1*resolution        signal_2        cor_value_sig2_lag1
  # 2*resolution        signal_2        cor_value_sig2_lag2
  # 3*resolution        signal_2        cor_value_sig2_lag3
  # .
  # .
  # .

  #cat('Generating autocorrelation plot...\n')

  # SECTION: figure generation:
  #p <- ggplot(?, aes(?, ?, group = ?, color= ?)) +
  #  geom_line() +
  #  theme(legend.title=element_blank()) +
  #  labs(fill = "", x = "distance (bp)", y ="autocorrelation") +
  #  theme_classic(base_size = fontSize) +
  #  theme(legend.title=element_blank())

  # If x_label is not NA, add a label to x axis for plot p

  # If y_label is not NA, add a label to y axis for plot p



  # saving the generated figure
  #ggsave(paste0('./', outdir, '/', 'autoCorrelation.png'),
  #       plot = ?,
  #       width = ?, height = ?, units = "cm",
  #       dpi = 300, limitsize = TRUE)
}

