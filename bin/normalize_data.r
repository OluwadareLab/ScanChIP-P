#######################################################################
# This code was taken and manipulated from HPTAD's GitHub repository  #
# HPTAD: A computational method to identify topologically associating #
#   domains from HiChIP and PLAC-seq datasets (Rosen, et al. 2023)    #
#######################################################################

library(VGAM)
library(MASS)
library(argparse)
library(data.table)
library(dplyr)
library(zoo)
library(HiCcompare)



main <- function() {
  parser <- setup_parser()
  args <- parse_arguments(parser) # uncomment this line if you intend to run from the command line (and comment the next line)
  for (chrom in args$chroms) {
    long_filename <- paste0(args$indir, '/', args$prefix, '.', chrom, '.long.intra.bedpe.gz')
    bedpe <- fread(long_filename, header = F)
    bedpe <- bin_bedpe(bedpe, args$binsize) # assign each readpair to a binpair - no counting here
    bedpe <- filter_by_distance(bedpe, args$binsize, args$upperlimit)
    bedpe <- count_pairs(bedpe)
    bedpe <- filter_by_region(bedpe, args$filter, args$binsize)
    bedpe <- annotate_chip_sets(bedpe, args$chip, args$binsize)
    short_filename <- paste0(args$indir, '/', args$prefix, '.', chrom, '.shrt.vip.bed.gz')
    shorts <- fread(short_filename, header = FALSE)

    bedpe <- merge_short_counts(bedpe, shorts, args$binsize)
    features <- fread(args$features, header = FALSE)
    bedpe <- bedpe %>% group_by(pair_type) %>% create_features(features, args$binsize, args$upperlimit)
    bedpe <- filter_by_feature(bedpe, condition_trues = c('effective_length_left>0', 'effective_length_right>0'))
    bedpe <- filter_by_feature(bedpe, condition_trues = c('mappability_left>=0.8', 'mappability_right>=0.8'))
    bedpe <- bedpe %>% group_by(pair_type) %>% run_regression()
    normalized_bedpe <- bedpe[,c('x1','y1','normalized')]
    write.table(normalized_bedpe,"output.txt",sep="\t",row.names=FALSE)
  }
}


##################################################
###              HPTAD Functions               ###
##################################################


# Bin long reads based on defined bin size
bin_bedpe <- function(bedpe, binsize) {
  bedpe$bin1 <- ((bedpe$V2 + bedpe$V3) / 2) %/% binsize * binsize
  bedpe$bin2 <- ((bedpe$V5 + bedpe$V6) / 2) %/% binsize * binsize
  bedpe <- bedpe[, c('V1', 'bin1', 'V4', 'bin2')]
  bedpe$x1 <- ifelse(bedpe$bin1 <= bedpe$bin2, bedpe$bin1, bedpe$bin2)
  bedpe$y1 <- ifelse(bedpe$bin1 > bedpe$bin2, bedpe$bin1, bedpe$bin2)
  bedpe$x2 <- bedpe$x1 + binsize
  bedpe$y2 <- bedpe$y1 + binsize
  bedpe <- bedpe[, c('V1', 'x1', 'x2', 'V4', 'y1', 'y2')]
  colnames(bedpe) <- c('chr1', 'x1', 'x2', 'chr2', 'y1', 'y2')
  return (bedpe)
}


# Filter files by distance (bin size to upper bound)
filter_by_distance <- function(bedpe, lower, upper) {
  bedpe <- bedpe %>% filter((y1 - x1 >= lower) & (y1 - x1 <= upper))
  return (bedpe)
}


# Compute bin counts
count_pairs <- function(bedpe) {
  bedpe <- bedpe %>% dplyr::add_count(chr1,x1,chr2,y1, name = "count")
  bedpe <- unique(bedpe)
}


# Filter out exclusion regions
filter_by_region <- function(bedpe, filter_file, binsize) {
  if (filter_file == 'None') {
    return(bedpe)
  } else {
    filter_bed <- fread(filter_file, header = FALSE)
    filter_bed <- bin_bed(filter_bed, binsize)
    filter_bed <- filter_bed[,c('chr','start')]
    filter_bed$overlap <- 1
    bedpe <- merge(bedpe, filter_bed, by.x = c('chr1', 'x1'), by.y = c('chr', 'start'), all.x = TRUE)
    bedpe <- merge(bedpe, filter_bed, by.x = c('chr2', 'y1'), by.y = c('chr', 'start'), all.x = TRUE, suffixes = c('_left', '_right'))
    bedpe <- bedpe[(is.na(bedpe$overlap_left)) & (is.na(bedpe$overlap_right)),]
    bedpe <- bedpe[, c('chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'count')]
    return (bedpe)
  }
}


# Annotate bin pairs as AND, XOR, or NOT
annotate_chip_sets <- function(bedpe, chip_file, binsize) {
  chips <- fread(chip_file, header = FALSE)
  chips <- bin_bed(chips, binsize)
  chips$chip <- 1
  bedpe <- merge(bedpe, chips, by.x = c('chr1', 'x1'), by.y = c('chr', 'start'), all.x = TRUE)
  bedpe <- merge(bedpe, chips, by.x = c('chr2', 'y1'), by.y = c('chr', 'start'), all.x = TRUE, suffixes = c('_left', '_right'))
  bedpe$pair_type <- ifelse(is.na(bedpe$chip_left) & is.na(bedpe$chip_right), "NOT",
                            ifelse(!is.na(bedpe$chip_left) & !is.na(bedpe$chip_right), 'AND', 'XOR'))
  bedpe <- bedpe[, c('chr1','x1','x2','chr2','y1','y2','pair_type','count')]
  return (bedpe)
}


# Binning function for annotation
bin_bed <- function(d, binsize) {
  colnames(d)[1:3] <- c('chr', 'start', 'end')
  d$start <- d$start %/% binsize * binsize
  d$end <- d$end %/% binsize * binsize
  mask <- which(d$start == d$end)
  d_binsize <- d[mask,]
  d_binsize_exceeding <- d[-mask,]
  if (nrow(d_binsize_exceeding) > 0) {
    d_binsize_exceeding <- rbindlist(apply(d_binsize_exceeding, 1, function(entry, binsize) {
      start_points <- seq(as.integer(entry[2]), as.integer(entry[3]) - 1, binsize)
      expanded_d <- data.frame(chr = rep(entry[1], length(start_points)), start = start_points)
      return(expanded_d)
    }, binsize))
    d_binsize <- as.data.table(d_binsize)
    d_binsize_exceeding <- as.data.table(d_binsize_exceeding)
    d <- unique(rbind(dplyr::select(d_binsize, c("chr", "start")), dplyr::select(d_binsize_exceeding, c("chr", "start"))))
  }
  d$end <- d$start + binsize
  return(d)
}


# Merge the short counts
merge_short_counts <- function(bedpe, shorts, binsize) {
  shorts <- count_shorts(shorts, binsize)
  bedpe <- merge(bedpe, shorts, by.x = c('chr1', 'x1'), by.y = c('chr', 'start'), all.x = TRUE)
  bedpe <- merge(bedpe, shorts, by.x = c('chr2', 'y1'), by.y = c('chr', 'start'), all.x = TRUE, suffixes = c('_left', '_right'))
  bedpe <- bedpe[!is.na(bedpe$short_count_left) & !is.na(bedpe$short_count_right),]
  return (bedpe)
}


# Sum short counts for merge
count_shorts <- function(shorts, binsize) {
  shorts$start <- ((shorts$V2 + shorts$V3) / 2) %/% binsize * binsize
  shorts <- shorts[,c('V1', 'start')]
  colnames(shorts)[1] <- 'chr'
  shorts <- shorts %>% dplyr::add_count(chr, start, name = "short_count")
  shorts <- unique(shorts)
  return (shorts)
}


# Create features from supplied feature file
create_features <- function(bedpe, features, binsize, binrange) {
  colnames(features) <- c("chr", "start", "end", "effective_length", "gc", "mappability")
  features$end <- NULL
  bedpe <- merge(bedpe, features, by.x = c("chr1", "x1"), by.y = c("chr", "start"), all.x = TRUE)
  bedpe <- merge(bedpe, features, by.x = c("chr2", "y1"), by.y = c("chr", "start"), all.x = TRUE, suffixes = c('_left', '_right'))
  dist <- abs(bedpe$x1 - bedpe$y1)
  bedpe$log_dist <- log((1 + dist) / binrange)
  bedpe$log_len <- log((bedpe$effective_length_left + 1) * (bedpe$effective_length_right + 1) / (binsize ** 2))
  bedpe$log_gc <- log(bedpe$gc_left * bedpe$gc_right)
  bedpe$log_map <- log(bedpe$mappability_left * bedpe$mappability_right)
  max_short_left <- max(bedpe$short_count_left)
  max_short_right <- max(bedpe$short_count_right)
  bedpe$log_short <- log(((bedpe$short_count_left + 1) * (bedpe$short_count_right + 1)) /
                           ((max_short_left + 1) * (max_short_right + 1)))
  bedpe <- bedpe[,c('chr1','x1','x2','chr2','y1','y2','pair_type','count',
                    'effective_length_left','effective_length_right',
                    'gc_left','gc_right','mappability_left','mappability_right',
                    'short_count_left','short_count_right',
                    'log_dist','log_len','log_gc','log_map','log_short')]
  return(bedpe)
}


# Filter based on feature values
filter_by_feature <- function(data, condition_trues) {
  for (condition in condition_trues) {
    data <- data[eval(parse(text=paste0('data$', condition))),]
  }
  return (data)
}


# Run positive poisson regression
run_regression <- function(data) {
  model <- vglm(count ~ 1 + log_len + log_gc + log_map + log_short + log_dist, family = 'pospoisson', data = data)
  data$expected <- model@fitted.values
  data$normalized <- data$count / data$expected

  data$pvalue <- ppois(data$count, data$expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, data$expected, lower.tail = FALSE, log.p = FALSE)
  data$fdr <- p.adjust(data$pvalue, method = 'fdr')
  return (data)
}

setup_parser <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument('-i', '--indir', help = 'input directory containing bedpe files', required = TRUE)
  parser$add_argument('-o', '--outdir', help = 'output directory', required = TRUE)
  parser$add_argument('-p', '--prefix', help = 'name of the dataset; should be the same prefix that input bedpe files start with')
  parser$add_argument('-C', '--chip', help = 'filepath to the ChIP peaks', required = TRUE) 
  parser$add_argument('-c', '--chromosomes', help = '"A comma-separated list of chromosomes', required = TRUE)
  parser$add_argument('-f', '--features', help = 'path of the file containing genomic features (mappability, gc content, effective length)', required = TRUE)
  parser$add_argument('-x', '--filter', help = '"None" or path to the bed file containing regions to be filtered', default = 'None', required = FALSE)
  parser$add_argument('-b', '--binsize', help = 'bin size, default = 40000', default = 40000, type = 'integer')
  parser$add_argument('-u', '--upperlimit', help = 'upper limit for distance between bins, default = 2000000', default = 2000000, type = 'integer')
  return (parser)
}

parse_arguments <- function(parser) {
  args <- parser$parse_args()
  message(paste('Reading input from directory:', args$indir))
  message(paste('Output will be written to:', args$outdir))
  message(paste('Genomic features file specified:', args$features))
  chromosomes <- strsplit(args$chromosomes, ',')[[1]]
  chromosomes <- paste0("chr", chromosomes)
  args$chroms <- chromosomes
  message('Analyzing')
  message(paste('   ', args$chroms))
  if(args$filter == 'None') {
    message('No filter region specified')
  } else {
    message('Filtered bins will be read from:', args$filter)
  }
  message(paste('Binsize:', args$binsize))
  message(paste('Prefix:', args$prefix))
  message(paste('Upperlimit:', args$upperlimit))
  return (args)
}

main()
