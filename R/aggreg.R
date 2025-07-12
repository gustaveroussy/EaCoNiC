## PACKAGES USED
## CRAN : tibble, dplyr, cairoDevice
## BIOCONDUCTOR : GenomeInfoDb, GenomicRanges, rtracklayer, ComplexHeatmap
## GITHUB : jchiquet/blockseg



# SEG.rds.list <- list.files(path = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/CNA', pattern = '*.SEG.*.RDS$', full.names = TRUE, recursive = TRUE)

## User function to set the cache directory
# EaCoNiC.set.cache <- function(cache.dir = NULL) {
#   if(is.null(cache.dir)) {
#     message('No directory specified :')
#     if ('EaCoNiC.cache' %in% utils::installed.packages()) {
#       message("\tAn installation of 'EaCoNiC.cache' was found and will be used.")
#       cache.dir <- paste0(base::tempdir(check = TRUE), '/EaCoNiC.cache/')
#     } else {
#       cache_dir <- tools::R_user_dir("recount3", "cache")
#       message("\tA cache folder will be used.")
#       BiocFileCache::BiocFileCache(cache_dir)
#     }
#   }
#   if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
#   message(paste0('EaCoNiC cache directory set to [', cache.dir, ']'))
#   # options(EACONIC_CACHE=cache.dir)
# }

## Get annotations from the UCSC and cache them locally (works for multiple genomes)
## USING https://github.com/jeffbhasin/goldmine
# EaCoNiC.cache.annotation <- function(genome = NULL, cache.dir = NULL, cache.update = FALSE) {
#   
#   ## Getting the list of available registered genomes @ UCSC (from GenomeInfoDb)
#   avail.genz <- sort(GenomeInfoDb::registered_UCSC_genomes()$genome)
#   
#   if (is.null(genome)) stop("A UCSC genome name is required. It should be one of : '", paste(avail.genz, collapse = "', '"), "'")
#   if(is.null(cache.dir)) stop('A cache directory to write files is required !')
#   if(!file.exists(cache.dir)) stop(paste0('Cache directory [', cache.dir, '] does not exist !'))
#   
#   message('An internet connexion is required when cached data are not available.')
#   
#   genome <- genome[genome %in% avail.genz]
#   if(length(genome) == 0) stop('None of the specified genome(s) is available @ UCSC !')
#   
#   annot.type <- c('cytoBandIdeo', 'cpgIslandExt', 'genomicSuperDups', 'dgvMerged')
#   
#   annot.obj <- list()
#   # annot.obj <- setNames(vector(mode = 'list', length = length(annot.type)), annot.type)
#   message(paste0('Getting annotations for ', genome, ' ...'))
#   
#   ## Trying to get chromosomes info
#   trycI <- try(cI <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome, assembled.molecules.only = TRUE, as.Seqinfo = TRUE), silent = TRUE)
#   ## Getting RefSeq genes
#   annot.obj$refGene <- goldmine::getGenes(geneset = 'refseq', genome = genome, cachedir = cache.dir)
#   ## Cycling in annotation tables
#   for (tt in annot.type) {
#     tryAn <- try(annot.obj[[tt]] <- goldmine::getUCSCTable(table = tt, genome = genome, cachedir = cache.dir, sync = cache.update), silent = TRUE)
#     if (!class(tryAn)[1] == 'try-error') annot.obj[[tt]][['bin']] <- NULL
#   }
#   ## Convert to GRanges
#   for (tt in seq_along(annot.obj)) {
#     annot.obj[[tt]] <- GenomicRanges::makeGRangesFromDataFrame(df = annot.obj[[tt]], keep.extra.columns = TRUE, ignore.strand = FALSE, seqnames.field = colnames(annot.obj[[tt]])[1], start.field = colnames(annot.obj[[tt]])[2], end.field = colnames(annot.obj[[tt]])[3], starts.in.df.are.0based = TRUE)
#     ## Adding seqinfo when available
#     if (!class(trycI)[1] == 'try-error') {
#       GenomeInfoDb::seqlevels(annot.obj[[tt]], pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(cI)
#       GenomeInfoDb::seqinfo(annot.obj[[tt]]) <- cI
#     }
#   }
#   
#   ## Merging refGenes
#   `%>%` <- dplyr::`%>%`
#   toto <- suppressMessages(annot.obj$refGene %>% tibble::as_tibble() %>% dplyr::group_by(seqnames, strand, name) %>% dplyr::summarize(start = min(start), end = max(end)) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
# 
#   ## Saving is there is something to save
#   # if(!all(sapply(X = annot.obj, is.null))) saveRDS(object = annot.obj, file = paste0(cache.dir, '/', genome, '.RDS'), compress = 'bzip2')
#   if(!all(sapply(X = annot.obj, is.null))) save(object = annot.obj, file = paste0(cache.dir, '/', genome, '.rda'), compress = 'bzip2')
#   
#   # return(annot.obj)
# }


## Get annotations from the UCSC
## USING https://github.com/jeffbhasin/goldmine
## Default (forced) annot.type are : c('cytoBandIdeo', 'cpgIslandExt', 'genomicSuperDups', 'dgvMerged'). Any UCSC table is accepted, as long as it does have the [chr], [start], [end] data in columns 1:3, or 2:4 if first column is [bin].
get.annotation <- function(genome = NULL, annot.type = c('cytoBandIdeo', 'cpgIslandExt', 'genomicSuperDups', 'dgvMerged'), cache.dir = NULL, cache.update = FALSE, return.data = TRUE) {
  
  ## Getting the list of available registered genomes @ UCSC (from GenomeInfoDb)
  # avail.genz <- sort(GenomeInfoDb::registered_UCSC_genomes()$genome)
  
  ## Getting the list of available registered genomes @ UCSC (from rtracklayer)
  mySession <- rtracklayer::browserSession()
  avail.genz <- rtracklayer::ucscGenomes()$db
  
  
  if (is.null(genome)) stop("A UCSC genome name is required. It should be one of : '", paste(avail.genz, collapse = "', '"), "'")
  
  if(is.null(cache.dir)) {
    message("No cache directory provided : annotation data won't be kept after usage ...")
  } else if(!file.exists(cache.dir)) {
    message(paste0('Provided cache directory [', cache.dir, '] does not exist. It will be created now as new cache. An internet connexion is required.'))
  } else if (cache.update) {
    message(paste0('Provided cache directory [', cache.dir, '] does exist. It will be USED and UPDATED. An internet connexion is required.'))
  } else {
    message(paste0('Provided cache directory [', cache.dir, '] does exist. It will be USED.'))
  }
  
  ## Check genome
  genome <- genome[genome %in% avail.genz]
  if(length(genome) == 0) stop('None of the specified genome(s) is available @ UCSC !')
  rtracklayer::genome(mySession) <- genome
  
  ## Get genome-cache dir
  gcd <- paste0(cache.dir, '/', genome)
  
  ## Create cache dir if needed
  if (!dir.exists(gcd)) dir.create(path = gcd, recursive = TRUE)
  
  ## Structure annot.type
  annot.type <- unique(c('ncbiRefSeq', annot.type))
  
  ## Prepare output
  annot.obj <- list()
  
  ## Trying to get chromosomes info
  trycI <- try(cI <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome, assembled.molecules.only = TRUE, as.Seqinfo = TRUE), silent = TRUE)
  
  ## Cycling in annotation tables
  message(paste0('Getting annotations for ', genome, ' ...'))
  for (tt in annot.type) {
    tt_file <- paste0(gcd, '/', tt, '.RDS')
    if (!file.exists(tt_file) | cache.update) {
      message('. Retrieving [', tt, '] from UCSC ...')
      tryAn <- try(annot.obj[[tt]] <- rtracklayer::getTable(rtracklayer::ucscTableQuery(x = mySession, table = tt)), silent = TRUE)
      if (!class(tryAn)[1] == 'try-error') {
        ### Remove the bin column (if available)
        annot.obj[[tt]][['bin']] <- NULL
        ### Handle the ncbiRefSeq special case
        if (tt == 'ncbiRefSeq') {
          annot.obj[[tt]] <- annot.obj[[tt]][,c('chrom', 'txStart', 'txEnd', 'strand', 'name', 'name2')]
          ### Merging refGenes
          `%>%` <- dplyr::`%>%`
          annot.obj[[tt]] <- suppressMessages(annot.obj[[tt]] %>% tibble::as_tibble() %>% dplyr::group_by(chrom, strand, name2) %>% dplyr::summarize(start = min(txStart), end = max(txEnd)) %>% as.data.frame())
          annot.obj[[tt]] <- annot.obj[[tt]][,c('chrom', 'start', 'end', 'strand', 'name2')]
          colnames(annot.obj[[tt]])[colnames(annot.obj[[tt]]) == 'name2'] <- 'name'
        }
        ### Convert to GRange
        annot.obj[[tt]] <- GenomicRanges::makeGRangesFromDataFrame(df = annot.obj[[tt]], keep.extra.columns = TRUE, ignore.strand = FALSE, seqnames.field = colnames(annot.obj[[tt]])[1], start.field = colnames(annot.obj[[tt]])[2], end.field = colnames(annot.obj[[tt]])[3], starts.in.df.are.0based = TRUE)
        ### Add chr data when available
        if (!class(trycI)[1] == 'try-error') {
          GenomeInfoDb::seqlevels(x = annot.obj[[tt]], pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(cI)
          GenomeInfoDb::seqinfo(annot.obj[[tt]]) <- cI
        }
        ### Save in cache
        message('.. Saving in cache [', tt_file, '] ...')
        saveRDS(object = annot.obj[[tt]], file = tt_file, compress = 'xz')
      }
    } else {
      message('Loading [', tt, '] from cache ...')
      annot.obj[[tt]] <- readRDS(tt_file)
    }
  }
  
  ## Return the annotation object
  if (return.data) return(annot.obj)
}


## Get annotations from the UCSC and cache them locally (works for multiple genomes)
## USING rtracklayer
# EaCoNiC.cache.annotation_old <- function(genome = NULL, cache.dir = NULL) {
#   
#   ## Getting the list of available registered genomes @ UCSC (from GenomeInfoDb)
#   avail.genz <- sort(GenomeInfoDb::registered_UCSC_genomes()$genome)
#   
#   if (is.null(genome)) stop("A UCSC genome name is required. It should be one of : '", paste(avail.genz, collapse = "', '"), "'")
#   if(is.null(cache.dir)) stop('A cache directory to write files is required !')
#   if(!file.exists(cache.dir)) stop(paste0('Cache directory [', cache.dir, '] does not exist !'))
#   
#   message('An internet connexion is required.')
# 
#   ## Setting the cache dir
#   # EaCoNiC.set.cache(cache.dir = cache.dir)
# 
#   genome <- genome[genome %in% avail.genz]
#   if(length(genome) == 0) stop('None of the specified genome(s) is available @ UCSC !')
#   
#   annot.type <- list(
#     cytoBandIdeo = c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'),
#     refGene = c('chrom', 'txStart', 'txEnd', 'strand', 'exonCount', 'name2'),
#     cpgIslandExt = c('chrom', 'chromStart', 'chromEnd', 'name', 'perCpG', 'perGc', 'obsExp'),
#     genomicSuperDups = c('chrom', 'chromStart', 'chromEnd', 'strand'),
#     dgvMerged = c('chrom', 'chromStart', 'chromEnd', 'strand', 'varType'))
#   
#   for (gen in genome) {
#     annot.obj <- setNames(vector(mode = 'list', length = length(annot.type)), names(annot.type))
#     message(paste0('Getting annotations for ', gen, ' ...'))
#     ## Trying to get chromosomes info
#     trycI <- try(cI <- GenomeInfoDb::getChromInfoFromUCSC(genome = gen, assembled.molecules.only = TRUE, as.Seqinfo = TRUE), silent = TRUE)
#     ## Cycling in annotation tables
#     for (tt in names(annot.type)) {
#       ## Querying
#       tryres <- try(query <- rtracklayer::ucscTableQuery(gen, table = tt), silent = TRUE)
#       if (class(tryres) == 'try-error') {
#         message(paste0("\tGenome [", gen, '] did NOT have ', tt, '.'))
#       } else {
#         df <- rtracklayer::getTable(query)
#         message(paste0("\tRetrieved ", tt, '.'))
#         ## Limiting to desired columns
#         df <- df[, colnames(df) %in% annot.type[[tt]]]
#         ## Converting to GRanges
#         gr <- GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, ignore.strand = FALSE, seqnames.field = annot.type[[tt]][1], start.field = annot.type[[tt]][2], end.field = annot.type[[tt]][3], starts.in.df.are.0based = TRUE)
#         ## Adding seqinfo when available
#         if (!class(trycI) == 'try-error') {
#           GenomeInfoDb::seqlevels(gr, pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(cI)
#           GenomeInfoDb::seqinfo(gr) <- cI
#         }
#         ## Merging UCSC refGenes
#         if (tt == 'refGene') {
#           `%>%` <- dplyr::`%>%`
#           gr <- suppressMessages(gr %>% tibble::as_tibble() %>% dplyr::group_by(seqnames, strand, name2) %>% dplyr::summarize(start = min(start), end = max(end)) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
#         }
#         ## Saving in the object
#         annot.obj[[tt]] <- gr
#       }
#     }
#     ## Saving is there is something to save
#     if(!all(sapply(X = annot.obj, is.null))) save(object = annot.obj, file = paste0(cache.dir, '/', gen, '.rda'), compress = 'bzip2')
#   }
# }

## [NOT EXPORTED] Get cached genome
# get.cached.genome <- function(genome = 'hg19') {
#   ### Getting cache directory if available in options()
#   cache.dir <- getOption('EACONIC_CACHE')
#   ### Checking if a cache file exists
#   cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.rda'))
#   if (is.null(cache.dir) | !cache.exists) {
#     ## If there is no cache directory in options()
#     message('No cache directory found, or no cache found for provided genome : Annotation will be retrieved to a temporary directory (an internet connection is required).')
#     cache.dir <- paste0(tempdir(check = TRUE), '/EaCoNiC_cache/')
#     dir.create(cache.dir)
#     EaCoNiC.cache.annotation(cache.dir = cache.dir, genome = genome)
#   }
#   ## Loading cache
#   message(paste0('Loading cache for genome ', genome, ' ...'))
#   load(paste0(cache.dir, '/', genome, '.rda'), envir = .GlobalEnv)
#   # annot.obj <- readRDS(paste0(cache.dir, '/', genome, '.RDS'), envir = .GlobalEnv)
#   # return(annot.obj)
# }

## [NOT EXPORTED] Convert a cached dataset of chromosome data from a specified genome to a more complex object
genome2cs <- function(genome = NULL) {
  ## Checks
  ### Getting cache directory if available in options()
  cache.dir <- getOption('EACONIC_CACHE')
  ### Checking if a cache file exists
  cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.rda'))
  if (is.null(cache.dir) | !cache.exists) {
    ## If there is no cache directory in options()
    message('No cache directory found, or no cache found for provided genome : Annotation will be retrieved to a temporary directory (an internet connection is required).')
    cache.dir <- paste0(tempdir(check = TRUE), '/EaCoNiC_cache/')
    dir.create(cache.dir)
    EaCoNiC.cache.annotation(cache.dir = cache.dir, genome = genome)
  }
  ## Loading cache
  message('Retrieving cytobands from cache ...')
  load(paste0(cache.dir, '/', genome, '.rda'))
  
  if (is.null(annot.obj$cytoBandIdeo)) {
    message("\tRequested genome did not have cytobands data !")
    return(NULL)
  }
}

## [NOT EXPORTED] Get cytobands (may require an internet connection if genome is not cached). Filtering of uncertain cytobands can be performed optionally.
# get.cytobands <- function(genome = NULL, uncertain.filter = TRUE) {
#   ## Checks
#   ### Getting cache directory if available in options()
#   cache.dir <- getOption('EACONIC_CACHE')
#   ### Checking if a cache file exists
#   # cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.RDS'))
#   cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.rda'))
#   if (is.null(cache.dir) | !cache.exists) {
#     ## If there is no cache directory in options()
#     message('No cache directory found, or no cache found for provided genome : Annotation will be retrieved to a temporary directory (an internet connection is required).')
#     cache.dir <- paste0(tempdir(check = TRUE), '/EaCoNiC_cache/')
#     dir.create(cache.dir)
#     EaCoNiC.cache.annotation(cache.dir = cache.dir, genome = genome)
#   }
#   ## Loading cache
#   message('Retrieving cytobands from cache ...')
#   # annot.obj <- readRDS(file = paste0(cache.dir, '/', genome, '.RDS'))
#   load(file = paste0(cache.dir, '/', genome, '.rda'), envir = .GlobalEnv)
#   if (is.null(annot.obj$cytoBandIdeo)) {
#     message("\tRequested genome did not have cytobands data !")
#     return(NULL)
#   }
#   gr <- annot.obj$cytoBandIdeo
#   # rm(annot.obj)
#   if(uncertain.filter) {
#     if (!'gieStain' %in% colnames(GenomicRanges::mcols(gr))) {
#       message("\tRequested genome did not have staining, so filering of uncertain cytobands was not possible.")
#     } else {
#       bandcolz <- biovizBase::getBioColor(type = 'CYTOBAND')
#       unc.bands <- gr$gieStain %in% names(bandcolz)[c(2,3,5)]
#       message(paste0("\tFiltering ", length(which(unc.bands)), ' uncertain bands ...'))
#       gr <- gr[!unc.bands,]
#     }
#   }
#   return(gr)
# }

## Simplify a vector of chr names by identifying and deleting its common root 
chr.simplifier <- function(chr.names = NULL) {
  message('Simplfying chromosome names ...')
  cs.decomp <- lapply(chr.names, function(x) { unlist(strsplit(x, split = '')) })
  cs.len <- sapply(cs.decomp, length)
  cs.com <- lapply(seq_len(min(cs.len)), function(x) {
    unique(vapply(cs.decomp, function(y) { y[x] }, 'a'))
  })
  com.key <- paste(vapply(cs.com, function(x) { if(length(x) == 1) x else ''}, 'a'), collapse = '')
  if (com.key == '') return(chr.names) else return(gsub(pattern = com.key, replacement = '', chr.names))
}

## Aggregate samples to form a cohort 
EaCoNiC.SEG.aggregate <- function(SEG.rds.list = list.files(path = getwd(), pattern = '*.SEG.*.RDS*', full.names = TRUE, recursive = TRUE), L2R = TRUE, BAF = TRUE, cohort.name = 'TEST', out_dir = getwd(), return.data = FALSE) {
  ## Early checks
  if (is.null(SEG.rds.list)) stop("A list of *SEG*.RDS files (results from EaCoN::Segment()) is required !")
  if (length(SEG.rds.list) < 2) stop("To build a cohort, at least two samples are required !")
  if (!any(c(L2R, BAF))) stop("At least one of 'L2R' and 'BAF' should be set to TRUE.")
  if (is.null(out_dir) & !return.data) stop("No output directory specified, and 'return.data' is set to FALSE : as nothing will be output, nothing will be performed !")
  if (is.null(cohort.name)) stop("A cohort name is required !")
  
  ## Loading data : First pass to control data compatibility
  message('Reading samples metadata to ensure compatibility ...')
  samp.meta.list <- lapply(seq_along(SEG.rds.list), function(x) {
    message(paste0("\tParsing metadata from file ", x, ' / ', length(SEG.rds.list), ' ...'))
    sobj <- readRDS(SEG.rds.list[x])
    return(list(sample.name = sobj$meta$basic$samplename, loss.cut = sobj$meta$eacon[['L2R-segments-loss-cutoff']], gain.cut = sobj$meta$eacon[['L2R-segments-gain-cutoff']], genome =  sobj$meta$basic$genome, segmenter = sobj$meta$eacon$segmenter))
  })
  samp.info <- as.data.frame(t(as.data.frame(lapply(samp.meta.list, unlist), check.names = FALSE)))
  samp.info$loss.cut <- as.numeric(samp.info$loss.cut)
  samp.info$gain.cut <- as.numeric(samp.info$gain.cut)
  rownames(samp.info) <- NULL
  
  ## Integrity checks
  ### Looking for replicated sample names
  message('Checking samples compatibility ...')
  if (any(duplicated(samp.info$sample.name))) stop('Some sample names are replicated ! Unique sample names are required.')
  my.genome <- unique(samp.info$genome)
  if (length(my.genome) > 1) stop('Some samples correspond to different genome builds ! Unique genome build is required.')
  if (length(table(samp.info$segmenter)) > 1) warning('Some samples correspond to different segmenters. While this does not prevent the aggregation into a cohort, this is not recommended.')
  
  ### TEST based on rtracklayer
  chr.gr <- get.cytobands(genome = my.genome, uncertain.filter = TRUE)
  
  ## Loading L2R segmentation data
  reg.l2r.df.filt.cut <- reg.l2r.df.filt <- NULL
  if (L2R) {
    message('Loading L2R segmentation results ...')
    l2r.gr.list <- lapply(X = seq_along(SEG.rds.list), function(x) {
      message(paste0("\tLoading segments from file ", x, ' / ', length(SEG.rds.list), ' ...'))
      sobj <- readRDS(SEG.rds.list[x])
      samp.l2r.gr <- GenomicRanges::makeGRangesFromDataFrame(df = data.frame(seqnames = sobj$data$chrs[sobj$cbs$nocut$Chr], start = sobj$cbs$nocut$Start, end = sobj$cbs$nocut$End, L2R = sobj$cbs$nocut$Log2Ratio), keep.extra.columns = TRUE, ignore.strand = TRUE)
      
      ## Adding seqinfo when possible
      if(!is.null(chr.gr)) {
        ref.seqinfo <- GenomeInfoDb::seqinfo(chr.gr)
        if (!is.null(ref.seqinfo)) {
          GenomeInfoDb::seqlevels(samp.l2r.gr, pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(chr.gr)
          GenomeInfoDb::seqinfo(samp.l2r.gr) <- GenomeInfoDb::seqinfo(chr.gr)
        }
      }
      return(samp.l2r.gr)
    })
    names(l2r.gr.list) <- samp.info$sample.name
    
    ## Getting all possible intervals from the cohort merge.
    if(is.null(chr.gr)) {
      samps.l2r.gr <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(l2r.gr.list)))
    } else {
      unl.l2r.gr <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(c(l2r.gr.list, chr = IRanges::reduce(chr.gr)))))
      samp.chr.ov <- GenomicRanges::findOverlaps(query = unl.l2r.gr, subject = IRanges::reduce(chr.gr), minoverlap = 1L)
      samps.l2r.gr <- unl.l2r.gr[S4Vectors::from(samp.chr.ov)]
    }
    
    ## Synching to a DF
    message('Merging all breakpoints ...')
    reg.l2r.df <- cbind(as.data.frame(samps.l2r.gr)[,1:3], as.data.frame(lapply(seq_along(l2r.gr.list), function(x) {
      message(paste0("\tAdding ", samp.info$sample.name[x], '...'))
      GenomicRanges::mcols(samps.l2r.gr) <- S4Vectors::DataFrame(L2R = NA)
      over.gr <- IRanges::findOverlaps(query = samps.l2r.gr, subject = l2r.gr.list[[x]], minoverlap = 1L)
      samps.l2r.gr$L2R[over.gr@from] <- l2r.gr.list[[x]]$L2R[over.gr@to]
      return(samps.l2r.gr$L2R)
    }), col.names = samp.info$sample.name, check.names = FALSE, stringsAsFactors = FALSE))
    colnames(reg.l2r.df)[1:3] <- c('Chr', 'Start', 'End')
    
    ## Filter intervals with only NAs
    message('Filtering empty intervals ...')
    reg.l2r.df.filt <- reg.l2r.df[!is.na(matrixStats::rowAnys(as.matrix(reg.l2r.df[,-c(1:3)]))),]
    
    ## Applying L2R cut
    message('Setting normal segments ...')
    reg.l2r.df.filt.cut <- reg.l2r.df.filt
    for (x in 4:ncol(reg.l2r.df.filt.cut)) {
      reg.l2r.df.filt.cut[,x][reg.l2r.df.filt.cut[,x] > as.numeric(samp.info$loss.cut[(x-3)]) & reg.l2r.df.filt.cut[,x] < as.numeric(samp.info$gain.cut[(x-3)])] <- 0
    }
  }
  
  ### BAF )
  reg.baf.df.filt <- NULL
  if (BAF) {
    ## Converting L2R (from cbs NOcut) to GRanges, combining them with chr.gr to a GRangesList
    message('Loading BAF segmentation results ...')
    baf.gr.list <- lapply(X = seq_along(SEG.rds.list), function(x) {
      message(paste0("\tLoading segments from file ", x, ' / ', length(SEG.rds.list), ' ...'))
      sobj <- readRDS(SEG.rds.list[x])
      baf.res <- data.frame(sobj$data$SNPpos[as.numeric(rownames(sobj$data$Tumor_BAF_segmented[[1]])),], baf = sobj$data$Tumor_BAF_segmented[[1]][,1])
      baf.tbl <- tibble::as_tibble(baf.res) 
      baf.tbl <- dplyr::group_by(baf.tbl, chrs, bafrle = with(rle(baf), rep(seq_along(lengths), lengths)))
      baf.df <- suppressMessages(as.data.frame(dplyr::summarise(baf.tbl, segStart = min(pos), segEnd = max(pos), BAF = min(baf))))
      baf.df$bafrle <- NULL
      samp.baf.gr <- GenomicRanges::makeGRangesFromDataFrame(df = baf.df, seqnames.field = 'chrs', start.field = 'segStart', end.field = 'segEnd', keep.extra.columns = TRUE, ignore.strand = TRUE)
      
      ## Adding seqinfo when possible
      if(!is.null(chr.gr)) {
        ref.seqinfo <- GenomeInfoDb::seqinfo(chr.gr)
        if (!is.null(ref.seqinfo)) {
          GenomeInfoDb::seqlevels(samp.baf.gr, pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(chr.gr)
          GenomeInfoDb::seqinfo(samp.baf.gr) <- GenomeInfoDb::seqinfo(chr.gr)
        }
      }
      return(samp.baf.gr)
    })
    names(baf.gr.list) <- samp.info$sample.name
    
    if(is.null(chr.gr)) {
      samps.baf.gr <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(baf.gr.list)))
    } else {
      unl.baf.gr <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(c(baf.gr.list, chr = IRanges::reduce(chr.gr)))))
      samp.chr.ov <- GenomicRanges::findOverlaps(query = unl.baf.gr, subject = IRanges::reduce(chr.gr), minoverlap = 1L)
      samps.baf.gr <- unl.baf.gr[S4Vectors::from(samp.chr.ov)]
    }
    
    ## Getting all possible intervals from the cohort merge.
    # unl.baf.gr <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(c(baf.gr.list, chr = chr.gr))))
    
    ## Synching to a DF
    message('Merging all breakpoints ...')
    reg.baf.df <- cbind(as.data.frame(samps.baf.gr)[,1:3], as.data.frame(lapply(seq_along(baf.gr.list), function(x) {
      message(paste0("\tAdding ", samp.info$sample.name[x], '...'))
      GenomicRanges::mcols(samps.baf.gr) <- S4Vectors::DataFrame(BAF = NA)
      over.gr <- IRanges::findOverlaps(query = samps.baf.gr, subject = baf.gr.list[[x]], minoverlap = 1L)
      samps.baf.gr$BAF[over.gr@from] <- baf.gr.list[[x]]$BAF[over.gr@to]
      return(samps.baf.gr$BAF)
    }), col.names = samp.info$sample.name, check.names = FALSE, stringsAsFactors = FALSE))
    colnames(reg.baf.df)[1:3] <- c('Chr', 'Start', 'End')
    
    ## Filter intervals with only NAs
    reg.baf.df.filt <- reg.baf.df[!is.na(matrixStats::rowAnys(as.matrix(reg.baf.df[,-c(1:3)]))),]
  }
  out.obj <- list(samples.info = samp.info, L2R.cut = reg.l2r.df.filt.cut, L2R.nocut = reg.l2r.df.filt, BAF = reg.baf.df.filt, meta = list(type = 'SEG', cohort.name = cohort.name, genome = my.genome, sample.files = SEG.rds.list))
  if (!is.null(out_dir)) saveRDS(object = out.obj, file = paste0(out_dir, '/', cohort.name, '.SEG.COHORT.RDS'), compress = 'bzip2')
  if (return.data) return(out.obj) else return(invisible(NULL))
}

# EaCoN.cache.annotation(cache.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/EaCoNiC_cache', genome = 'hg19')

# EaCoN.set.cache(cache.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/EaCoN_cache')

# EaCoNiC.SEG.aggregate(SEG.rds.list = list.files(path = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/CNA', pattern = '*.SEG.*.RDS$', full.names = TRUE, recursive = TRUE), out_dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev', L2R = TRUE, BAF = TRUE, cohort.name = 'TEST2', return.data = FALSE)

## Annotate genomic regions
regions.annotate <- function(regions_df = NULL, seqnames = 'Chr', startname = 'Start', endname = 'End', annot.type = c('ncbiRefSeq', 'cytoBandIdeo', 'cpgIslandExt', 'genomicSuperDups', 'dgvMerged'), annot.match = c('in', 'out', 'in', 'in', 'in'), genome = 'hg38', cache.dir = NULL) {
  
  if (is.null(regions_df)) message("No coordinates data.frame provided : can't annotate !")
  
  if (is.null(cache.dir)) message('No cache dir provided : annotation data will be retrieved ...') else message('Annotation data will be retrieved from cache ...')
  
  if (length(annot.match) != length(annot.type)) stop("'annot.type' and 'annot.match' should have the same length !")
  # 
  if (!all(tolower(annot.match) %in% c('in', 'out'))) stop("'annot.match' should be one of 'in' or 'out'")
  
  ## Get annotation
  annot.obj <- get.annotation(genome = genome, annot.type = annot.type, cache.dir = cache.dir, cache.update = FALSE, return.data = TRUE)
  ## Check remaining annotation types
  match_sn <- setNames(object = annot.match, nm = annot.type)
  annot.type <- names(annot.obj)
  match_sn <- match_sn[annot.type]
  
  ## Add width
  regions_df$Width <- regions_df[[endname]] - regions_df[[startname]] +1
  
  ## Cycle annotations
  for (ann in annot.type) {
    message('. ', ann)
    if (unname(match_sn[ann] == 'out')) {
      ## Start
      starts_gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions_df[, !colnames(regions_df) %in% endname], keep.extra.columns = FALSE, ignore.strand = TRUE, seqnames.field = seqnames, start.field = startname, end.field = startname)
      starts_ol <- GenomicRanges::findOverlaps(query = starts_gr, subject = annot.obj[[ann]], minoverlap = 1L)
      regions_df[[paste0(ann, '.Start')]][S4Vectors::from(starts_ol)] <- annot.obj[[ann]]$name[S4Vectors::to(starts_ol)]
      rm(starts_gr, starts_ol)
      ## End
      ends_gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions_df[, !colnames(regions_df) %in% startname], keep.extra.columns = FALSE, ignore.strand = TRUE, seqnames.field = seqnames, start.field = endname, end.field = endname)
      ends_ol <- GenomicRanges::findOverlaps(query = ends_gr, subject = annot.obj[[ann]], minoverlap = 1L)
      regions_df[[paste0(ann, '.End')]][S4Vectors::from(ends_ol)] <- annot.obj[[ann]]$name[S4Vectors::to(ends_ol)]
      rm(ends_gr, ends_ol)
    } else {
      regions_gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions_df, keep.extra.columns = FALSE, ignore.strand = TRUE, seqnames.field = seqnames, start.field = startname, end.field = endname)
      
      regions_df[[paste0(ann, '.N')]] <- GenomicRanges::countOverlaps(query = regions_gr, subject = annot.obj[[ann]])
      
      in_oldf <- as.data.frame(GenomicRanges::findOverlaps(query = regions_gr, subject = annot.obj[[ann]], minoverlap = 1L))
      regions_df[[ann]] <- NA
      unique.hits <- unique(in_oldf$queryHits)
      in.list <- sapply(unique.hits, function(x) {
        paste(annot.obj[[ann]]$name[in_oldf$subjectHits[in_oldf$queryHits == x]], collapse = ',')
      }, simplify = TRUE)
      regions_df[[ann]][unique.hits] <- in.list
      
      rm(regions_gr, in_oldf, unique.hits, in.list)
    }
  }
  return(regions_df)
}


## Perform differential analysis on segmentation regions (... correspond to additional parameters to read.table, to parse the annotation table)
## SEG.cohort.obj     Robject                 The output from the EaCoNiC.SEG.aggregate() function
## factors.df         dataframe               A data frame that contains at least two columns : one with samplenames, the other(s) being factors to use to perform the diff tests
## samples.colname   character               Column name of factors.file to identify samples
## factors.colname       vector of characters    Column name(s) of factors.file to specify classes to compare
## test.2             character               Type of test to use for factors with 2 classes. 'W' for Wilcoxon, 'T' for Student's t-test.
## test.mult          character               Type of test to use for factors with >2 classes. 'ANOVA' is self-explanatory, 'KW' for Kruskal-Wallis.
## data.type          character               Type of data to use from the SEG.cohort.obj. 'L2R.cut' for log2ratio values with normal segments set to 0, 'L2R.nocut' for log2ratio without normal filtering, 'BAF' for B-allelic frequency values.
## annotate           logical                 If TRUE, and if such data exist in the genome cache, use it to annotate the regions content (Refseq genes for most genomes, DGV and more for homo sapiens)
## out_dir            character               Path to a folder where results will be written
## ...                ...                     Any additional parameter to read.table() to properly parse the factors file.
EaCoNiC.DiffReg <- function(SEG.cohort.obj = NULL, factors.df = NULL, samples.colname = 'Sample', factors.colname = NULL, test.2 = 'W', test.mult = 'KW', data.type = 'L2R.cut', annotate = TRUE, out_dir = getwd(), ...) {
  ## Basic checks
  if(is.null(SEG.cohort.obj)) stop('A SEG cohort object is required !')
  if(is.null(factors.df)) stop('An data frame with factors to compare is required !')
  # if(!file.exists(factors.file)) stop('Could not find the provided factors file !')
  if(length(factors.colname) == 0) stop('At least one factor to analyze is required !')
  if(!is.character(samples.colname)) stop('A column name to read sample names from the factors table is required !')
  if (!toupper(test.2) %in% c('W', 'T')) stop("Test types to analyze 2-classes factors are 'W' (wilcoxon) or 'T' (t-test).")
  if (!toupper(test.mult) %in% c('KW', 'ANOVA')) stop("Test types to analyze factors with more than 2 classes are 'KW' (Kruskal-Wallis) or 'ANOVA' (self-explanatory).")
  avail.datatypes <- c('L2R.cut', 'L2R.nocut', 'BAF')
  if (!data.type %in% avail.datatypes) stop(paste0("Data type to use should be one of : '", paste(avail.datatypes, collapse = "', '"), "' !"))
  
  ## Loading the factors file
  # factors.df <- read.table(file = factors.file, ...)
  # factors.df <- data.table::fread(input = factors.file, ...)
  
  ## Secondary checks
  if(!samples.colname %in% colnames(factors.df)) stop(paste0("Could not find the samples column name '", samples.colname, "' in the provided annotation file !"))
  sn.in.annot <- SEG.cohort.obj$samples.info$sample.name %in% factors.df[[samples.colname]]
  if(!any(sn.in.annot)) stop("The annotation table does not contain any sample name in common with the SEG cohort object !")
  fac.in.annot <- factors.colname %in% colnames(factors.df)
  if (!any(fac.in.annot)) stop("The annotation table does not contain any of the provided differential factors !")
  if(!all(fac.in.annot)) {
    message(paste0("The annotation table does not contain some of the provided differential factors : '", paste(factors.colname[!fac.in.annot], collapse = "', '"), "' !"))
    factors.colname <- factors.colname[fac.in.annot]
  }
  
  main.dir <- paste0(out_dir, '/Differential_Analysis/')
  ## Looping through factors
  for (my.fac in factors.colname) {
    message(paste0('Evaluating factor [', my.fac, '] ...'))
    ## Local copy of objects (to reduce if needed)
    cur.seg <- SEG.cohort.obj[[data.type]]
    cur.annot <- factors.df
    ## Synching objects
    seg.sn <- colnames(cur.seg)[-c(1:3)]
    seg.order <- order(seg.sn)
    cur.seg <- cur.seg[, c(1:3, seg.order+3)]
    seg.sn <- seg.sn[seg.order]
    cur.annot <- factors.df[order(cur.annot[[samples.colname]]),]
    cur.annot <- factors.df[cur.annot[[samples.colname]] %in% seg.sn,]
    cur.seg <- cur.seg[,c(rep(TRUE, 3), seg.sn %in% cur.annot[[samples.colname]])]
    cur.mat <- as.matrix(cur.seg[,-c(1:3)])
    
    ## Creating levels combinations
    cur.annot[[my.fac]] <- as.factor(cur.annot[[my.fac]])
    mylevels <- levels(cur.annot[[my.fac]])
    if(nlevels(cur.annot[[my.fac]]) < 2) stop('Provided factor has less than 2 classes !')
    combn.res <- combn(mylevels, 2)
    # all.combz <- lapply(1:ncol(combn.res), function(x) { list(combn.res[1,x], combn.res[2,x])})
    all.combz <- lapply(1:ncol(combn.res), function(x) { setNames(list(combn.res[1,x], combn.res[2,x]), nm = combn.res[,x]) })
    names(all.combz) <- vapply(1:ncol(combn.res), function(x) { paste(combn.res[, x, drop = TRUE], collapse = '_vs_') }, 'a')
    
    ### Calculating class-level median values
    class.tits <- paste(levels(cur.annot[[my.fac]]), ' (', table(cur.annot[[my.fac]]), ')', sep = '')
    class.med <- setNames(as.data.frame(lapply(levels(cur.annot[[my.fac]]), function(x) { matrixStats::rowMedians(cur.mat[,cur.annot[[my.fac]] == x], na.rm = TRUE) })), nm = paste0(levels(cur.annot[[my.fac]]), '.Med'))
    ### Calculating class-level sd
    class.sd <- setNames(as.data.frame(lapply(levels(cur.annot[[my.fac]]), function(x) { matrixStats::rowSds(cur.mat[,cur.annot[[my.fac]] == x], na.rm = TRUE) })), nm = paste0(levels(cur.annot[[my.fac]]), '.Sd'))
    
    if(length(mylevels) > 2) {
      ## Adding 'vs Other' combinations
      XvsO.combz <- sapply(as.character(mylevels), function(x) { setNames(list(x, mylevels[!mylevels == x]), nm = c(x, 'Other')) }, simplify = FALSE)
      names(XvsO.combz) <- vapply(XvsO.combz, function(x) { paste0(names(x), collapse = '_vs_') }, 'a')
      all.combz <- c(all.combz, XvsO.combz)
      ## Performing the mult test
      message(paste0("\tTesting all classes with ", test.mult, ' ...'))
      out.df <- difftest(my.mat = cur.mat, my.factor = cur.annot[[my.fac]], test.type = test.mult)
      
      ## Performing annotation
      if (annotate) {
        final.df <- EaCoNiC.regions.annotate(regions.df = out.df, genome = SEG.cohort.obj$meta$genome)
      }
      
      ## Writing results
      diff.dir <- paste0(c(main.dir, my.fac), collapse = '/')
      dir.create(diff.dir, recursive = TRUE)
      final.df <- cbind(cur.seg[,c(1:3)], Width = cur.seg$End - cur.seg$Start +1, class.med, class.sd, out.df)
      write.table(final.df, file = paste0(diff.dir, '/', paste(c(my.fac, data.type, test.mult), collapse = '_'), '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
      
    }
    
    ## Looping through combinations
    for (my.comb in names(all.combz)) {
      message(paste0("\tTesting ", my.comb, ' ...'))
      ## Reducing objects if needed
      mini.classes <- names(all.combz[[my.comb]])
      mini.samp.idx <- cur.annot[[my.fac]] %in% unlist(all.combz[[my.comb]])
      mini.annot <- cur.annot[mini.samp.idx,]
      mini.seg <- cur.seg[,c(rep(TRUE, 3), colnames(cur.seg)[-c(1:3)] %in% mini.annot[[samples.colname]])]
      mini.mat <- cur.mat[, mini.samp.idx]
      ## Handling 'Other'
      if (mini.classes[2] == 'Other') {
        newlevel <- rep('Other', nlevels(mini.annot[[my.fac]]))
        newlevel[levels(mini.annot[[my.fac]]) == mini.classes[1]] <- mini.classes[1]
        levels(mini.annot[[my.fac]]) <- newlevel
      }
      ## Computing median L2Rs and Sds
      mini.med <- setNames(as.data.frame(lapply(mini.classes, function(x) { matrixStats::rowMedians(mini.mat[,mini.annot[[my.fac]] == x], na.rm = TRUE) })), nm = paste0(mini.classes, '.Med'))
      mini.sd <- setNames(as.data.frame(lapply(mini.classes, function(x) { matrixStats::rowSds(mini.mat[,mini.annot[[my.fac]] == x], na.rm = TRUE) })), nm = paste0(mini.classes, '.Sd'))
      
      ## Performing the statistical test
      out.df <- difftest(my.mat = mini.mat, my.factor = mini.annot[[my.fac]], test.type = test.2)
      message(str(out.df))
      final.df <- cbind(cur.seg[,c(1:3)], Width = cur.seg$End - cur.seg$Start +1, mini.med, DIFF = mini.med[,1] - mini.med[,2], mini.sd, out.df)
      colnames(final.df)[colnames(final.df) == 'DIFF'] <- paste0(c(my.comb, data.type, 'Med.Diff'), collapse ='.')
      
      ## Performing annotation
      if (annotate) {
        # final.df <- EaCoN::EaCoN.region.annotate(my.df = final.df, genome = SEG.cohort.obj$meta$genome)
        final.df <- EaCoNiC.regions.annotate(regions.df = final.df, genome = SEG.cohort.obj$meta$genome)
      }
      
      
      # diff.dir <- paste0(c(main.dir, my.fac), collapse = '/')
      # dir.create(diff.dir, recursive = TRUE)
      # final.df <- cbind(cur.seg[,c(1:3)], Width = cur.seg$End - cur.seg$Start +1, class.med, class.sd, out.df)
      # write.table(final.df, file = paste0(diff.dir, '/', paste(c(my.fac, data.type, test.mult), collapse = '_'), '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
      # 
      
      ## Writing results
      diff.dir <- paste0(c(main.dir, my.fac, my.comb), collapse = '/')
      dir.create(diff.dir, recursive = TRUE)
      write.table(final.df, file = paste0(diff.dir, '/', paste(c(my.fac, my.comb, data.type, test.2), collapse = '_'), '.tsv'), sep = "\t", quote = FALSE, row.names = FALSE)
      ### Filtered
      out.df <- final.df[!is.na(final.df$AdjP),]
      write.table(out.df, file = paste0(diff.dir, '/', paste(c(my.fac, my.comb, data.type, test.2, 'filtered'), collapse = '_'), '.tsv'), sep = "\t", quote = FALSE, row.names = FALSE)
      
      ## Volcano plot
      pcuts <- c(.05, .01)
      png(filename = paste0(diff.dir, '/', paste(c(my.fac, my.comb, data.type, test.2, 'volcano.png'), collapse = '_')), width = 1024, height = 768)
      # plot(final.df[,7], -log10(final.df$RawP), xlim = max(abs(final.df[,7]))*c(-1,1), pch = 20, col = densCols(final.df[,7], -log10(final.df$RawP), bandwidth = .05), main = sub(pattern = '_vs_', replacement = ' vs ', x = paste0(my.fac, ' : ', my.comb, "\nVolcano plot")), xlab = 'L2R Difference', ylab = '-log10(RawP)')
      plot(out.df[,7], -log10(out.df$RawP), xlim = max(abs(out.df[,7]))*c(-1,1), pch = 20, col = densCols(out.df[,7], -log10(out.df$RawP), bandwidth = .05), main = sub(pattern = '_vs_', replacement = ' vs ', x = paste0(my.fac, ' : ', my.comb, "\nVolcano plot")), xlab = 'L2R Difference', ylab = '-log10(RawP)')
      abline(v = 0, col = 4, lwd = 3, lty = 2)
      abline(h = -log10(pcuts), col = c('brown', 2), lwd = 3, lty = 3)
      # text(max(abs(final.df[,7]))*.95, (-log10(pcuts)), labels = c(length(which(final.df$RawP < pcuts[1])), length(which(final.df$RawP < pcuts[2]))), col = c('brown', 2), pos = 3)
      text(max(abs(out.df[,7]))*.95, (-log10(pcuts)), labels = c(length(which(out.df$RawP < pcuts[1])), length(which(out.df$RawP < pcuts[2]))), col = c('brown', 2), pos = 3)
      dev.off()
    }
  }
}


## [NOT EXPORTED] Perform the differential test
difftest <- function(my.mat = NULL, my.factor = NULL, test.type = 'W') {
  
  test.type <- toupper(test.type)
  
  test2function <- list(
    'ANOVA' = c('stats::lm', 'F'),
    'KW' = c('stats::kruskal.test', 'KW'),
    'T' =  c('stats::t.test', 'T'),
    'W' = c('stats::wilcox.test', 'W'))
  
  ## Check
  if (!test.type %in% names(test2function)) stop('Unsupported test !')
  table(my.factor)
  
  stat.slotname <- if(test.type == 'ANOVA') 'F value' else 'statistic'
  pval.slotname <- if(test.type == 'ANOVA') 'Pr(>F)' else 'p.value'
  
  ## Getting the test function
  func.split <- unlist(strsplit(test2function[[test.type]][1], '::'))
  test.function <- base::get(func.split[2], envir = base::loadNamespace(func.split[1]))
  
  ## Performing the test
  out.df <- data.frame(Stat = rep(NA, nrow(my.mat)))
  out.df$RawP <- out.df$Stat
  for (k in 1:nrow(my.mat)) {
    # message(k)
    tt <- try(test.res <- test.function(as.numeric(my.mat[k,]) ~ my.factor), silent = TRUE)
    # message(k, is(tt))
    if (is(tt) != 'try-error') {
      if (test.type == 'ANOVA') {
        test.res <- stats::anova(test.res)
      }
      out.df$Stat[k] <- test.res[[stat.slotname]][1]
      out.df$RawP[k] <- test.res[[pval.slotname]][1]
    }
  }
  colnames(out.df)[1] <- test2function[[test.type]][2]
  out.df$AdjP <- stats::p.adjust(p = out.df$RawP, method = 'BH')
  return(out.df)
}

## Clusterized Heatmap
## SEG.cohort.obj       Robject                         The output from the EaCoNiC.SEG.aggregate() function.
## data.type            character                       The data slot to use. Should be one of 'L2R.cut', 'L2R.nocut', 'BAF'.
## chr.discard          vector of characters (or NULL)  Chromosomes to discard before clustering (useful to get rid of gonosomes for genomes with gender)
## dist.method          character                       The method to compute the distance bewteen samples. Should be one of amap::Dist().
## agglo.method         character                       The agglomeration method to apply to the distance matrix. Should be one of stats::hclust(). 
## factors.df           data frame (or NULL)             A data frame that contains at least two columns : one with samplenames, the other(s) being factors to use to perform the diff tests
## samples.colnames     character (or NULL)             Column name of factors.file to identify samples
## factors.colname      vector of characters (or NULL)  Column name(s) of factors.file to specify classes to compare
## chr.simplify         logical                         Simplify the chromosome names for plots by removing a potential common root (by example, 'chr' for homo sapiens).
## target.bins          numeric (or NULL)               To simulate the effect of the size of the region in the distance computation, larger regions are repeated multiple times (bins). The value of target.bins is the total number of bins : the function will automatically multiply the regions to 
EaCoNiC.heatmap <- function(SEG.cohort.obj = NULL, data.type = 'L2R.cut', chr.discard = NULL, dist.method = 'pearson', agglo.method = "ward.D", amp.cut = 1.5, del.cut = -1.5, factors.df = NULL, samples.colname = 'SampleName', factors.colname = NULL, chr.simplify = TRUE, target.bins = NULL, out_dir = getwd(), ...) {
  ## Basic checks
  ### SEG object
  if(is.null(SEG.cohort.obj)) stop('A SEG cohort object is required !')
  ### Factors file
  if(is.null(factors.df)) message('No factor table provided.') else {
    if(length(factors.colname) == 0) stop('A factor table was provided, but no factor to plot !')
    if(!is.character(samples.colname)) stop('A factor table was provided, but no column name to identify sample names !')
  }
  ### Data type
  avail.datatypes <- c('L2R.cut', 'L2R.nocut', 'BAF')
  if(!data.type %in% avail.datatypes) stop(paste0("Data type to use should be one of : '", paste(avail.datatypes, collapse = "', '"), "' !"))
  ### Chromosomes to discard
  avail.chr <- levels(SEG.cohort.obj[[data.type]]$Chr)
  if(!is.null(chr.discard)) if (!all(chr.discard %in% avail.chr)) stop('At least one of requested chromosomes to discard is not available for this cohort !')
  ## Output dir
  if(!dir.exists(out_dir)) stop('The specified output directory does not exist !')
  ## Target bins
  # if(!is.null(target.bins) & !is.numeric(target.bins)) stop('The target number of width bins should be an integer !')
  if(!is.null(target.bins) & !is.numeric(target.bins)) stop('Target number of width bins should be NULL or a positive integer')
  if(target.bins < 1) stop('The target number of width bins should be a positive integer !')
  
  ## Loading the factors file
  if (!is.null(factors.df)) {
    # factors.df <- read.table(file = factors.file, ...)
    ## Factors-specific checks
    if(!samples.colname %in% colnames(factors.df)) stop(paste0("Could not find the samples column name '", samples.colname, "' in the provided annotation file !"))
    sn.in.annot <- SEG.cohort.obj$samples.info$sample.name %in% factors.df[[samples.colname]]
    if(!any(sn.in.annot)) stop("The annotation table does not contain any sample name in common with the SEG cohort object !")
    fac.in.annot <- factors.colname %in% colnames(factors.df)
    if (!any(fac.in.annot)) stop("The annotation table does not contain any of the provided differential factors !")
    if(!all(fac.in.annot)) {
      message(paste0("The annotation table does not contain some of the provided differential factors : '", paste(factors.colname[!fac.in.annot], collapse = "', '"), "' !"))
      factors.colname <- factors.colname[fac.in.annot]
    }
  }
  
  my.df <- SEG.cohort.obj[[data.type]]
  
  ## Limiting chromosomes
  if(!is.null(chr.discard)) {
    message(paste0('Discarding chromosome(s) [', paste(chr.discard, collapse = ', '), '] ...'))
    my.df <- my.df[!my.df$Chr %in% chr.discard,]
  }
  nreg.ori <- nrow(my.df)
  
  ## Weighting by segment width
  if (!is.null(target.bins)) {
    message("Weighting by region width ...")
    if (nreg.ori >= target.bins) message(paste0('Target width bins (', target.bins, ') is smaller or equal than identified regions (', nreg.ori, '), thus no width weighting could be applied.')) else {
      my.w <- my.df$End - my.df$Start +1
      gen.length <- sum(vapply(unique(my.df$Chr), function(x) { max(my.df$End[my.df$Chr == x]) }, 1L))
      bin.size.first <- round(gen.length / target.bins)
      smalliz <- my.w <= bin.size.first
      smalliz.wsum <- sum(my.w[smalliz])
      target.bins <- target.bins + length(which(!smalliz))
      bin.size.second <- round((gen.length - smalliz.wsum) / target.bins)
      seg.weight <- ceiling(my.w / bin.size.second)
      my.df <- my.df[rep(seq_len(nrow(my.df)), seg.weight),]
      message(paste0("\tGenerated ", sum(seg.weight), ' bins from ', nreg.ori, " regions.\n\t", length(which(smalliz)), ' of these regions were not split to bins, due to their size below the average bin size (', bin.size.second, ').'))
    }
  } else {
    message('No target.bins value specified : all genomic regions will be used as is.')
    target.bins <- nreg.ori
  }
  
  ## Converting df to a matrix
  my.mat <- as.matrix(my.df[,-c(1:3)])
  
  
  ## Preparing graphical split on chromosomes
  chr.split <- as.factor(my.df$Chr)
  if (chr.simplify) levels(chr.split) <- chr.simplifier(levels(chr.split))
  
  ## Palette
  zlim = c(-1, 0, 1)
  myPalette <- c("red", "grey85", "blue")
  myRamp <- circlize::colorRamp2(zlim, myPalette)
  
  
  ## Looping on distances
  for (myd in dist.method) {
    
    ## Output dir
    hch.dir <- paste0(out_dir, '/Clusterized_Heatmap/', data.type, '/', myd, '_', agglo.method, if(!is.null(chr.discard)) paste0('/NO.', paste(chr.discard, collapse = '.'), if(is.null(target.bins)) '' else '_ww', '/') else '/ALL_chromosomes/')
    dir.create(hch.dir, recursive = TRUE)
    
    ## Clustering
    message('Performing the hierarchical clustering ...')
    my.dist <- amap::Dist(x = t(my.mat), method = myd, nbproc = 1)
    ### Dealing with methods that dislike missing values
    if(any(is.na(my.dist))) {
      my.mat.tmp <- my.mat
      ### Replacing NA values by row median
      for(x in seq_len(nrow(my.mat))) my.mat.tmp[x,is.na(my.mat[x,])] <- median(my.mat[x,], na.rm = TRUE)
      my.dist <- amap::Dist(x = t(my.mat.tmp), method = dist.method, nbproc = 1)
    }
      
    my.hc <- stats::hclust(d = my.dist, method = agglo.method)
    
    ## Heatmap (landscape)
    message('Plotting the landscape heatmap ...')
    hal <- if (!is.null(factors.df)) ComplexHeatmap::rowAnnotation(df = factors.df[,factors.colname], show_annotation_name = FALSE) else NULL
    hml <- suppressMessages(ComplexHeatmap::Heatmap(
      matrix = t(my.mat)
      , name = data.type
      , show_column_names = FALSE
      , col = myRamp
      , cluster_rows = my.hc
      , cluster_columns = FALSE
      , row_title = SEG.cohort.obj$meta$cohort.name
      , row_names_max_width = ComplexHeatmap::max_text_width(rownames(my.mat))*10
      , left_annotation = hal
      , column_split = chr.split
      , raster_device = 'png'
      , row_dend_width = grid::unit(4, "cm")
      , column_dend_height = grid::unit(4, "cm")
    ))
    ## Writing
    
    png(filename = paste0(hch.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, data.type, myd, agglo.method, 'heatmap'), collapse = '_'), '.png'),
        # width = nrow(my.mat)/3 + if(!is.null(factors.colname)) 300 else 0,
        width = max(nreg.ori/3 + if(!is.null(factors.colname)) 300 else 0, 1800),
        height = max(ncol(my.mat)*35, 700)
    )
    ComplexHeatmap::draw(hml)
    dev.off()
    
    ## HC split membership
    hc.split <- cutree(my.hc, k = 2:ncol(my.mat))
    colnames(hc.split) <- paste0('HC.', colnames(hc.split))
    write.table(data.frame(Sample = rownames(hc.split), hc.split), file = paste0(hch.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, data.type, myd, agglo.method, 'membership'), collapse = '_'), '.txt'), quote = FALSE, row.names = FALSE, sep = "\t")
  }
}


## Auto-correlation map
EaCoNiC.autocor <- function(SEG.cohort.obj = NULL, data.type = 'L2R.cut', cor.method = 'spearman', chr.discard = NULL, chr.simplify = TRUE, out_dir = getwd()) {
  ## Basic checks
  if(is.null(SEG.cohort.obj)) stop('A SEG cohort object is required !')
  ### Data type
  avail.datatypes <- c('L2R.cut', 'L2R.nocut', 'BAF')
  if(!data.type %in% avail.datatypes) stop(paste0("Data type to use should be one of : '", paste(avail.datatypes, collapse = "', '"), "' !"))
  ### Chromosomes to discard
  avail.chr <- levels(SEG.cohort.obj[[data.type]]$Chr)
  if(!is.null(chr.discard)) if (!all(chr.discard %in% avail.chr)) stop('At least one of requested chromosomes to discard is not available for this cohort !')
  ## Output dir
  if(!dir.exists(out_dir)) stop('The specified output directory does not exist !')
  
  ## Getting data
  my.df <- SEG.cohort.obj[[data.type]]
  
  ## Limiting chromosomes
  if(!is.null(chr.discard)) {
    message(paste0('Discarding chromosome(s) [', paste(chr.discard, collapse = ', '), '] ...'))
    my.df <- my.df[!my.df$Chr %in% chr.discard,]
  }
  
  ## Converting to matrix
  my.mat <- as.matrix(my.df[,-c(1:3)])
  rownames(my.mat) <- seq_len(nrow(my.mat))
  
  ## Removing genomic regions with no variance
  rV <- matrixStats::rowVars(x = my.mat, na.rm = TRUE)
  rV.idx <- rV > 0
  my.mat <- my.mat[rV.idx,]
  ## Computing correlation matrix
  message('Computing region correlations ...')
  mat2cor <- my.mat
  if(any(is.na(matrixStats::colVars(my.mat)))) {
    ### Dealing with methods that dislike missing values
    message('Missing value(s) found. Cheating with random noise matrix...')
    mat.sd <- sd(as.vector(my.mat), na.rm = TRUE)
    noise.mat <- matrix(data = rnorm(n = length(my.mat), mean = 0, sd = mat.sd/1000), nrow = nrow(my.mat), ncol = ncol(my.mat))
    mat2cor[is.na(mat2cor)] <- 0
    mat2cor <- mat2cor + noise.mat
  }
  my.cor <- stats::cor(t(mat2cor), method = cor.method)
  
  # ## Handling NA values
  # if(any(is.na(my.cor))) {
  #   my.mat.tmp <- my.mat
  #   for(x in seq_len(nrow(my.mat))) my.mat.tmp[x,is.na(my.mat[x,])] <- median(my.mat[x,], na.rm = TRUE)
  #   my.cor <- stats::cor(t(my.mat.tmp), method = cor.method)
  #   rm(my.mat.tmp)
  # }
  
  ## Chromosomes cut breakpoints
  kbr <- vapply(sort(unique(my.df$Chr[rV.idx])), function(k) max(which(my.df$Chr[rV.idx] == k)), 1)
  ## Correlation matrix segmentation
  message('Segmenting the correlation map ...')
  max.blocks <- min(floor(nrow(my.cor)/20), 500)
  bs.res <- blockseg::blockSeg(Y = my.cor, max.break = max.blocks, verbose = FALSE)
  ## Computing median correlation scores per block
  ### Getting actual breaks
  my.breaks <- sort(unique(c(bs.res@RowBreaks[[length(bs.res@RowBreaks)]], bs.res@ColBreaks[[length(bs.res@ColBreaks)]], kbr, nrow(my.cor))))
  
  ##Filtering contiguous breaks
  message('Filtering contiguous breakpoints ...')
  runs <- unname(split(seq_along(my.breaks), cumsum(c(0, diff(my.breaks) > 1))))
  streams.idx <- unlist(lapply(runs[lengths(runs) > 1], range))
  unik.idx <- unlist(runs[which(vapply(runs, length, 1L) ==1)])
  my.breaks2 <- sort(unique(my.breaks[sort(c(streams.idx, unik.idx))], nrow(my.cor)))
  
  
  ## Preparing graphical split on chromosomes
  chr.split <- as.factor(my.df$Chr[rV.idx])
  if (chr.simplify) levels(chr.split) <- chr.simplifier(levels(chr.split))
  
  ## Preparing graphical split correlation map blocks
  mb.diff <- diff(my.breaks2)
  mb.diff <- c(my.breaks2[1], mb.diff)
  blocks.split <- as.factor(rep(seq_along(mb.diff), times = mb.diff))
  
  ## Generating output dir
  cor.dir <- paste0(out_dir, '/Auto_Correlation/', data.type, '/')
  dir.create(path = cor.dir, recursive = TRUE)
  ## Saving the correlation map
  saveRDS(my.cor, file = paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'cormap'), collapse = '_'), '.RDS'), compress = 'bzip2')
  # writexl::write_xlsx(x = as.data.frame(my.cor), path = paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'cormap'), collapse = '_'), '.xlsx'), col_names = FALSE, format_headers = FALSE)
  
  ## Palette
  zlim = c(-1, 0, 1)
  myPalette <- c("red", "grey85", "blue")
  myRamp <- circlize::colorRamp2(zlim, myPalette)
  
  ## Plotting the correlation map
  set.seed(12345)
  hacR <- ComplexHeatmap::rowAnnotation(df = data.frame(Chr = chr.split), name = 'Chr')
  set.seed(12345)
  hacC <- ComplexHeatmap::columnAnnotation(df = data.frame(Chr = chr.split), name = 'Chr')
  
  hmp <- suppressMessages(ComplexHeatmap::Heatmap(matrix = my.cor
                                                  , name = cor.method
                                                  , show_row_names = FALSE
                                                  , show_column_names = FALSE
                                                  , col = myRamp
                                                  , cluster_rows = FALSE
                                                  , cluster_columns = FALSE
                                                  , row_split = blocks.split
                                                  , column_split = blocks.split
                                                  , column_gap = grid::unit(.3, 'mm')
                                                  , row_gap = grid::unit(.3, 'mm')
                                                  , right_annotation = hacR
                                                  , bottom_annotation = hacC
                                                  , border = TRUE
                                                  , raster_device = 'png'
  ))
  
  png(paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'cormap_blocks'), collapse = '_'), '.png'), width = 2000 + (length(my.breaks2)), height = 2000)
  ComplexHeatmap::draw(hmp)
  dev.off()
  
  ## Aggregate blocks
  `%do%` <- foreach::`%do%`
  agreg.tmp <- t(foreach:::foreach(b = levels(blocks.split), .combine = 'cbind') %do% { as.matrix(aggregate(my.cor[,blocks.split == b], list(block = blocks.split), median)[,-1]) })
  agreg.res <- pracma::rot90(foreach:::foreach(b = seq_len(ncol(agreg.tmp)), .combine = 'cbind') %do% { as.matrix(aggregate(agreg.tmp[,b], list(block = blocks.split), median)[,-1]) })
  start.idx <- c(1, my.breaks2[-length(my.breaks2)]+1)
  agreg.pos <- data.frame(Chr = my.df$Chr[rV.idx][my.breaks2], Start = my.df$Start[rV.idx][start.idx], End = my.df$End[rV.idx][my.breaks2])
  colnames(agreg.res) <- rownames(agreg.res) <- paste0(agreg.pos$Chr, ':', agreg.pos$Start, '-', agreg.pos$End)
  
  # agreg.annot <- EaCoN.regions.annotate(regions.df = agreg.pos, genome = 'mm10')
  # agreg.out <- data.frame(agreg.pos, agreg.res, agreg.annot[,-c(1:3)], check.names = FALSE)
  
  agreg.out <- data.frame(agreg.pos, agreg.res, check.names = FALSE)
  writexl::write_xlsx(x = agreg.out, path = paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'cormap_blocks'), collapse = '_'), '.xlsx'), col_names = TRUE, format_headers = TRUE)
}


## Subset
EaCoNiC.SEG.subset <- function(SEG.cohort.obj = NULL, samples.keep = NULL, chr.keep = NULL, cohort.name = NULL, out_dir = getwd(), return.data = FALSE){
  
}


## test
# EaCoNiC.set.cache(cache.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/EaCoNiC_cache')

# EaCoNiC.DiffReg(SEG.cohort.obj = readRDS('/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/TEST2.SEG.COHORT.RDS'), factors.file = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/samples_info.txt', samples.colname = 'sn', factors.colname = c('Disease_control_benefit', 'BOR'), test.2 = 'W', test.mult = 'KW', data.type = 'L2R.cut', annotate = TRUE, out_dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev', header = TRUE, sep = "\t")

# EaCoNiC.heatmap(SEG.cohort.obj = readRDS('/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/B20067_FRDA_05.SEG.COHORT.RDS'), data.type = 'L2R.cut', factors.file = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/samples_info.txt', samples.colname = 'sn', factors.colname = c('Disease_control_benefit', 'BOR'), dist.method = 'pearson', agglo.method = 'ward.D', amp.cut = 1.5, loss.cut = -1.5, out_dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev', header = TRUE, sep = "\t", chr.discard = c('chrX', 'chrY'), target.bins = 1E+05, chr.simplify = TRUE)

# EaCoNiC.autocor <- function(SEG.cohort.obj = readRDS('/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/B20067_FRDA_05.SEG.COHORT.RDS'), data.type = 'L2R.cut', chr.discard = NULL, cor.method = 'spearman', chr.simplify = TRUE, out_dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev')
