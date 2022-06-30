## PACKAGES USED
## CRAN : tibble, dplyr, cairoDevice
## BIOCONDUCTOR : GenomeInfoDb, GenomicRanges, rtracklayer, ComplexHeatmap
## GITHUB : jchiquet/blockseg



# SEG.rds.list <- list.files(path = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/CNA', pattern = '*.SEG.*.RDS$', full.names = TRUE, recursive = TRUE)

## User function to set the cache directory
EaCoN.set.cache <- function(cache.dir = NULL) {
  if(is.null(cache.dir)) {
    message('No directory specified :')
    if ('EaCoN.cache' %in% utils::installed.packages()) {
      message("\tAn installation of 'EaCoN.cache' was found and will be used.")
      cache.dir <- paste0(base::tempdir(check = TRUE), '/EaCoN.cache/')
    } else {
      cache_dir <- tools::R_user_dir("recount3", "cache")
      message("\tA cache folder will be used.")
      BiocFileCache::BiocFileCache(cache_dir)
    }
  }
  if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
  message(paste0('EaCoN cache directory set to [', cache.dir, ']'))
  # options(EACON_CACHE=cache.dir)
}

## Get annotations from the UCSC and cache them locally (works for multiple genomes)
EaCoN.cache.annotation <- function(genome = NULL, cache.dir = NULL) {
  message('An internet connexion is required.')
  if(is.null(cache.dir)) stop('A cache directory to write files is required !')
  if(!file.exists(cache.dir)) stop(paste0('Cache directory [', cache.dir, '] does not exist !'))
  
  ## Setting the cache dir
  # EaCoN::EaCoN.set.cache(cache.dir = cache.dir)
  EaCoN.set.cache(cache.dir = cache.dir)
  
  ## Getting the list of available registered genomes @ UCSC (from GenomeInfoDb)
  avail.genz <- sort(GenomeInfoDb::registered_UCSC_genomes()$genome)
  
  if (is.null(genome)) {
    message(paste0('No genome specified : all ', length(avail.genz), ' registered UCSC genomes will be used. This may take a while !'))
    genome <- avail.genz
  }
  genome <- genome[genome %in% avail.genz]
  if(length(genome) == 0) stop('None of the specified genome(s) is available @ UCSC !')
  
  annot.type <- list(
    cytoBandIdeo = c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'),
    refGene = c('chrom', 'txStart', 'txEnd', 'strand', 'exonCount', 'name2'),
    cpgIslandExt = c('chrom', 'chromStart', 'chromEnd', 'name', 'perCpG', 'perGc', 'obsExp'),
    genomicSuperDups = c('chrom', 'chromStart', 'chromEnd', 'strand'),
    dgvMerged = c('chrom', 'chromStart', 'chromEnd', 'strand', 'varType'))
  
  for (gen in genome) {
    annot.obj <- setNames(vector(mode = 'list', length = length(annot.type)), names(annot.type))
    message(paste0('Getting annotations for ', gen, ' ...'))
    ## Trying to get chromosomes info
    trycI <- try(cI <- GenomeInfoDb::getChromInfoFromUCSC(genome = gen, assembled.molecules.only = TRUE, as.Seqinfo = TRUE), silent = TRUE)
    ## Cycling in annotation tables
    for (tt in names(annot.type)) {
      ## Querying
      tryres <- try(query <- rtracklayer::ucscTableQuery(gen, table = tt), silent = TRUE)
      if (class(tryres) == 'try-error') {
        message(paste0("\tGenome [", gen, '] did NOT have ', tt, '.'))
      } else {
        df <- rtracklayer::getTable(query)
        message(paste0("\tRetrieved ", tt, '.'))
        ## Limiting to desired columns
        df <- df[, colnames(df) %in% annot.type[[tt]]]
        ## Converting to GRanges
        gr <- GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, ignore.strand = FALSE, seqnames.field = annot.type[[tt]][1], start.field = annot.type[[tt]][2], end.field = annot.type[[tt]][3], starts.in.df.are.0based = TRUE)
        ## Adding seqinfo when available
        if (!class(trycI) == 'try-error') {
          GenomeInfoDb::seqlevels(gr, pruning.mode = 'coarse') <- GenomeInfoDb::seqlevels(cI)
          GenomeInfoDb::seqinfo(gr) <- cI
        }
        ## Merging UCSC refGenes
        if (tt == 'refGene') {
          `%>%` <- dplyr::`%>%`
          gr <- suppressMessages(gr %>% tibble::as_tibble() %>% dplyr::group_by(seqnames, strand, name2) %>% dplyr::summarize(start = min(start), end = max(end)) %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE))
        }
        ## Saving in the object
        annot.obj[[tt]] <- gr
      }
    }
    ## Saving is there is something to save
    if(!all(sapply(X = annot.obj, is.null))) save(object = annot.obj, file = paste0(cache.dir, '/', gen, '.rda'), compress = 'bzip2')
  }
}

## [NOT EXPORTED] Get cached genome
get.cached.genome <- function(genome = 'hg19') {
  ### Getting cache directory if available in options()
  cache.dir <- getOption('EACON_CACHE')
  ### Checking if a cache file exists
  cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.rda'))
  if (is.null(cache.dir) | !cache.exists) {
    ## If there is no cache directory in options()
    message('No cache directory found, or no cache found for provided genome : Annotation will be retrieved to a temporary directory (an internet connection is required).')
    cache.dir <- paste0(tempdir(check = TRUE), '/EaCoN_cache/')
    dir.create(cache.dir)
    # EaCoN::EaCoN.cache.annotation(cache.dir = cache.dir, genome = genome)
    EaCoN.cache.annotation(cache.dir = cache.dir, genome = genome)
  }
  ## Loading cache
  message(paste0('Loading cache for genome ', genome, ' ...'))
  load(paste0(cache.dir, '/', genome, '.rda'), envir = .GlobalEnv)
}

## [NOT EXPORTED] Convert a cached dataset of chromosome data from a specified genome to a more complex object
genome2cs <- function(genome = NULL) {
  ## Checks
  ### Getting cache directory if available in options()
  cache.dir <- getOption('EACON_CACHE')
  ### Checking if a cache file exists
  cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.rda'))
  if (is.null(cache.dir) | !cache.exists) {
    ## If there is no cache directory in options()
    message('No cache directory found, or no cache found for provided genome : Annotation will be retrieved to a temporary directory (an internet connection is required).')
    cache.dir <- paste0(tempdir(check = TRUE), '/EaCoN_cache/')
    dir.create(cache.dir)
    # EaCoN::EaCoN.cache.annotation(cache.dir = cache.dir, genome = genome)
    EaCoN.cache.annotation(cache.dir = cache.dir, genome = genome)
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
get.cytobands <- function(genome = NULL, uncertain.filter = TRUE) {
  ## Checks
  ### Getting cache directory if available in options()
  cache.dir <- getOption('EACON_CACHE')
  ### Checking if a cache file exists
  cache.exists <- file.exists(paste0(cache.dir, '/', genome, '.rda'))
  if (is.null(cache.dir) | !cache.exists) {
    ## If there is no cache directory in options()
    message('No cache directory found, or no cache found for provided genome : Annotation will be retrieved to a temporary directory (an internet connection is required).')
    cache.dir <- paste0(tempdir(check = TRUE), '/EaCoN_cache/')
    dir.create(cache.dir)
    # EaCoN::EaCoN.cache.annotation(cache.dir = cache.dir, genome = genome)
    EaCoN.cache.annotation(cache.dir = cache.dir, genome = genome)
  }
  ## Loading cache
  message('Retrieving cytobands from cache ...')
  load(paste0(cache.dir, '/', genome, '.rda'))
  if (is.null(annot.obj$cytoBandIdeo)) {
    message("\tRequested genome did not have cytobands data !")
    return(NULL)
  }
  gr <- annot.obj$cytoBandIdeo
  rm(annot.obj)
  if(uncertain.filter) {
    if (!'gieStain' %in% colnames(GenomicRanges::mcols(gr))) {
      message("\tRequested genome did not have staining, so filering of uncertain cytobands was not possible.")
    } else {
      bandcolz <- biovizBase::getBioColor(type = 'CYTOBAND')
      unc.bands <- gr$gieStain %in% names(bandcolz)[c(2,3,5)]
      message(paste0("\tFiltering ", length(which(unc.bands)), ' uncertain bands ...'))
      gr <- gr[!unc.bands,]
    }
  }
  return(gr)
}

## Simplify a vector of chr names by identifying and deleting its common root 
chr.simplifier <- function(chr.names = NULL) {
  message('Simplfying chromosome names ...')
  cs.decomp <- lapply(chr.names, function(x) { unlist(strsplit(x, split = '')) })
  cs.len <- sapply(cs.decomp, length)
  cs.com <- lapply(seq_len(min(cs.len)), function(x) {
    unique(vapply(cs.decomp, function(y) { y[x] }, 'a'))
  })
  com.key <- paste(vapply(cs.com, function(x) { if(length(x) == 1) x else ''}, 'a'), collapse = '')
  if (com.key == '') {
    return(chr.names)
  } else return(gsub(pattern = com.key, replacement = '', levels(chr.split)))
}

## Aggregate samples to form a cohort 
EaCoNiC.SEG.aggregate <- function(SEG.rds.list = list.files(path = getwd(), pattern = '*.SEG.*.RDS*', full.names = TRUE, recurisve = TRUE), L2R = TRUE, BAF = TRUE, cohort.name = 'TEST', out.dir = getwd(), return.data = FALSE) {
  ## Early checks
  if (is.null(SEG.rds.list)) stop("A list of *SEG*.RDS files (results from EaCoN::Segment()) is required !")
  if (length(SEG.rds.list) < 2) stop("To build a cohort, at least two samples are required !")
  if (!any(c(L2R, BAF))) stop("At least one of 'L2R' and 'BAF' should be set to TRUE.")
  if (is.null(out.dir) & !return.data) stop("No output directory specified, and 'return.data' is set to FALSE : as nothing will be output, nothing will be performed !")
  if (is.null(cohort.name)) stop("A cohort name is required !")
  
  ## Loading data : First pass to control data compatibility
  message('Reading samples metadata to ensure compatibility ...')
  samp.meta.list <- lapply(seq_along(SEG.rds.list), function(x) {
    message(paste0("\tParsing metadata from file ", x, ' / ', length(SEG.rds.list), ' ...'))
    sobj <- readRDS(SEG.rds.list[x])
    return(list(sample.name = sobj$meta$basic$samplename, loss.cut = sobj$meta$eacon[['L2R-segments-loss-cutoff']], gain.cut = sobj$meta$eacon[['L2R-segments-gain-cutoff']], genome =  sobj$meta$basic$genome, segmenter = sobj$meta$eacon$segmenter))
  })
  samp.info <- as.data.frame(t(as.data.frame(lapply(samp.meta.list, unlist), check.names = FALSE)))
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
      samp.chr.ov <- GenomicRanges::findOverlaps(query = unl.l2r.gr, subject = IRanges::reduce(chr.gr))
      samps.l2r.gr <- unl.l2r.gr[S4Vectors::from(samp.chr.ov)]
    }
    
    ## Synching to a DF
    message('Merging all breakpoints ...')
    reg.l2r.df <- cbind(as.data.frame(samps.l2r.gr)[,1:3], as.data.frame(lapply(seq_along(l2r.gr.list), function(x) {
      message(paste0("\tAdding ", samp.info$sample.name[x], '...'))
      GenomicRanges::mcols(samps.l2r.gr) <- S4Vectors::DataFrame(L2R = NA)
      over.gr <- IRanges::findOverlaps(samps.l2r.gr, l2r.gr.list[[x]])
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
      reg.l2r.df.filt.cut[,x][reg.l2r.df.filt.cut[,x] > samp.info$loss.cut[(x-3)] & reg.l2r.df.filt.cut[,x] < samp.info$gain.cut[(x-3)]] <- 0
    }
  }
  
  ### BAF )
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
      samp.chr.ov <- GenomicRanges::findOverlaps(query = unl.baf.gr, subject = IRanges::reduce(chr.gr))
      samps.baf.gr <- unl.baf.gr[S4Vectors::from(samp.chr.ov)]
    }
    
    ## Getting all possible intervals from the cohort merge.
    # unl.baf.gr <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(c(baf.gr.list, chr = chr.gr))))
    
    ## Synching to a DF
    message('Merging all breakpoints ...')
    reg.baf.df <- cbind(as.data.frame(samps.baf.gr)[,1:3], as.data.frame(lapply(seq_along(baf.gr.list), function(x) {
      message(paste0("\tAdding ", samp.info$sample.name[x], '...'))
      GenomicRanges::mcols(samps.baf.gr) <- S4Vectors::DataFrame(BAF = NA)
      over.gr <- IRanges::findOverlaps(samps.baf.gr, baf.gr.list[[x]])
      samps.baf.gr$BAF[over.gr@from] <- baf.gr.list[[x]]$BAF[over.gr@to]
      return(samps.baf.gr$BAF)
    }), col.names = samp.info$sample.name, check.names = FALSE, stringsAsFactors = FALSE))
    colnames(reg.baf.df)[1:3] <- c('Chr', 'Start', 'End')
    
    ## Filter intervals with only NAs
    reg.baf.df.filt <- reg.baf.df[!is.na(matrixStats::rowAnys(as.matrix(reg.baf.df[,-c(1:3)]))),]
  }
  out.obj <- list(samples.info = samp.info, L2R.cut = reg.l2r.df.filt, L2R.nocut = reg.l2r.df.filt.cut, BAF = reg.baf.df.filt, meta = list(type = 'SEG', cohort.name = cohort.name, genome = my.genome, sample.files = SEG.rds.list))
  if (!is.null(out.dir)) saveRDS(object = out.obj, file = paste0(out.dir, '/', cohort.name, '.SEG.COHORT.RDS'), compress = 'bzip2')
  if (return.data) return(out.obj) else return(invisible(NULL))
}

# EaCoN.cache.annotation(cache.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/EaCoNiC_cache', genome = 'hg19')

# EaCoN.set.cache(cache.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/EaCoN_cache')

# EaCoNiC.SEG.aggregate(SEG.rds.list = list.files(path = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/CNA', pattern = '*.SEG.*.RDS$', full.names = TRUE, recursive = TRUE), out.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev', L2R = TRUE, BAF = TRUE, cohort.name = 'TEST2', return.data = FALSE)

## Annotate genomic regions
EaCoN.regions.annotate <- function(regions.df = NULL, genome = 'hg19') {
  
  ## Get genome cache
  # EaCoN::get.cached.genome(genome)
  get.cached.genome(genome)
  
  ## Cytobands
  if (!is.null(annot.obj$cytoBandIdeo)) {
    message('Annotating regions with cytologic bands ...')
    regions.df$Cyto.End <- regions.df$Cyto.Start <- NA
    ## Start
    fin.cstart.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = 'Chr', start.field = 'Start', end.field = 'Start', keep.extra.columns = FALSE)
    cyto.start.ol <- GenomicRanges::findOverlaps(query = fin.cstart.gr, subject = annot.obj$cytoBandIdeo)
    regions.df$Cyto.Start[S4Vectors::from(cyto.start.ol)] <- annot.obj$cytoBandIdeo$name[S4Vectors::to(cyto.start.ol)]
    rm(fin.cstart.gr, cyto.start.ol)
    ## End
    fin.cend.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, seqnames.field = 'Chr', start.field = 'End', end.field = 'End', keep.extra.columns = FALSE)
    cyto.end.ol <- GenomicRanges::findOverlaps(query = fin.cend.gr, subject = annot.obj$cytoBandIdeo)
    regions.df$Cyto.End[S4Vectors::from(cyto.end.ol)] <- annot.obj$cytoBandIdeo$name[S4Vectors::to(cyto.end.ol)]
    rm(fin.cend.gr, cyto.end.ol)
  }
  
  ## DGV
  if (!is.null(annot.obj$dgvMerged)) {
    message('Annotating regions with DGV CNV ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, keep.extra.columns = TRUE)
    regions.df$DGV.CNV <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$dgvMerged)
    rm(fin.gr)
  }
  
  ## GenomicSuperDups
  if (!is.null(annot.obj$genomicSuperDups)) {
    message('Annotating regions with UCSC genomic superduplications ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, keep.extra.columns = TRUE)
    regions.df$GenomicSuperDups <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$genomicSuperDups)
    rm(fin.gr)
  }
  
  ## CpGislands
  if (!is.null(annot.obj$cpgIslandExt)) {
    message('Annotating regions with UCSC CpG islands ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, keep.extra.columns = TRUE)
    regions.df$CpGislands <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$cpgIslandExt)
    rm(fin.gr)
  }
  
  ## Refseq Genes
  if (!is.null(annot.obj$refGene)) {
    message('Annotating regions with RefSeq genes ...')
    fin.gr <- GenomicRanges::makeGRangesFromDataFrame(df = regions.df, keep.extra.columns = TRUE)
    regions.df$RefSeq.Genes <- GenomicRanges::countOverlaps(query = fin.gr, subject = annot.obj$refGene)
    genes.oldf <- as.data.frame(GenomicRanges::findOverlaps(query = fin.gr, subject = annot.obj$refGene))
    regions.df$Symbols <- NA
    unique.hits <- unique(genes.oldf$queryHits)
    genes.list <- sapply(unique.hits, function(x) {
      paste(annot.obj$refGene$name2[genes.oldf$subjectHits[genes.oldf$queryHits == x]], collapse = ',')
    }, simplify = TRUE)
    regions.df$Symbols[unique.hits] <- genes.list
    rm(fin.gr, genes.oldf, unique.hits, genes.list)
  }
  return(regions.df)
}


## Perform differential analysis on segmentation regions (... correspond to additional parameters to read.table, to parse the annotation table)
## SEG.cohort.obj     Robject                 The output from the EaCoNiC.SEG.aggregate() function
## factors.file         character             A path to a table that contains at least two columns : one with samplenames, the other(s) being factors to use to perform the diff tests
## samples.colnames   character               Column name of factors.file to identify samples
## factors.name       vector of characters    Column name(s) of factors.file to specify classes to compare
## test.2             character               Type of test to use for factors with 2 classes. 'W' for Wilcoxon, 'T' for Student's t-test.
## test.mult          character               Type of test to use for factors with >2 classes. 'ANOVA' is self-explanatory, 'KW' for Kruskal-Wallis.
## data.type          character               Type of data to use from the SEG.cohort.obj. 'L2R.cut' for log2ratio values with normal segments set to 0, 'L2R.nocut' for log2ratio without normal filtering, 'BAF' for B-allelic frequency values.
## annotate           logical                 If TRUE, and if such data exist in the genome cache, use it to annotate the regions content (Refseq genes for most genomes, DGV and more for homo sapiens)
## out.dir            character               Path to a folder where results will be written
## ...                ...                     Any additional parameter to read.table() to properly parse the factors file.
EaCoNiC.DiffReg <- function(SEG.cohort.obj = NULL, factors.file = NULL, samples.colname = 'Sample', factors.name = NULL, test.2 = 'W', test.mult = 'KW', data.type = 'L2R.cut', annotate = TRUE, out.dir = getwd(), ...) {
  ## Basic checks
  if(is.null(SEG.cohort.obj)) stop('A SEG cohort object is required !')
  if(is.null(factors.file)) stop('An table with factors to compare is required !')
  if(!file.exists(factors.file)) stop('Could not find the provided factors file !')
  if(length(factors.name) == 0) stop('At least one factor to analyze is required !')
  if(!is.character(samples.colname)) stop('A column name to read sample names from the factors table is required !')
  if (!toupper(test.2) %in% c('W', 'T')) stop("Test types to analyze 2-classes factors are 'W' (wilcoxon) or 'T' (t-test).")
  if (!toupper(test.mult) %in% c('KW', 'ANOVA')) stop("Test types to analyze factors with more than 2 classes are 'KW' (Kruskal-Wallis) or 'ANOVA' (self-explanatory).")
  avail.datatypes <- c('L2R.cut', 'L2R.nocut', 'BAF')
  if (!data.type %in% avail.datatypes) stop(paste0("Data type to use should be one of : '", paste(avail.datatypes, collapse = "', '"), "' !"))
  
  ## Loading the factors file
  factors.df <- read.table(file = factors.file, ...)
  
  ## Secondary checks
  if(!samples.colname %in% colnames(factors.df)) stop(paste0("Could not find the samples column name '", samples.colname, "' in the provided annotation file !"))
  sn.in.annot <- SEG.cohort.obj$samples.info$sample.name %in% factors.df[[samples.colname]]
  if(!any(sn.in.annot)) stop("The annotation table does not contain any sample name in common with the SEG cohort object !")
  fac.in.annot <- factors.name %in% colnames(factors.df)
  if (!any(fac.in.annot)) stop("The annotation table does not contain any of the provided differential factors !")
  if(!all(fac.in.annot)) {
    message(paste0("The annotation table does not contain some of the provided differential factors : '", paste(factors.name[!fac.in.annot], collapse = "', '"), "' !"))
    factors.name <- factors.name[fac.in.annot]
  }
  
  main.dir <- paste0(out.dir, '/Differential_Analysis/')
  ## Looping through factors
  for (my.fac in factors.name) {
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
        # final.df <- EaCoN::EaCoN.region.annotate(my.df = final.df, genome = SEG.cohort.obj$meta$genome)
        final.df <- EaCoN.regions.annotate(regions.df = final.df, genome = SEG.cohort.obj$meta$genome)
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
      mini.seg <- cur.seg[,c(1:3, colnames(cur.seg) %in% mini.annot[[samples.colname]])]
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
      
      final.df <- cbind(cur.seg[,c(1:3)], Width = cur.seg$End - cur.seg$Start +1, mini.med, DIFF = mini.med[,1] - mini.med[,2], mini.sd, out.df)
      colnames(final.df)[colnames(final.df) == 'DIFF'] <- paste0(c(my.comb, data.type, 'Med.Diff'), collapse ='.')
      
      ## Performing annotation
      if (annotate) {
        # final.df <- EaCoN::EaCoN.region.annotate(my.df = final.df, genome = SEG.cohort.obj$meta$genome)
        final.df <- EaCoN.regions.annotate(regions.df = final.df, genome = SEG.cohort.obj$meta$genome)
      }
      
      ## Writing results
      diff.dir <- paste0(c(out.dir, my.fac, my.comb), collapse = '/')
      dir.create(diff.dir, recursive = TRUE)
      write.table(final.df, file = paste0(diff.dir, '/', paste(c(my.fac, my.comb, data.type, test.2), collapse = '_'), '.txt'), sep = "\t", quote = FALSE, row.names = FALSE)
      
      ## Volcano plot
      pcuts <- c(.05, .01)
      png(filename = paste0(diff.dir, '/', paste(c(my.fac, my.comb, data.type, test.2, 'volcano.png'), collapse = '_')), width = 1024, height = 768)
      plot(final.df[,7], -log10(final.df$RawP), xlim = max(abs(final.df[,7]))*c(-1,1), pch = 20, col = densCols(final.df[,7], -log10(final.df$RawP), bandwidth = .05), main = sub(pattern = '_vs_', replacement = ' vs ', x = paste0(my.fac, ' : ', my.comb, "\nVolcano plot")), xlab = 'L2R Difference', ylab = '-log10(RawP)')
      abline(v = 0, col = 4, lwd = 3, lty = 2)
      abline(h = -log10(pcuts), col = c('brown', 2), lwd = 3, lty = 3)
      text(max(abs(final.df[,7]))*.95, (-log10(pcuts)), labels = c(length(which(final.df$RawP < pcuts[1])), length(which(final.df$RawP < pcuts[2]))), col = c('brown', 2), pos = 3)
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
  test.function <- base::get(func.split[2], envir = loadNamespace(func.split[1]))
  
  ## Performing the test
  out.df <- data.frame(Stat = rep(NA, nrow(my.mat)))
  out.df$RawP <- out.df$Stat
  for (k in 1:nrow(my.mat)) {
    test.res <- test.function(as.numeric(my.mat[k,]) ~ my.factor)
    if (test.type == 'ANOVA') {
      test.res <- stats::anova(test.res)
    }
    out.df$Stat[k] <- test.res[[stat.slotname]][1]
    out.df$RawP[k] <- test.res[[pval.slotname]][1]
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
## factors.file         character (or NULL)             A path to a table that contains at least two columns : one with samplenames, the other(s) being factors to use to perform the diff tests
## samples.colnames     character (or NULL)             Column name of factors.file to identify samples
## factors.name         vector of characters (or NULL)  Column name(s) of factors.file to specify classes to compare
## chr.simplify         logical                         Simplify the chromosome names for plots by removing a potential common root (by example, 'chr' for homo sapiens).
## target.bins          numeric (or NULL)               To simulate the effect of the size of the region in the distance computation, larger regions are repeated multiple times (bins). The value of target.bins is the total number of bins : the function will automatically multiply the regions to 
EaCoNiC.heatmap <- function(SEG.cohort.obj = NULL, data.type = 'L2R.cut', chr.discard = NULL, dist.method = 'pearson', agglo.method = "ward.D", amp.cut = 1.5, del.cut = -1.5, factors.file = NULL, samples.colname = 'SampleName', factors.name = NULL, chr.simplify = TRUE, target.bins = NULL, out.dir = getwd(), ...) {
  ## Basic checks
  ### SEG object
  if(is.null(SEG.cohort.obj)) stop('A SEG cohort object is required !')
  ### Factors file
  if(is.null(factors.file)) message('No factor table provided.') else {
    if(!file.exists(factors.file)) stop('Could not find the provided factors file !')
    if(length(factors.name) == 0) stop('A factor table was provided, but no factor to plot !')
    if(!is.character(samples.colname)) stop('A factor table was provided, but nocolumn name identify sample names !')
  }
  ### Data type
  avail.datatypes <- c('L2R.cut', 'L2R.nocut', 'BAF')
  if(!data.type %in% avail.datatypes) stop(paste0("Data type to use should be one of : '", paste(avail.datatypes, collapse = "', '"), "' !"))
  ### Chromosomes to discard
  avail.chr <- levels(SEG.cohort.obj[[data.type]]$Chr)
  if(!is.null(chr.discard)) if (!all(chr.discard %in% avail.chr)) stop('At least one of requested chromosomes to discard is not available for this cohort !')
  ## Output dir
  if(!dir.exists(out.dir)) stop('The specified output directory does not exist !')
  ## Target bins
  if(!is.null(target.bins) & !is.numeric(target.bins)) stop('The target number of width bins should be an integer !')
  if(target.bins < 1) stop('The target number of width bins should be a positive integer !')
  
  ## Loading the factors file
  if (!is.null(factors.file)) {
    factors.df <- read.table(file = factors.file, ...)
    ## Factors-specific checks
    if(!samples.colname %in% colnames(factors.df)) stop(paste0("Could not find the samples column name '", samples.colname, "' in the provided annotation file !"))
    sn.in.annot <- SEG.cohort.obj$samples.info$sample.name %in% factors.df[[samples.colname]]
    if(!any(sn.in.annot)) stop("The annotation table does not contain any sample name in common with the SEG cohort object !")
    fac.in.annot <- factors.name %in% colnames(factors.df)
    if (!any(fac.in.annot)) stop("The annotation table does not contain any of the provided differential factors !")
    if(!all(fac.in.annot)) {
      message(paste0("The annotation table does not contain some of the provided differential factors : '", paste(factors.name[!fac.in.annot], collapse = "', '"), "' !"))
      factors.name <- factors.name[fac.in.annot]
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
  }
  
  ## Converting df to a matrix
  my.mat <- as.matrix(my.df[,-c(1:3)])
  
  ## Clustering
  message('Performing the hierarchical clustering ...')
  my.dist <- amap::Dist(x = t(my.mat), method = dist.method, nbproc = 1)
  ### Dealing with methods that dislike missing values
  if(any(is.na(my.dist))) {
    my.mat.tmp <- my.mat
    ### Replacing NA values by row median
    for(x in seq_len(nrow(my.mat))) my.mat.tmp[x,is.na(my.mat[x,])] <- median(my.mat[x,], na.rm = TRUE)
    my.dist <- amap::Dist(x = t(my.mat.tmp), method = dist.method, nbproc = 1)
  }
  my.hc <- stats::hclust(d = my.dist, method = agglo.method)
  
  ## Preparing graphical split on chromosomes
  chr.split <- as.factor(my.df$Chr)
  if (chr.simplify) levels(chr.split) <- chr.simplifier(levels(chr.split))
  
  ## Palette
  zlim = c(-1, 0, 1)
  myPalette <- c("red", "grey85", "blue")
  myRamp <- circlize::colorRamp2(zlim, myPalette)
  
  ## Output dir
  hch.dir <- paste0(out.dir, '/Clusterized_Heatmap/', data.type, '/', dist.method, '_', agglo.method, if(!is.null(chr.discard)) paste0('/NO.', paste(chr.discard, collapse = '.'), if(is.null(target.bins)) '' else '_ww', '/') else '/ALL_chromosomes/')
  dir.create(hch.dir, recursive = TRUE)
  
  ## Heatmap (portrait)
  message('Plotting the portrait heatmap ...')
  # Creating sample annotation
  hap <- if (!is.null(factors.file)) ComplexHeatmap::HeatmapAnnotation(df = factors.df[,factors.name]) else NULL
  # Plotting
  hmp <- suppressMessages(ComplexHeatmap::Heatmap(my.mat, name = data.type,
                                                  show_row_names = FALSE,
                                                  col = myRamp,
                                                  column_title = SEG.cohort.obj$meta$cohort.name,
                                                  cluster_columns = my.hc,
                                                  cluster_rows = FALSE,
                                                  top_annotation = hap,
                                                  cell_fun = function(j, i, x, y, width, height, fill) {
                                                    if(my.mat >= amp.cut) grid.circle(x = x, y = y, r = min(unit.c(width, height)) * 2, gp = gpar(fill = 'cyan4', col = NA))
                                                    if(my.mat <= del.cut) grid.circle(x = x, y = y, r = min(unit.c(width, height)) * 2, gp = gpar(fill = 'orangered4', col = NA))
                                                  },
                                                  row_split = chr.split))
  ## Writing
  png(filename = paste0(hch.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, data.type, 'heatmap', 'P'), collapse = '_'), '.png'),
      width = ncol(my.mat)*35 + if(!is.null(factors.name)) 200 else 0,
      # height = nrow(my.mat)/3 + 300
      height = nreg.ori/3 + 300
  )
  ComplexHeatmap::draw(hmp)
  dev.off()
  
  ## Heatmap (landscape)
  message('Plotting the landscape heatmap ...')
  hal <- if (!is.null(factors.file)) ComplexHeatmap::rowAnnotation(df = factors.df[,factors.name], show_annotation_name = FALSE) else NULL
  hml <- suppressMessages(ComplexHeatmap::Heatmap(t(my.mat), name = data.type,
                                                  show_column_names = FALSE,
                                                  col = myRamp,
                                                  cluster_rows = my.hc,
                                                  cluster_columns = FALSE,
                                                  row_title = SEG.cohort.obj$meta$cohort.name,
                                                  left_annotation = hal,
                                                  column_split = chr.split))
  ## Writing
  png(filename = paste0(hch.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, data.type, 'heatmap', 'L'), collapse = '_'), '.png'),
      # width = nrow(my.mat)/3 + if(!is.null(factors.name)) 300 else 0,
      width = nreg.ori/3 + if(!is.null(factors.name)) 300 else 0,
      height = ncol(my.mat)*35)
  ComplexHeatmap::draw(hml)
  dev.off()
  
  ## HC split membership
  hc.split <- cutree(my.hc, k = 2:ncol(my.mat))
  colnames(hc.split) <- paste0('HC.', colnames(hc.split))
  write.table(data.frame(Sample = rownames(hc.split), hc.split), file = paste0(hch.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, data.type, 'membership'), collapse = '_'), '.txt'), quote = FALSE, row.names = FALSE, sep = "\t")
}


## Auto-correlation map
EaCoNiC.autocor <- function(SEG.cohort.obj = NULL, data.type = 'L2R.cut', cor.method = 'spearman', chr.discard = NULL, chr.simplify = TRUE, out.dir = getwd()) {
  ## Basic checks
  if(is.null(SEG.cohort.obj)) stop('A SEG cohort object is required !')
  ### Data type
  avail.datatypes <- c('L2R.cut', 'L2R.nocut', 'BAF')
  if(!data.type %in% avail.datatypes) stop(paste0("Data type to use should be one of : '", paste(avail.datatypes, collapse = "', '"), "' !"))
  ### Chromosomes to discard
  avail.chr <- levels(SEG.cohort.obj[[data.type]]$Chr)
  if(!is.null(chr.discard)) if (!all(chr.discard %in% avail.chr)) stop('At least one of requested chromosomes to discard is not available for this cohort !')
  ## Output dir
  if(!dir.exists(out.dir)) stop('The specified output directory does not exist !')
  
  ## Getting data
  my.df <- SEG.cohort.obj[[data.type]]
  
  ## Limiting chromosomes
  if(!is.null(chr.discard)) {
    message(paste0('Discarding chromosome(s) [', paste(chr.discard, collapse = ', '), '] ...'))
    my.df <- my.df[!my.df$Chr %in% chr.discard,]
  }
  
  ## Converting to matrix
  my.mat <- as.matrix(my.df[,-c(1:3)])
  
  ## Computing correlation matrix
  message('Computing region correlations ...')
  my.cor <- stats::cor(t(my.mat), method = cor.method)
  
  ## Handling NA values
  if(any(is.na(my.cor))) {
    my.mat.tmp <- my.mat
    for(x in seq_len(nrow(my.mat))) my.mat.tmp[x,is.na(my.mat[x,])] <- median(my.mat[x,], na.rm = TRUE)
    my.cor <- stats::cor(t(my.mat.tmp), method = cor.method)
    rm(my.mat.tmp)
  }
  
  ## Correlation matrix segmentation
  message('Segmenting the correlation map ...')
  max.blocks <- min(floor(nrow(my.cor)/20), 500)
  bs.res <- blockseg::blockSeg(Y = my.cor, max.break = max.blocks, verbose = FALSE)
  ## Computing median correlation scores per block
  ### Getting actual breaks
  my.breaks <- sort(unique(c(bs.res@RowBreaks[[length(bs.res@RowBreaks)]], bs.res@ColBreaks[[length(bs.res@ColBreaks)]], nrow(my.cor))))
  ##Filtering contiguous breaks
  message('Filtering contiguous breakpoints ...')
  runs <- unname(split(seq_along(my.breaks), cumsum(c(0, diff(my.breaks) > 1))))
  streams.idx <- unlist(lapply(runs[lengths(runs) > 1], range))
  unik.idx <- unlist(runs[which(vapply(runs, length, 1L) ==1)])
  my.breaks2 <- sort(unique(my.breaks[sort(c(streams.idx, unik.idx))], nrow(my.cor)))
  
  # png('test.png', width = 2000, height = 2000)
  # image(x = seq_len(ncol(my.cor)), y = seq_len(nrow(my.cor)), z = t(my.cor))
  # abline(v = my.breaks2, col = 4, lwd = 2)
  # abline(h = my.breaks2, col = 4, lwd = 2)
  # dev.off()
  # rotate90 <- function(mat) t(mat[nrow(mat):1,,drop=FALSE])
  # rotate180 <- function(x) { x[] <- rev(x); return(x) }
  
  ## Preparing graphical split on chromosomes
  chr.split <- as.factor(my.df$Chr)
  if (chr.simplify) levels(chr.split) <- chr.simplifier(levels(chr.split))
  
  ## Preparing graphical split correlation map blocks
  mb.diff <- diff(my.breaks2)
  mb.diff[1] <- mb.diff[1]+1
  blocks.split <- as.factor(rep(seq_along(mb.diff), times = mb.diff))
  
  ## Generating output dir
  cor.dir <- paste0(out.dir, '/Auto_Correlation/', data.type, '/')
  dir.create(cor.dir)
  ## Saving the correlation map
  saveRDS(my.cor, file = paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'cormap'), collapse = '_'), '.RDS'), compress = 'bzip2')
  
  ## Palette
  zlim = c(-1, 0, 1)
  myPalette <- c("red", "grey85", "blue")
  myRamp <- circlize::colorRamp2(zlim, myPalette)
  
  ## Plotting the correlation map
  hac <- ComplexHeatmap::HeatmapAnnotation(df = chr.split, which = "row", name = 'Chr')
  
  hmp <- suppressMessages(ComplexHeatmap::Heatmap(my.cor, name = cor.method,
                                                  show_row_names = FALSE,
                                                  show_column_names = FALSE,
                                                  col = myRamp,
                                                  cluster_rows = FALSE,
                                                  cluster_columns = FALSE,
                                                  row_split = blocks.split,
                                                  column_split = blocks.split,
                                                  row_gap = grid::unit(.3, 'mm'),
                                                  column_gap = grid::unit(.3, 'mm'),
                                                  right_annotation = hac))
  
  png(paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'cormap_blocks'), collapse = '_'), '.png'), width = 2000 + (length(my.breaks2)), height = 2000)
  ComplexHeatmap::draw(hmp)
  dev.off()
  
  
  ## Saving results
  saveRDS(bs.res, file = paste0(cor.dir, '/', paste(c(SEG.cohort.obj$meta$cohort.name, cor.method, 'corseg'), collapse = '_'), '.RDS'), compress = 'bzip2')
  dev.off()
}


## Subset
EaCoNiC.SEG.subset <- function(SEG.cohort.obj = NULL, samples.keep = NULL, chr.keep = NULL, cohort.name = NULL, out.dir = getwd(), return.data = FALSE){
  
}


## test
EaCoN.set.cache(cache.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/EaCoN_cache')

# EaCoNiC.DiffReg(SEG.cohort.obj = readRDS('/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/TEST2.SEG.COHORT.RDS'), factors.file = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/samples_info.txt', samples.colname = 'sn', factors.name = c('Disease_control_benefit', 'BOR'), test.2 = 'W', test.mult = 'KW', data.type = 'L2R.cut', annotate = TRUE, out.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev', header = TRUE, sep = "\t")

# EaCoNiC.heatmap(SEG.cohort.obj = readRDS('/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/B20067_FRDA_05.SEG.COHORT.RDS'), data.type = 'L2R.cut', factors.file = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/samples_info.txt', samples.colname = 'sn', factors.name = c('Disease_control_benefit', 'BOR'), dist.method = 'pearson', agglo.method = 'ward.D', amp.cut = 1.5, loss.cut = -1.5, out.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev', header = TRUE, sep = "\t", chr.discard = c('chrX', 'chrY'), target.bins = 1E+05, chr.simplify = TRUE)

# EaCoNiC.autocor <- function(SEG.cohort.obj = readRDS('/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev/B20067_FRDA_05.SEG.COHORT.RDS'), data.type = 'L2R.cut', chr.discard = NULL, cor.method = 'spearman', chr.simplify = TRUE, out.dir = '/home/job/WORKSPACE/B20067_FRDA_05/ANALYSIS/P31_AUMA/dev')
