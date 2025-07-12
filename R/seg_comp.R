
## Factorization of the process

## FUNCTION to load FACET's VCF
facets_vcf_reader <- function(vcf_file = NULL, seg_colname = 'CNLR_MEDIAN_CLUST') {
  if (!file.exists(vcf_file)) stop("Can't find vcf_file [", vcf_file, ']')
  vdata <- bedr::read.vcf(x = vcf_file, split.info = TRUE)
  return(list(seg = data.frame(chr = vdata$vcf$CHROM, start = vdata$vcf$POS, end = vdata$vcf$POS + vdata$vcf$SVLEN, value = vdata$vcf[[seg_colname]]), meta = list(ploidy = as.numeric(vdata$header$ploidy), purity = as.numeric(vdata$header$purity))))
}

## FUNCTION ::Create a chromosomes-like object from a loaded BSgenome object
chromobjector <- function(BSg = NULL, chrs = NULL) {
  if (is.null(BSg)) stop("NULL object !", call. = FALSE)
  chromobj <- list(species = GenomeInfoDb::organism(BSg), genomebuild = metadata(BSg)$genome)
  chromdf <- data.frame(chrom = BSgenome::seqnames(BSg), chrN = seq_along(BSgenome::seqnames(BSg)), chr.length = GenomeInfoDb::seqlengths(BSg), stringsAsFactors = FALSE)
  
  ## Restrict ?
  ### Check
  if (!is.null(chrs)) {
    if (!all(chrs %in% chromdf$chrom)) stop('"chrs" provided but not all given values found in the provided BSgenome !') else chromdf <- chromdf[chromdf$chrom %in% chrs,]
  }
  
  chromdf$chr.length.sum <- cumsum(as.numeric(chromdf$chr.length))
  chromdf$chr.length.toadd <- c(0, chromdf$chr.length.sum[-nrow(chromdf)])
  chromdf$mid.chr <- round(diff(c(0, chromdf$chr.length.sum)) /2)
  chromdf$mid.chr.geno <- chromdf$mid.chr + chromdf$chr.length.toadd
  chromobj$chromosomes <- chromdf
  rm(chromdf)
  chromobj$chrom2chr <- sapply(chromobj$chromosomes$chrom, function(k) { chromobj$chromosomes$chrN[chromobj$chromosomes$chrom == k]}, simplify = FALSE)
  chromobj$chr2chrom <- sapply(chromobj$chromosomes$chrN, function(k) { chromobj$chromosomes$chrom[chromobj$chromosomes$chrN == k]}, simplify = FALSE)
  names(chromobj$chr2chrom) <- chromobj$chromosomes$chrN
  chromobj$genome.length <- sum(as.numeric(chromobj$chromosomes$chr.length), na.rm = TRUE)
  return(chromobj)
}

## FUNCTION to perform the segments comparison.
##  Can optionally take a bed file (for WES/targeted data) to visualize the areas of coverage.
##  Can optionally take a UCSC cytobandideo table ty visualize centromeres.
##  Can optionally use the profiles' ploidy to adapt Y plot limits
##  Requires an installed BSgenome package name to build chr data (for genomic coordinates)
##  2 modes of action : 'L2R' or 'TCN' (adapts Y plot limits and types of correlation computed)
##  For Pearson correlations, abs(L2R) are capped to avoid artificial increase of the correlation score.
##  FACETS segmentation profiles are not well centered, so an option to roughly recenter it is proposed. Only used in L2R mode, it will subtract the median of segments within the [-0.5:0.5] range.
##  For very similar values in the correlated profiles (especially for TCN), a random value is added, from a random distribution of median 0 and epsilon as its standard deviation.
seg_comp <- function(prof1_df = NULL, prof2_df = NULL, prof1_ploidy = NULL, prof2_ploidy = NULL, prof1_name = 'Profile1', prof2_name = 'Profile2', diff_max = 0.1, genome_package = 'BSgenome.Hsapiens.UCSC.hg38', prof_mod = 'L2R', prof1_l2r_recenter = FALSE, prof2_l2r_recenter = FALSE, capping = 1.5, epsilon = 1E-16, plot = TRUE, prof1_col = 'blue', prof2_col = 'red', cytobandideo_file = NULL, capture_bed_file = NULL, out_dir = getwd()) {
  
  valid_modes <- c('L2R', 'TCN')
  prof_mod <- toupper(prof_mod)
  
  ## Checks
  if(is.null(prof1_df)) stop('A segmentation table for the 1st profile is required !')
  if(is.null(prof2_df)) stop('A segmentation table for the 2nd profile is required !')
  if(!is.null(prof1_ploidy) & !is.numeric(prof1_ploidy)) stop('Ploidy for the 1st profile should be numeric !')
  if(!is.null(prof2_ploidy) & !is.numeric(prof2_ploidy)) stop('Ploidy for the 1st profile should be numeric !')
  if(!is.null(cytobandideo_file) & !file.exists(cytobandideo_file)) stop('Could not find the UCSC cytobandideo table ["', cytobandideo_file, '"] !')
  if(!is.null(capture_bed_file) & !file.exists(capture_bed_file)) stop('Could not find the capture bed ["', capture_bed_file, '"] !')
  if (!prof_mod %in% valid_modes) stop('Mode ["', prof_mod, '"] not supported ! Please use any of ["', paste(valid_modes, collapse = '", "'), '"]')
  if(!dir.exists(out_dir)) stop('Output directory ["', out_dir, '"] does not exist !')
  if (!'BSgenome' %in% utils::installed.packages()) stop('An installed BSgenome is required. For this, the "BSgenome" CRAN package should be installed !')
  if (!genome_package %in% BSgenome::installed.genomes()) stop('The ["', genome_package, '"] BSgenome package is not installed !')
  
  if (prof_mod == 'TCN') {
    if (any(c(prof1_l2r_recenter, prof1_l2r_recenter))) {
      message('Log2ratio recentering requested but TCN mode active : no recentring, as not needed.')
      prof2_l2r_recenter <- prof1_l2r_recenter <- FALSE
    }
    if (!is.null(capping)) {
      if (round(capping) != capping) stop('TCN mode active, and provided capping value is not an integer ! Are you sure you did not provided a capping value compatible for L2R mode ?')
    }
    # if (all(is.null(c(prof1_ploidy, prof2_ploidy)))) message("TCN mode active, but no ploidy given : automatic Y-axis scaling won't be performed, and segment blocks won't get drawn...")
  }
  
  ## RUN
  
  ## Create chromobjector object
  message('Loading [', genome_package, '] genome info ...')
  suppressPackageStartupMessages(require(genome_package, character.only = TRUE))
  BSg.obj <- getExportedValue(genome_package, genome_package)
  ### Restrict to available chromosomes
  avail_chrs <- unique(c(prof1_df$chr, prof2_df$chr))
  ### Create cs object
  cs <- chromobjector(BSg = BSg.obj, chrs = avail_chrs)

  ## Load cytobandideo (for plots)
  if (!is.null(cytobandideo_file)) {
    message('Loading cytobandideo table ...')
    cbi <- data.table::fread(file = cytobandideo_file, data.table = FALSE)
    colnames(cbi)[1] <- 'chrom'
    ### Restrict
    cbi <- cbi[cbi$chrom %in% avail_chrs,]
    ### Get centromere pos
    cpos_tbl <- tibble::as_tibble(x = cbi)
    cpos_tbl <- cpos_tbl[cpos_tbl$gieStain %in% c('acen', 'gvar', 'stalk'),]
    cpos_tbl <- dplyr::group_by(.data = cpos_tbl, chrom, gieStain)
    suppressMessages(cpos_tbl <- dplyr::summarize(.data = cpos_tbl, start = min(chromStart), end = max(chromEnd)))
    cpos_tbl$Gstart <- cpos_tbl$start + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[cpos_tbl$chrom])]
    cpos_tbl$Gend <- cpos_tbl$end + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[cpos_tbl$chrom])]
    rm(cbi)
  }
  
  ## Load capture bed (for plots)
  if (!is.null(capture_bed_file)) {
    message('Loading capture bed data ...')
    bedgr <- data.table::fread(file = capture_bed_file, sep = '\t', data.table = FALSE)[,1:3]
    colnames(bedgr) <- c('chr', 'start', 'end')
    ### Restrict to available chromosomes
    bedgr <- bedgr[bedgr$chr %in% avail_chrs,]
    ## Add genomic coordinates
    bedgr$Gstart <- bedgr$start + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[bedgr$chr])]
    bedgr$Gend <- bedgr$end + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[bedgr$chr])]
  }
  
  ## TCN : if no ploidy given, try to infer it as the the most populated CN level
  if (prof_mod == 'TCN' & is.null(prof1_ploidy)) {
    message('TCN mode active but no ploidy given for first profile. Inferring it as the most populated TCN level ...')
    prof1_tbl <- tibble::as_tibble(x = prof1_df)
    prof1_tbl <- dplyr::group_by(.data = prof1_tbl, value)
    levlen <- dplyr::summarize(.data = prof1_tbl, width = sum(value))
    prof1_ploidy <- levlen$value[which.max(levlen$width)]
    rm(prof1_tbl, levlen)
  }
  if (prof_mod == 'TCN' & is.null(prof2_ploidy)) {
    message('TCN mode active but no ploidy given for second profile. Inferring it as the most populated TCN level ...')
    prof2_tbl <- tibble::as_tibble(x = prof2_df)
    prof2_tbl <- dplyr::group_by(.data = prof2_tbl, value)
    levlen <- dplyr::summarize(.data = prof2_tbl, width = sum(value))
    prof2_ploidy <- levlen$value[which.max(levlen$width)]
    rm(prof2_tbl, levlen)
  }
  
  ## Convert position tables to GRanges
  gen_info <- GenomeInfoDb::seqinfo(BSg.obj)
  ### Restricting to chromosomes available in provided profiles
  gen_info <- gen_info[unique(prof1_df$chr),]
  gen_info <- gen_info[unique(prof2_df$chr),]
  ### Conversion
  grg <- GenomicRanges::makeGRangesFromDataFrame(df = data.frame(chr = GenomeInfoDb::seqnames(gen_info), start = 1, end = GenomeInfoDb::seqlengths(gen_info)))
  gr1 <- GenomicRanges::makeGRangesFromDataFrame(df = prof1_df, keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = 'chr', start.field = 'start', end.field = 'end', seqinfo = gen_info)
  gr2 <- GenomicRanges::makeGRangesFromDataFrame(df = prof2_df, keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = 'chr', start.field = 'start', end.field = 'end', seqinfo = gen_info)
  
  ## Pool genome info and breakpoints
  inters <- IRanges::disjoin(unlist(GenomicRanges::GRangesList(list(grg, gr1, gr2))))
  
  ## Add genomic pos
  inters$Gstart <- inters@ranges@start + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[as.vector(inters@seqnames)])]
  inters$Gend <- inters@ranges@start + inters@ranges@width + cs$chromosomes$chr.length.toadd[unlist(cs$chrom2chr[as.vector(inters@seqnames)])]
  
  ## Overlaps
  ov1 <- GenomicRanges::findOverlaps(query = gr1, subject = inters, minoverlap = 1L)
  ov2 <- GenomicRanges::findOverlaps(query = gr2, subject = inters, minoverlap = 1L)
  
  ## Merge
  inters$PROF2 <- inters$PROF1 <- NA
  inters$PROF1[ov1@to] <- gr1$value[ov1@from]
  inters$PROF2[ov2@to] <- gr2$value[ov2@from]
  
  ## Convert to df
  inters <- as.data.frame(inters)
  colnames(inters)[1] <- 'chr'
  
  ## Recenter ?
  if (prof_mod == 'L2R' & prof1_l2r_recenter) {
    message('Recentering profile 1')
    # m1 <- median(inters$PROF1[inters$PROF1 > -.5 & inters$PROF1 < .5], na.rm = TRUE)
    # inters$PROF1 <- inters$PROF1 - p1_tab
    inters_tbl <- tibble::as_tibble(inters)
    p1_tbl <- dplyr::group_by(.data = inters_tbl, PROF1)
    levlen <- dplyr::summarize(.data = p1_tbl, width = sum(width))
    inters$PROF1 <- inters$PROF1 - levlen$PROF1[which.max(levlen$width)]
    rm(inters_tbl, p1_tbl, levlen)
  }
  if (prof_mod == 'L2R' & prof2_l2r_recenter) {
    message('Recentering profile 2')
    # m2 <- median(inters$PROF2[inters$PROF2 > -.5 & inters$PROF2 < .5], na.rm = TRUE)
    # inters$PROF2 <- inters$PROF2 - m2
    inters_tbl <- tibble::as_tibble(inters)
    p2_tbl <- dplyr::group_by(.data = inters_tbl, PROF2)
    levlen <- dplyr::summarize(.data = p2_tbl, width = sum(width))
    inters$PROF2 <- inters$PROF2 - levlen$PROF2[which.max(levlen$width)]
    rm(inters_tbl, p2_tbl, levlen)
  }
  
  ## Diff
  inters$Diff <- inters$PROF1 - inters$PROF2
  
  ## Create output dir
  sn <- paste(c(prof1_name, prof2_name), collapse = '_vs_')
  sn_dir <- paste0(out_dir, '/', paste(c(sn, prof_mod), collapse = '_'))
  dir.create(path = sn_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## Save comp table
  message('Saving the comparison table ...')
  out_df <- as.data.frame(inters)[, c(1:4,8:10)]
  colnames(out_df)[5:6] <- c(prof1_name, prof2_name)
  data.table::fwrite(x = out_df, file = paste0(sn_dir, '/', sn, '_comp.table.tsv'), sep = '\t', na = NA)
  rm(out_df)
  
  ## Compute correlations/metrics
  message('Computing the width-weighted correlations & metrics ...')
  cor_methods <- c('pearson', 'spearman', if(prof_mod == 'TCN') 'polychoric')
  metrics_list <- list()
  ### Genome
  nona <- !is.na(inters$PROF1) & !is.na(inters$PROF2)
  x2cor <- inters$PROF1[nona]
  y2cor <- inters$PROF2[nona]
  metrics_list$Genome <- vapply(cor_methods, function(cmi) {
    if (all(x2cor == y2cor)) return(1)
    ## Correlation (100 iterations to compensate for the epsilon effect)
    xwcres <- vapply(1:100, function(x) {
      if (sd(x2cor, na.rm = TRUE) == 0 & !is.null(epsilon)) x2cor <- x2cor + rnorm(n = length(x2cor), mean = 0, sd = epsilon)
      if (sd(y2cor, na.rm = TRUE) == 0 & !is.null(epsilon)) y2cor <- y2cor + rnorm(n = length(y2cor), mean = 0, sd = epsilon)
      # if (!is.null(epsilon)) {
      #   x2cor <- x2cor + rnorm(n = length(x2cor), mean = 0, sd = epsilon)
      #   y2cor <- y2cor + rnorm(n = length(y2cor), mean = 0, sd = epsilon)
      # }
      if (cmi == 'pearson') {
        ## Capping L2R for pearson method
        sat1L <- x2cor < -capping
        if(any(sat1L)) x2cor[sat1L] <- -capping
        sat1U <- x2cor > capping
        if(any(sat1U)) x2cor[sat1U] <- capping
        sat2L <- y2cor < -capping
        if(any(sat2L)) y2cor[sat2L] <- -capping
        sat2U <- y2cor > capping
        if(any(sat2U)) y2cor[sat2U] <- capping
      }
      xwc <- try(wCorr::weightedCorr(x = x2cor, y = y2cor, method = cmi, weights = inters$width[nona], ML = FALSE, fast = TRUE), silent = TRUE)
      if (inherits(x = xwc, 'try-error')) xwc <- NaN
      return(xwc)
    }, .1)
    wcres <- if (all(is.nan(xwcres))) NaN else median(x = xwcres, na.rm = TRUE)
  }, .1)
  names(metrics_list$Genome) <- paste0('ww', toupper(names(metrics_list$Genome)))
  
  ### Add identity fraction (including normal)
  nodiff <- abs(inters$Diff) < diff_max
  metrics_list$Genome <- c(metrics_list$Genome, setNames(object = sum(inters$width[c(which(nodiff), which(is.na(nodiff)))], na.rm = TRUE) / sum(inters$width), nm = 'Id_frac.noNA'), setNames(object = sum(inters$width[nodiff], na.rm = TRUE) / sum(inters$width), nm = 'Id_frac.wNA'))
  
  # ### Add identity fraction (excluding normal)
  # nodiff <- if (prof_mod == 'L2R') inters$Diff < diff_max & inters$PROF1 != 0 & inters$PROF2 != 0 else if (prof_mod == 'TCN') inters$Diff < diff_max & inters$PROF1 != prof1_ploidy & inters$PROF2 != prof2_ploidy
  # metrics_list$Genome <- c(metrics_list$Genome, setNames(object = sum(inters$width[c(which(nodiff), which(is.na(nodiff)))], na.rm = TRUE) / sum(inters$width), nm = 'Id_frac.NN.noNA'), setNames(object = sum(inters$width[nodiff], na.rm = TRUE) / sum(inters$width), nm = 'Id_frac.NN.wNA'))
  # 
  ## Chromosomes
  for (k in unique(inters$chr)) {
    nonak <- inters$chr == k & !is.na(inters$PROF1) & !is.na(inters$PROF2)
    x2cor <- inters$PROF1[nonak]
    y2cor <- inters$PROF2[nonak]
    metrics_list[[k]] <- vapply(cor_methods, function(cmi) {
      if (length(x2cor) == 1) return(NaN)
      if (all(x2cor == y2cor)) return(1)
      xwcres <- vapply(1:100, function(x) {
        if (sd(x2cor, na.rm = TRUE) == 0 & !is.null(epsilon)) x2cor <- x2cor + rnorm(n = length(x2cor), mean = 0, sd = epsilon)
        if (sd(y2cor, na.rm = TRUE) == 0 & !is.null(epsilon)) y2cor <- y2cor + rnorm(n = length(y2cor), mean = 0, sd = epsilon)
        # if (!is.null(epsilon)) {
        #   x2cor <- x2cor + rnorm(n = length(x2cor), mean = 0, sd = epsilon)
        #   y2cor <- y2cor + rnorm(n = length(y2cor), mean = 0, sd = epsilon)
        # }
        if (cmi == 'pearson') {
          ## Capping L2R for pearson method
          sat1L <- x2cor < -capping
          if(any(sat1L)) x2cor[sat1L] <- -capping
          sat1U <- x2cor > capping
          if(any(sat1U)) x2cor[sat1U] <- capping
          sat2L <- y2cor < -capping
          if(any(sat2L)) y2cor[sat2L] <- -capping
          sat2U <- y2cor > capping
          if(any(sat2U)) y2cor[sat2U] <- capping
        }
        xwc <- try(wCorr::weightedCorr(x = x2cor, y = y2cor, method = cmi, weights = inters$width[nonak], ML = FALSE, fast = TRUE), silent = TRUE)
        if (inherits(x = xwc, 'try-error')) xwc <- NaN
        return(xwc)
      }, .1)
      wcres <- if (all(is.nan(xwcres))) NaN else median(x = xwcres, na.rm = TRUE)
    }, .1)
    names(metrics_list[[k]]) <- paste0('ww', toupper(names(metrics_list[[k]])))
    ### Add identity fraction (including normal)
    # inters_k <- inters[inters$chr == k,]
    inters_k <- inters[nonak & inters$chr == k & (!is.na(inters$PROF1) | !is.na(inters$PROF2)),]
    nodiff <- abs(inters_k$Diff) < diff_max
    metrics_list[[k]] <- c(metrics_list[[k]], setNames(object = sum(inters_k$width[nodiff]) / sum(inters$width[nonak]), nm = 'Id_frac.noNA'), setNames(object = sum(inters_k$width[nodiff]) / sum(inters$width[inters$chr == k]), nm = 'Id_frac.wNA'))
    # ### Add identity fraction (excluding normal)
    # inters_k <- if (prof_mod == 'L2R') inters_k[!inters_k$PROF1 %in% 0 & !inters_k$PROF2 %in% 0,] else if (prof_mod == 'TCN') inters_k[!inters_k$PROF1 %in% prof1_ploidy & !inters_k$PROF2 %in% prof2_ploidy,]
    # nodiff <- inters_k$Diff < diff_max
    # metrics_list[[k]] <- c(metrics_list[[k]], setNames(object = sum(inters_k$width[nodiff]) / sum(inters$width[nonak]), nm = 'Id_frac.noNA'), setNames(object = sum(inters_k$width[nodiff]) / sum(inters$width[inters$chr == k]), nm = 'Id_frac.wNA'))
    
  }
  metrics_mat <- t(as.data.frame(metrics_list))
  metrics_df <- data.frame(Genomic_element = rownames(metrics_mat), metrics_mat)
  rm(metrics_list, metrics_mat)
  
  ## Save correlation table
  message('Saving the metrics ...')
  data.table::fwrite(x = metrics_df, file = paste0(sn_dir, '/', sn, '_metrics.table.tsv'), sep = '\t', na = NA)
  
  ## Plots ?
  if (plot) {
    ## Adjusting Y axis range
    y_range <- if (prof_mod == 'TCN') c(0,5) else if (prof_mod == 'L2R') c(-2,2) else c(0,0) 
    if (prof_mod == 'TCN' & !all(is.null(c(prof1_ploidy, prof2_ploidy)))) {
      y_range[2] <- (max(ceiling(c(prof1_ploidy, prof2_ploidy))) * 2) +1
    }
    
    if(!is.null(prof1_ploidy)) p1_ybasis <- round(prof1_ploidy)
    if(!is.null(prof2_ploidy)) p2_ybasis <- round(prof1_ploidy)
    
    ## Genomic view
    message("Generating pan-genomic plot ...")
    png(filename = paste0(sn_dir, '/', paste(c(sn, prof_mod, 'pangeno'), collapse = '_'), '.png'), width = 2160, height = 1080)
    par(mar = c(5,5,4,2))
    plot(x = 0, y = 0, ylim = y_range, xlim = c(inters$Gstart[1], rev(inters$Gend)[1]), type = 'n', xaxs = 'i', yaxs = 'i', xlab = 'Position', ylab = prof_mod, main = paste0(prof1_name, if(prof1_l2r_recenter) ' (recentered)', ' (', prof1_col, ') vs ', prof2_name, if(prof2_l2r_recenter) ' (recentered)', ' (', prof2_col, ') : segmented ', prof_mod, '\n', paste(vapply(colnames(metrics_df)[-1], function(cm) { paste0(cm, ' = ', round(metrics_df[[cm]][metrics_df$Genomic_element == 'Genome'], digits = 3)) }, 'a'), collapse = ' // ')), cex.axis = 2, cex.lab = 2, cex.main = 2)
    ## BED ?
    if (!is.null(capture_bed_file)) abline(v = bedgr$Gstart, col = 'grey95', lwd = 1)
    ## Medline
    if (prof_mod == 'L2R') abline(h = 0, lty = 3, lwd = 3)
    abline(v = cs$chromosomes$chr.length.sum, lty = 3, lwd = 3)
    ## Cytobandideo?
    if (!is.null(cytobandideo_file)) {
      acen_idx <- cpos_tbl$gieStain == 'acen'
      gvar_idx <- cpos_tbl$gieStain == 'gvar'
      stalk_idx <- cpos_tbl$gieStain == 'stalk'
      if (any(acen_idx)) rect(xleft = cpos_tbl$Gstart[acen_idx], ybottom = y_range[1], xright = cpos_tbl$Gend[acen_idx], ytop = y_range[2], border = NA, col = adjustcolor(col = 'yellow', alpha.f = .3))
      if (any(gvar_idx)) rect(xleft = cpos_tbl$Gstart[gvar_idx], ybottom = y_range[1], xright = cpos_tbl$Gend[gvar_idx], ytop = y_range[2], border = NA, col = adjustcolor(col = 'darkgreen', alpha.f = .3))
      if (any(stalk_idx)) rect(xleft = cpos_tbl$Gstart[stalk_idx], ybottom = y_range[1], xright = cpos_tbl$Gend[stalk_idx], ytop = y_range[2], border = NA, col = adjustcolor(col = 'brown', alpha.f = .3))
    }
    ## PROF1
    rect(xleft = inters$Gstart, ybottom = if (prof_mod == 'L2R') 0 else if(!is.null(prof1_ploidy)) p1_ybasis else inters$PROF1, xright = inters$Gend, ytop = inters$PROF1, border = adjustcolor(col = prof1_col, alpha.f = .3), col = adjustcolor(col = prof1_col, alpha.f = .3), lwd = 2)
    segments(x0 = inters$Gstart, y0 = inters$PROF1, x1 = inters$Gend, y1 = inters$PROF1, lwd = 15, col = prof1_col)
    ## PROF2
    rect(xleft = inters$Gstart, ybottom = if (prof_mod == 'L2R') 0 else if(!is.null(prof2_ploidy)) p2_ybasis else inters$PROF2, xright = inters$Gend, ytop = inters$PROF2, border = adjustcolor(col = prof2_col, alpha.f = .3), col = adjustcolor(col = prof2_col, alpha.f = .3), lwd = 2)
    segments(x0 = inters$Gstart, y0 = inters$PROF2, x1 = inters$Gend, y1 = inters$PROF2, lwd = 8, col = prof2_col)
    ## Chr names
    yrange = abs(diff(y_range))
    m.pos <- c(y_range[2] - (sign(y_range[2]) * yrange * .05), y_range[1] + (sign(y_range[2]) * yrange * .05))
    m.mod <- -(cs$chromosomes$chrN %% 2) + 2
    try(text(x = cs$chromosomes$mid.chr.geno, y = m.pos[m.mod], labels = cs$chromosomes$chrom, srt = 90, cex = 2))
    dev.off()
    
    ## Chromosomes
    k_dir <- paste0(sn_dir, '/chromosomes')
    dir.create(path = k_dir, recursive = TRUE, showWarnings = FALSE)
    
    message("Generating per-chromosome plots ...")
    for (k in unique(inters$chr)) {
      ### Plot k
      png(filename = paste0(k_dir, '/', paste(c(sn, prof_mod, k), collapse = '_'), '.png'), width = 2160, height = 1080)
      par(mar = c(5,5,4,2))
      plot(x = 0, y = 0, ylim = y_range, xlim = c(0, cs$chromosomes$chr.length[unlist(cs$chrom2chr[k])]), type = 'n', xaxs = 'i', yaxs = 'i', xlab = 'Position', ylab = prof_mod, main = paste0(prof1_name, if(prof1_l2r_recenter) ' (recentered)', ' (', prof1_col, ') vs ', prof2_name, if(prof2_l2r_recenter) ' (recentered)', ' (', prof2_col, ') : segmented ', prof_mod, ' // ', k, '\n', paste(vapply(colnames(metrics_df)[-1], function(cm) { paste0(cm, ' = ', round(metrics_df[[cm]][metrics_df$Genomic_element == k], digits = 3)) }, 'a'), collapse = ' // ')), cex.axis = 2, cex.lab = 2, cex.main = 2)
      ## BED ?
      if (!is.null(capture_bed_file)) abline(v = bedgr$start[bedgr$chr == k], col = 'grey95', lwd = 2)
      ## Medline
      if (prof_mod == 'L2R') abline(h = 0, lty = 3, lwd = 3)
      ## Cytobandideo?
      if (!is.null(cytobandideo_file)) {
        acen_idx <- cpos_tbl$chrom == k & cpos_tbl$gieStain == 'acen'
        gvar_idx <- cpos_tbl$chrom == k & cpos_tbl$gieStain == 'gvar'
        stalk_idx <- cpos_tbl$chrom == k & cpos_tbl$gieStain == 'stalk'
        if (any(acen_idx)) rect(xleft = cpos_tbl$start[acen_idx], ybottom = y_range[1], xright = cpos_tbl$end[acen_idx], ytop = y_range[2], border = 'black', lty = 3, lwd = 3, col = adjustcolor(col = 'yellow', alpha.f = .3))
        if (any(gvar_idx)) rect(xleft = cpos_tbl$start[gvar_idx], ybottom = y_range[1], xright = cpos_tbl$end[gvar_idx], ytop = y_range[2], border = 'black', lty = 3, lwd = 3, col = adjustcolor(col = 'darkgreen', alpha.f = .3))
        if (any(stalk_idx)) rect(xleft = cpos_tbl$start[stalk_idx], ybottom = y_range[1], xright = cpos_tbl$end[stalk_idx], ytop = y_range[2], border = 'black', lty = 3, lwd = 3, col = adjustcolor(col = 'brown', alpha.f = .3))
      }
      ## PROF1
      rect(xleft = inters$start[inters$chr == k], ybottom = if (prof_mod == 'L2R') 0 else if(!is.null(prof1_ploidy)) p1_ybasis else inters$PROF1[inters$chr == k], xright = inters$end[inters$chr == k], ytop = inters$PROF1[inters$chr == k], border = adjustcolor(col = prof1_col, alpha.f = .3), col = adjustcolor(col = prof1_col, alpha.f = .3), lwd = 2)
      segments(x0 = inters$start[inters$chr == k], y0 = inters$PROF1[inters$chr == k], x1 = inters$end[inters$chr == k], y1 = inters$PROF1[inters$chr == k], lwd = 15, col = prof1_col)
      ## PROF2
      rect(xleft = inters$start[inters$chr == k], ybottom = if (prof_mod == 'L2R') 0 else if(!is.null(prof2_ploidy)) p2_ybasis else inters$PROF2[inters$chr == k], xright = inters$end[inters$chr == k], ytop = inters$PROF2[inters$chr == k], border = adjustcolor(col = prof2_col, alpha.f = .3), col = adjustcolor(col = prof2_col, alpha.f = .3), lwd = 2)
      segments(x0 = inters$start[inters$chr == k], y0 = inters$PROF2[inters$chr == k], x1 = inters$end[inters$chr == k], y1 = inters$PROF2[inters$chr == k], lwd = 8, col = prof2_col)
      dev.off()
    }
  }
}







## DEMO

## FULL COHORT (using design file)

### COMMON
genome_package = 'BSgenome.Hsapiens.UCSC.hg38'
prof1_col <- 'blue'
prof2_col <- 'red'
cytobandideo_file <- '/home/job/WORKSPACE/ElsaLAB/CNA/EACON_TESTS/RESOURCES/cytoband.hg38'
capture_bed_file <- '/home/job/WORKSPACE/ElsaLAB/CNA/EACON_TESTS/RESOURCES/BEDs/twist_human_core_exome_and_integragen_custom_v2.bed'
out_dir <- "/home/job/WORKSPACE/ElsaLAB/CNA/seg_comp/RESULTS"
suppressPackageStartupMessages(require(genome_package, character.only = TRUE))
BSg.obj <- getExportedValue(genome_package, genome_package)


#########################
## FACETS vs FACETS (WES vs WGS)
#########################

### Common 2
techno1 <- 'WES'
techno2 <- 'WGS'
des_file <- '/home/job/WORKSPACE/ElsaLAB/CNA/seg_comp/FACETS.WES_vs_FACETS.WGS_design.tsv'

  ### Mode : L2R
  prof_mod = 'L2R'
  seg_colname <- 'CNLR_MEDIAN_CLUST'
  prof1_l2r_recenter <- TRUE
  prof2_l2r_recenter <- TRUE
  capping = 1.5

  ### Mode : TCN
  prof_mod = 'TCN'
  seg_colname <- 'TCN_EM'
  prof1_l2r_recenter <- FALSE
  prof2_l2r_recenter <- FALSE
  capping = 6

### RUN
des_df <- data.table::fread(file = des_file, sep = '\t', data.table = FALSE)
od2 <- paste(c(out_dir, prof_mod, paste0('FACETS.', techno1, '_vs_FACETS.', techno2)), collapse = '/')
dir.create(path = od2, recursive = TRUE)
### Looping on files
for (d in seq_len(nrow(des_df))) {
  message(des_df$prof1_name[d], ' vs ', des_df$prof2_name[d])
  p1 <- facets_vcf_reader(vcf_file = des_df$prof1_file[d], seg_colname = seg_colname)
  p2 <- facets_vcf_reader(vcf_file = des_df$prof2_file[d], seg_colname = seg_colname)
  prof1_ploidy <- if (prof_mod == 'L2R') NULL else p1$meta$ploidy
  prof2_ploidy <- if (prof_mod == 'L2R') NULL else p2$meta$ploidy
  seg_comp(prof1_df = p1$seg, prof2_df = p2$seg, prof1_ploidy = prof1_ploidy, prof2_ploidy = prof1_ploidy, prof1_name = des_df$prof1_name[d], prof2_name = des_df$prof2_name[d], genome_package = genome_package, prof_mod = prof_mod, prof1_l2r_recenter = prof1_l2r_recenter, prof2_l2r_recenter = prof2_l2r_recenter, capping = capping, epsilon = 1E-16, plot = TRUE, prof1_col = prof1_col, prof2_col = prof2_col, cytobandideo_file = cytobandideo_file, capture_bed_file = capture_bed_file, out_dir = od2)
}


#########################
## EACON.ASCAT.WES vs EACON.ASCAT.WGS
#########################

### Common 2
techno1 <- 'WES'
techno2 <- 'WGS'
prof1_l2r_recenter <- FALSE
prof2_l2r_recenter <- FALSE

  ### Mode : L2R
  prof_mod = 'L2R'
  capping = 1.5
  
  ### Mode : TCN
  prof_mod = 'TCN'
  capping = 6

## RUN
des_file <- paste0('/home/job/WORKSPACE/ElsaLAB/CNA/seg_comp/EACON.ASCAT.', techno1, '_vs_EACON.ASCAT.', techno2, '_design_', prof_mod, '.tsv')
des_df <- data.table::fread(file = des_file, sep = '\t', data.table = FALSE)
od2 <- paste(c(out_dir, prof_mod, paste0('EACON.ASCAT.', techno1, '_vs_EACON.ASCAT.', techno2)), collapse = '/')
dir.create(path = od2, recursive = TRUE)
for (d in seq_len(nrow(des_df))) {
  cs <- chromobjector(BSg = BSg.obj)
  message(des_df$prof1_name[d], ' vs ', des_df$prof2_name[d])
  ## Load EaCoN data 1
  p1obj <- readRDS(des_df$prof1_file[d])
  p1seg <- if (prof_mod == 'L2R') data.frame(chr = unlist(cs$chr2chrom[p1obj$cbs$cut$Chr]), start = p1obj$cbs$cut$Start, end = p1obj$cbs$cut$End, value = p1obj$cbs$cut$Log2Ratio) else data.frame(chr = p1obj$segments$chr, start = p1obj$segments$startpos, end = p1obj$segments$endpos, value = p1obj$segments$nMajor + p1obj$segments$nMinor)
  prof1_ploidy <- if (prof_mod == 'L2R') NULL else p1obj$ploidy$most.width
  ## Load EaCoN data 2
  p2obj <- readRDS(des_df$prof2_file[d])
  p2seg <- if (prof_mod == 'L2R') data.frame(chr = unlist(cs$chr2chrom[p2obj$cbs$cut$Chr]), start = p2obj$cbs$cut$Start, end = p2obj$cbs$cut$End, value = p2obj$cbs$cut$Log2Ratio) else data.frame(chr = p2obj$segments$chr, start = p2obj$segments$startpos, end = p2obj$segments$endpos, value = p2obj$segments$nMajor + p2obj$segments$nMinor)
  prof2_ploidy <- if (prof_mod == 'L2R') NULL else p2obj$ploidy$most.width
  ## Synch chromosomes
  p1seg <- p1seg[p1seg$chr %in% unique(p2seg$chr),]
  p2seg <- p2seg[p2seg$chr %in% unique(p1seg$chr),]
  ## Restrict cs
  avail_chr <- unique(c(p1seg$chr, p2seg$chr))
  cs <- chromobjector(BSg = BSg.obj, chrs = avail_chr)
  ## Common variables
  seg_comp(prof1_df = p1seg, prof2_df = p2seg, prof1_ploidy = prof1_ploidy, prof2_ploidy = prof2_ploidy, prof1_name = des_df$prof1_name[d], prof2_name = des_df$prof2_name[d], genome_package = genome_package, prof_mod = prof_mod, prof1_l2r_recenter = prof1_l2r_recenter, prof2_l2r_recenter = prof2_l2r_recenter, capping = capping, epsilon = 1E-16, plot = TRUE, prof1_col = prof1_col, prof2_col = prof2_col, cytobandideo_file = cytobandideo_file, capture_bed_file = capture_bed_file, out_dir = od2)
}


#########################
## EACON.ASCAT vs FACETS
#########################

  ### Mode : EACON.ASCAT.WES vs FACETS.WES // L2R
  techno1 <- 'WES'
  techno2 <- 'WES'
  seg_colname <- 'CNLR_MEDIAN_CLUST'
  prof_mod = 'L2R'
  capping = 1.5
  prof1_l2r_recenter <- FALSE
  prof2_l2r_recenter <- TRUE
  
  ### Mode : EACON.ASCAT.WES vs FACETS.WES // TCN
  techno1 <- 'WES'
  techno2 <- 'WES'
  seg_colname <- 'TCN_EM'
  prof_mod = 'TCN'
  capping = 6
  prof1_l2r_recenter <- FALSE
  prof2_l2r_recenter <- FALSE

  ### Mode : EACON.ASCAT.WGS vs FACETS.WGS // L2R
  techno1 <- 'WGS'
  techno2 <- 'WGS'
  seg_colname <- 'CNLR_MEDIAN_CLUST'
  prof_mod = 'L2R'
  capping = 1.5
  prof1_l2r_recenter <- FALSE
  prof2_l2r_recenter <- TRUE
  
  ### Mode : EACON.ASCAT.WGS vs FACETS.WGS // TCN
  techno1 <- 'WGS'
  techno2 <- 'WGS'
  seg_colname <- 'TCN_EM'
  prof_mod = 'TCN'
  capping = 6
  prof1_l2r_recenter <- FALSE
  prof2_l2r_recenter <- FALSE
  
## RUN
des_file <- paste0('/home/job/WORKSPACE/ElsaLAB/CNA/seg_comp/EACON.ASCAT.', techno1, '_vs_FACETS.', techno2, '_design_', prof_mod, '.tsv')
des_df <- data.table::fread(file = des_file, sep = '\t', data.table = FALSE)
od2 <- paste(c(out_dir, prof_mod, paste0('EACON.ASCAT.', techno1, '_vs_FACETS.', techno2)), collapse = '/')
dir.create(path = od2, recursive = TRUE)
for (d in seq_len(nrow(des_df))) {
  cs <- chromobjector(BSg = BSg.obj)
  message(des_df$prof1_name[d], ' vs ', des_df$prof2_name[d])
  ## Load EaCoN data 1
  p1obj <- readRDS(des_df$prof1_file[d])
  p1seg <- if (prof_mod == 'L2R') data.frame(chr = unlist(cs$chr2chrom[p1obj$cbs$cut$Chr]), start = p1obj$cbs$cut$Start, end = p1obj$cbs$cut$End, value = p1obj$cbs$cut$Log2Ratio) else data.frame(chr = p1obj$segments$chr, start = p1obj$segments$startpos, end = p1obj$segments$endpos, value = p1obj$segments$nMajor + p1obj$segments$nMinor)
  prof1_ploidy <- if (prof_mod == 'L2R') NULL else p1obj$ploidy$most.width
  ## Load FACETS data 2
  p2 <- facets_vcf_reader(vcf_file = des_df$prof2_file[d], seg_colname = seg_colname)
  prof2_ploidy <- if (prof_mod == 'L2R') NULL else p2$meta$ploidy
  ## Synch chromosomes
  p1seg <- p1seg[p1seg$chr %in% unique(p2$seg$chr),]
  p2$seg <- p2$seg[p2$seg$chr %in% unique(p1seg$chr),]
  ## Restrict cs
  avail_chr <- unique(c(p1seg$chr, p2$seg$chr))
  cs <- chromobjector(BSg = BSg.obj, chrs = avail_chr)
  ### Go
  seg_comp(prof1_df = p1seg, prof2_df = p2$seg, prof1_ploidy = prof1_ploidy, prof2_ploidy = prof2_ploidy, prof1_name = des_df$prof1_name[d], prof2_name = des_df$prof2_name[d], genome_package = genome_package, prof_mod = prof_mod, prof1_l2r_recenter = prof1_l2r_recenter, prof2_l2r_recenter = prof2_l2r_recenter, capping = capping, epsilon = 1E-16, plot = TRUE, prof1_col = prof1_col, prof2_col = prof2_col, cytobandideo_file = cytobandideo_file, capture_bed_file = capture_bed_file, out_dir = od2)
}










