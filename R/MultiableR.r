library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(limma)
library(ProteoMM)
library(xcms)
library(CAMERA)
library(metid)

#' @export
runXCMS <- function (mzMLs, 
                     sample.info, 
                     polarity, 
                     save.xdata = FALSE, 
                     save.location = NA, 
                     cwp = NA, 
                     owp = NA, 
                     pdp = NA,
                     camera.sigma = 6,
                     camera.perfwhm = 0.6,
                     camera.intval = "into", 
                     camera.ppm = 10,
                     camera.mzabs = 0.015,
                     camera.maxcharge = 3,
                     camera.maxiso = 5
) {
  if (is.na(cwp)) {
    cwp <- xcms::CentWaveParam(ppm = 10,
                         snthr = 6,
                         peakwidth = c(10, 60), 
                         mzdiff = 0.01,
                         noise = 0,
                         prefilter = c(3, 500))
  } else if (!isClass(cwp, "CentWaveParam")) {
    stop("cwp must be class CentWaveParam!")
  }
  
  if (is.na(owp)) {
    owp <- xcms::ObiwarpParam(binSize = 0.5)
  } else if (!isClass(owp, "ObiwarpParam")) {
    stop("owp must be class ObiwarpParam!")
  }
  
  if(is.null(sample.info$sample)) {
    stop("sample.info must contain column `sample`!")
  }else if(is.null(sample.info$group)) {
    stop("sample.info must contain column `group`!")
  } else if (!all(sample.info$group %in% c("sample", "blank"))) {
    stop('sample.info$group column contain labels other than "sample" or "blank"!')
  }
  
  if (is.na(pdp)) {
    pdp <- xcms::PeakDensityParam(sampleGroups = sample.info$group,
                            bw = 5,
                            binSize = 0.025,
                            minFraction = 0.5,
                            minSamples = 1)
  } else if (!isClass(pdp, "PeakDensityParam")) {
    stop("pdp must be class PeakDensityParam!")
  }
  
  if(save.xdata && is.na(save.location)) {
    stop("`save.xdata`` is set but no location is found. Please provide location at `save.location`.")
  } else if (save.xdata) {
    if (file.exists(save.location)) {
      loop <- TRUE
      continue_script <- FALSE
      while(loop) {
        rl <- readline(prompt="Warning: directory already exists. Proceed and overwrite existing xdata? [y/n]") 
        if (tolower(rl) == "y") {
          continue_script <- TRUE
          loop <- FALSE
        } else if (tolower(tl == "n")) {
          continue_script <- FALSE
          loop <- FALSE
        } else {
          loop <- TRUE
        }
      }
      if (!continue_script) {
        stop("Manual stop because directory already exists.")
      }
    } else {
      dir.create(save.location, recursive = TRUE)
    }
  }
  
  # reading MSData
  raw_data <- MSnbase::readMSData(files = mzMLs, pdata = new("NAnnotatedDataFrame", sample.info), mode = "onDisk")
  
  # detecting chromatographic peaks
  xdata <- xcms::findChromPeaks(raw_data, param = cwp)
  # aligning retention time
  xdata2 <- xcms::adjustRtime(xdata, param = owp)
  # grouping chrompeaks using sample vs blank
  xdata3 <- xcms::groupChromPeaks(xdata2, param = pdp)
  
  if (save.xdata) {
    xdata.path <- file.path(save.location, "xdata.Rdata")
    xdata2.path <- file.path(save.location, "xdata2.Rdata")
    xdata3.path <- file.path(save.location, "xdata3.Rdata")
    
    save(xdata, file = xdata.path)
    save(xdata2, file = xdata2.path)
    save(xdata3, file = xdata3.path)
  }
  
  polarity <- tolower(polarity)
  primary.adducts.file <- paste0("rules/primary_adducts_", stringr::str_extract(polarity, "^.{3}"),".csv")
  xsa <- CAMERA::xsAnnotate(as(MSnbase::filterMsLevel(xdata3, msLevel. = 1), "xcmsSet"), polarity = polarity)
  anF <- CAMERA::groupFWHM(xsa, sigma = camera.sigma, perfwhm = camera.perfwhm, intval = camera.intval)
  anI <- CAMERA::findIsotopes(anF, ppm = camera.ppm, mzabs = camera.mzabs, maxcharge = camera.maxcharge, maxiso = camera.maxiso, intval = camera.intval)
  anIC <- CAMERA::groupCorr(anI)
  anFA <- CAMERA::findAdducts(anIC, polarity=polarity, rules = read.csv(system.file(primary.adducts.file, package = "CAMERA")))
  peaklist <- CAMERA::getPeaklist(anFA)
  rownames(peaklist) <- rownames(xcms::featureDefinitions(xdata3))
  
  xdata.table <- xcms::featureValues(xdata3) %>% 
    tidyr::as_tibble(rownames = "Feature.ID") %>% 
    tidyr::pivot_longer(-Feature.ID) %>% 
    dplyr::mutate(name = stringr::str_remove(name, pattern = ".mzML")) %>%
    tidyr::pivot_wider()
  xdata.ftInfo <- xcms::featureDefinitions(xdata3) %>% 
    tidyr::as_tibble(rownames = "Feature.ID") %>% 
    dplyr::select(Feature.ID, mzmed, rtmed)
  camera.table <- peaklist %>% tidyr::as_tibble(rownames = "Feature.ID") %>% 
    dplyr::select(Feature.ID, isotopes, adduct, pcgroup)
  annotated.table <- xdata.ftInfo %>% dplyr::inner_join(xdata.table, by = "Feature.ID") %>%
    dplyr::inner_join(camera.table, by = "Feature.ID") %>%
    dplyr::select(Feature.ID, mzmed, rtmed, dplyr::filter(sample.info, group == "sample")$sample,
           dplyr::filter(sample.info, group == "blank")$sample, isotopes, adduct, pcgroup)
  
  return(annotated.table)
}

#' @export
metid_annotate <- function (df, polarity, omic = "Mx", database = NA, ms1.match.ppm = 15,
                            rt.match.tol = 1000000) {
  load("database.RData")
  if (!is.na(database)) {
    if (class(database) != "databaseClass") {
      stop("Database must be class of metid::databaseClass!")
    }
  } else {
    if (tolower(omic) == "lx") {
      data("hmdb.lipids.massdb")
      database <- hmdb.lipids.massdb
    } else if (tolower(omic) == "mx") {
      data("hmdb.mets.massdb")
      database <- hmdb.mets.massdb
    }
    else {
      stop("Unidentified omic: Lx or Mx; or provide manual database")
    }
  }
  if (!("Feature.ID" %in% colnames(df))) {
    stop("Cannot find column `Feature.ID` in dataframe!")
  }
  if ("mzmed" %in% colnames(df)) {
    df <- df %>% dplyr::rename(mz = mzmed)
  } else if (!("mz" %in% colnames(df))) {
    stop("Cannot find column `mz` or `mzmed` in dataframe!")
  }
  if ("rtmed" %in% colnames(df)) {
    df <- df %>% dplyr::rename(rt = rtmed)
  } else if (!("rt" %in% colnames(df))) {
    stop("Cannot find column `rt` or `rtmed` in dataframe!")
  }
  tmp_file <- tempfile("metid", fileext = ".csv")
  for.metid.file <- basename(tmp_file)
  on.exit(unlink(for.metid.file), add = TRUE)
  df %>% dplyr::select(Feature.ID, mz, rt) %>% readr::write_csv(tmp_file)
  
  polarity <- tolower(polarity)
  annotate_result <- metid::identify_metabolites(
    ms1.data = for.metid.file,
    ms1.match.ppm = ms1.match.ppm, 
    rt.match.tol = rt.match.tol,
    path = tempdir(),
    polarity = polarity,
    database = database
  )
  
  annotate_table <- annotate_result %>% metid::get_identification_table(type = "new", candidate.num = 1) %>%
    dplyr::rename(Feature.ID = name)
  
  return (annotate_table)
}

#' @export
write_lipidfinder_csv <- function(df, file.name) {
  lf_output <- df %>% replace(is.na(.), 0) %>%
    dplyr::mutate(Feature.ID = as.numeric(stringr::str_extract(Feature.ID, "([0-9]+)", group = 1))) %>%
    dplyr::rename(FeatureID = Feature.ID)
  
  readr::write_csv(lf_output, file = file.name)
}

#' Read csv files from LipidFinder output
#'
#' Convert LipidFinder output into a dataframe
#' The function also automatically process the sample name
#' File format: Please make sure the column "Polarity" is included
#' and a column named "FeatureID" for indexing the column
#' @param file.path The file location
#' @param omic omic the file belongs to: accept Mx / Lx
#' @return a tibble dataframe
#' @export
read_LipidFinder_csv <- function(file.path, omic = "Mx") {
    df <- readr::read_csv(file.path)
    if (omic == "Lx") {
        df <- df %>% dplyr::mutate(Omic = "Lx")
    } else if (omic == "Mx") {
        df <- df %>% dplyr::mutate(Omic = "Mx")
    } else {
        stop("Unknown omic profile. Omic options:[\"Lx\", \"Mx\"]")
    }
    fid_length <- nchar(max(df$FeatureID))
    df <- df %>%
      dplyr::mutate(FeatureID = paste0(stringr::str_extract(Omic, "^(.)"),
                                       stringr::str_extract(Polarity, "^(.)"),
                                       stringr::str_pad(FeatureID, fid_length, side = "left", pad = "0"))) %>%
      dplyr::select(Omic, Polarity, tidyselect::everything()) %>%
      dplyr::rename(Feature.ID = FeatureID) %>%
      tidyr::pivot_longer(ends_with("mean"), names_to = "sample") %>%
      dplyr::mutate(sample = stringr::str_replace(sample, "_mean", "")) %>%
      tidyr::pivot_wider(names_from = "sample")
    return(df)
}

#' Data frame preprocesssing
#'
#' The function preprocess the data by filtering, imputing, and transforming
#' @param df the tibble containing only Feature.IDs and sample abundance
#' @param filter_threshold minimum number of non NAs to keep the feature
#' @param impute imputation method. "gm2" = half of global minimum, "gm5" = a fifth of global minimum
#' @param log_base log transform base. Default = 2
#' @return processed and transformed tibble data frame
#' @export
tbManipulate <- function(df, filter_threshold, impute = "gm2", log_base = 2) {
    df.m <- df %>% tibble::column_to_rownames("Feature.ID") %>% as.matrix()
    thres <- filter_threshold
    filtered.df <- df.m[rowSums(!is.na(df.m)) > thres,]

    if(impute == "gm2" || impute == "gm5") {
        global.min <- min(filtered.df, na.rm = TRUE)
        if (impute == "gm2") {
            impute.value <- global.min / 2
        }
        if (impute == "gm5") {
            impute.value <- global.min / 5
        }
        filtered.df[is.na(filtered.df)] <- impute.value
        imputed.df <- filtered.df
    }

    log.df <- log(imputed.df, base = log_base)

    return(tibble::as_tibble(log.df, rownames = "Feature.ID"))
}

#' EigenMS normalisation
#'
#' Use the ProteoMM package to normalise the dataframe using EigenMS
#' @param df the tibble containing only Feature.IDs and sample abundance
#' @param sample.info sample info. Must contain a column named "sample"
#' for the sample information and a column named "group" for the experimental group
#' @return a EigenMS normalised tibble dataframe
#' @export
normEigenMS <- function(df, sample.info) {
    if (!all(sort(sample.info$sample) == sort(colnames(df %>% dplyr::select(-Feature.ID))))) {
         stop('sample name not match between sample.info and df')
    }

    m_prot.info <- df %>% dplyr::select(Feature.ID) %>% as.data.frame()

    sample.treatment <- as.factor(sample.info$group)
    m_logInts <- df %>% dplyr::select(sample.info$sample) %>% as.data.frame()

    hs_m_ints_eig1 <- ProteoMM::eig_norm1(m=m_logInts,treatment=sample.treatment,prot.info=m_prot.info)
    if(hs_m_ints_eig1$h.c > 0) {
        hs_m_ints_norm <- ProteoMM::eig_norm2(rv=hs_m_ints_eig1)
        output <- tibble::as_tibble(hs_m_ints_norm$norm_m, rownames = "Feature.ID")
    } else {
        message("WARNING: No bias trends detected!")
        output <- df
    }
    return(output)
}

#' Median normalisation
#'
#' Median normalise the dataframe
#' @param df the tibble containing only Feature.IDs and sample abundance
#' @param mean the average value to normalise the dataframe to. Default as NA
#' which normalise the mean of the dataframe to the mean of the medians
#' @return a Median normalised tibble dataframe
#' @export
normMedian <- function(df, mean = NA) {
    rowMeans <- df %>% dplyr::summarise_if(is.numeric, median, na.rm = TRUE) %>% tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "median")
    if (is.na(mean)) {
        overallMean <- rowMeans %>% dplyr::summarise(mean = mean(median)) %>% .$mean
    } else {
        overallMean <- mean
    }
    median.df <- df %>% tidyr::pivot_longer(-Feature.ID) %>% dplyr::inner_join(rowMeans, by = "name") %>%
                        dplyr::mutate(value = value - median + overallMean) %>% dplyr::select(-median) %>% tidyr::pivot_wider()
    return(median.df)
}

#' Limma differential analysis
#'
#' Perform limma differential analysis
#' @param df the tibble containing only Feature.IDs and sample abundance
#' @param design.matrix limma design matrix. See limma
#' @param cont.matrix limma cont matrix. See limma
#' @return The log-fold change and p value of each feature per comparison
#' @export
runLimma <- function(df, design.matrix, cont.matrix) {
    limma.fit <- df %>% tibble::column_to_rownames(var = "Feature.ID") %>% as.matrix() %>%
      limma::lmFit(design.matrix) %>%
      limma::contrasts.fit(cont.matrix) %>%
      limma::eBayes()

    comparisons <- colnames(cont.matrix)
    limma.tb.list <- list()
    for(i in 1:length(comparisons)) {
        limma.tb <- limma.fit %>% limma::topTable(coef = comparisons[i], adjust.method = "BH", number = Inf) %>% tidyr::as_tibble(rownames = "Feature.ID") %>% dplyr::mutate(coef = comparisons[i])
        limma.tb.list[[i]] <- limma.tb
    }
    limma.summary <- dplyr::bind_rows(limma.tb.list) %>% dplyr::select(Feature.ID, coef, logFC, adj.P.Val) %>%
            tidyr::pivot_longer(c(-Feature.ID, -coef)) %>% dplyr::mutate(colname = paste(name, coef, sep = ".")) %>%
            dplyr::select(Feature.ID, colname, value) %>% tidyr::pivot_wider(names_from = "colname")
    return(limma.summary)
}

#' @export
group_annotation <- function(df, annotated.df, match.id = "HMDB.ID") {
  if (!(match.id %in% colnames(annotated.df))) {
    stop(paste0("Column `", match.id, "` not found in annotation dataframe!"))
  }
  
  grouped.df = annotated.df %>% 
    dplyr::rename(Match.ID = tidyselect::matches(paste0("^",match.id,"$"))) %>% 
    dplyr::select(Omic, Feature.ID, Match.ID) %>% 
    dplyr::inner_join(df %>% dplyr::ungroup() %>% dplyr::select(-Omic, -Polarity), by = "Feature.ID") %>%
    dplyr::filter(!is.na(Match.ID)) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(rowmed = median(dplyr::c_across(tidyselect::where(is.numeric))),
                  rowmean = mean(dplyr::c_across(tidyselect::where(is.numeric))),
                  rowmax = max(dplyr::c_across(tidyselect::where(is.numeric)))) %>% 
    dplyr::group_by(Match.ID) %>% 
    dplyr::slice_max(rowmed) %>% 
    dplyr::slice_max(rowmean) %>%
    dplyr::slice_max(rowmax, n = 1) %>%
    dplyr::select(-c(rowmean, rowmed, rowmax, Feature.ID)) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(!!match.id := Match.ID)
  return(grouped.df)
}
