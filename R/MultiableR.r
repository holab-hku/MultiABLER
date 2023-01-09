library(tidyverse)
library(limma)
library(ProteoMM)
library(xcms)
library(CAMERA)
library(metid)

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
    cwp <- CentWaveParam(ppm = 10,
                         snthr = 6,
                         peakwidth = c(10, 60), 
                         mzdiff = 0.01,
                         noise = 0,
                         prefilter = c(3, 500))
  } else if (!isClass(cwp, "CentWaveParam")) {
    stop("cwp must be class CentWaveParam!")
  }
  
  if (is.na(owp)) {
    owp <- ObiwarpParam(binSize = 0.5)
  } else if (!isClass(owp, "ObiwarpParam")) {
    stop("owp must be class ObiwarpParam!")
  }
  
  if(is.null(sample.info$group)) {
    stop("sample.info must contain column `group`!")
  }
  
  if (is.na(pdp)) {
    pdp <- PeakDensityParam(sampleGroups = sample.info$group,
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
  
  raw_data <- readMSData(files = mzMLs, pdata = new("NAnnotatedDataFrame", sample.info), mode = "onDisk")
  
  xdata <- findChromPeaks(raw_data, param = cwp)
  xdata2 <- adjustRtime(xdata, param = owp)
  xdata3 <- groupChromPeaks(xdata2, param = pdp)
  
  if (save.xdata) {
    xdata.path <- file.path(save.location, "xdata.Rdata")
    xdata2.path <- file.path(save.location, "xdata2.Rdata")
    xdata3.path <- file.path(save.location, "xdata3.Rdata")
    
    save(xdata, file = xdata.path)
    save(xdata2, file = xdata2.path)
    save(xdata3, file = xdata3.path)
  }
  
  primary_adducts <- paste0("rules/primary_adducts_", str_extract(polarity, "^.{3}"),".csv")
  xsa <- xsAnnotate(as(filterMsLevel(xdata3, msLevel. = 1), "xcmsSet"), polarity = polarity)
  anF <- groupFWHM(xsa, sigma = camera.sigma, perfwhm = camera.perfwhm, intval = camera.intval)
  anI <- findIsotopes(anF, ppm = camera.ppm, mzabs = camera.mzabs, maxcharge = camera.maxcharge, maxiso = camera.maxiso, intval = camera.intval)
  anIC <- groupCorr(anI)
  anFA <- findAdducts(anIC, polarity=polarity, rules = read.csv(system.file(primary_adducts, package = "CAMERA")))
  peaklist <- getPeaklist(anFA)
  rownames(peaklist) <- rownames(featureDefinitions(xdata3))
  
  xdata.table <- featureValues(xdata3) %>% as_tibble(rownames = "Feature.ID")
  xdata.table <- xdata.table %>% pivot_longer(-Feature.ID) %>% mutate(name = str_remove(name, pattern = ".mzML")) %>% pivot_wider()
  xdata.ftInfo <- featureDefinitions(xdata3) %>% as_tibble(rownames = "Feature.ID")
  camera.table <- peaklist %>% as_tibble(rownames = "Feature.ID")
  annotated.table <- xdata.ftInfo %>% select(Feature.ID, mzmed, rtmed) %>% inner_join(xdata.table, by = "Feature.ID") %>%
    inner_join(camera.table %>% select(Feature.ID, isotopes, adduct, pcgroup), by = "Feature.ID") %>%
    select(Feature.ID, mzmed, rtmed, sample.info %>% filter(group == "sample") %>% .$file,
           sample.info %>% filter(group == "blank") %>% .$file, isotopes, adduct, pcgroup)
  
  return(annotated.table)
  
}

metid_annotate <- function (df, polarity, omic = "Mx") {
  load("database.RData")
  database <- if (tolower(omic) == "lx") hmdb_lipids_metid else hmdb_metabo_metid
  tmp_file <- tempfile("metid", fileext = ".csv")
  for.metid.file <- basename(tmp_file)
  df %>% select(Feature.ID, mz, rt) %>% write_csv(tmp_file)
  
  annotate_result <- metid::identify_metabolites(
    ms1.data = for.metid.file,
    ms1.match.ppm = 15, 
    rt.match.tol = 1000000,
    path = tempdir(),
    polarity = polarity,
    database = database
  )
  
  annotate_table <- annotate_result %>% get_identification_table(type = "new", candidate.num = 1) %>% rename(FeatureID = name)
  
  return (annotate_table)
}

write_lipidfinder_csv <- function(df, file.name) {
  lf_output <- df %>% replace(is.na(.), 0) %>%
    mutate(Feature.ID = as.numeric(str_replace(Feature.ID, "FT", ""))) %>%
    rename(FeatureID = Feature.ID)
  
  write_csv(lf_output, file = file.name)
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
    df <- read_csv(file.path)
    if (omic == "Lx") {
        df <- df %>% mutate(Omic = "Lx")
    } else if (omic == "Mx") {
        df <- df %>% mutate(Omic = "Mx")
    } else {
        stop("Unknown omic profile. Omic options:[\"Lx\", \"Mx\"]")
    }
    fid_length <- nchar(max(df$FeatureID))
    df <- df %>%
            mutate(FeatureID = paste0(str_extract(Omic, "^(.)"),
                                      str_extract(Polarity, "^(.)"),
                                      str_pad(FeatureID, fid_length, side = "left", pad = "0"))) %>%
            select(Omic, Polarity, everything()) %>%
            pivot_longer(ends_with("mean"), names_to = "sample") %>%
            mutate(sample = str_replace(sample, "_mean", "")) %>%
            pivot_wider(names_from = "sample")
    return(df)
}

#' Data frame preprocesssing
#'
#' The function preprocess the data by filtering, imputing, and transforming
#' @param df the tibble containing only FeatureIDs and sample abundance
#' @param filter_threshold minimum number of non NAs to keep the feature
#' @param impute imputation method. "gm2" = half of global minimum, "gm5" = a fifth of global minimum
#' @param log_base log transform base. Default = 2
#' @return processed and transformed tibble data frame
#' @export
tbManipulate <- function(df, filter_threshold, impute = "gm2", log_base = 2) {
    df.m <- df %>% column_to_rownames("FeatureID") %>% as.matrix()
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

    return(as_tibble(log.df, rownames = "FeatureID"))
}

#' EigenMS normalisation
#'
#' Use the ProteoMM package to normalise the dataframe using EigenMS
#' @param df the tibble containing only FeatureIDs and sample abundance
#' @param sample.info sample info. Must contain a column named "sample"
#' for the sample information and a column named "group" for the experimental group
#' @return a EigenMS normalised tibble dataframe
#' @export
normEigenMS <- function(df, sample.info) {
    if (!all(sort(sample.info$sample) == sort(colnames(df %>% select(-FeatureID))))) {
         stop('sample name not match between sample.info and df')
    }

    m_prot.info <- df %>% select(FeatureID) %>% as.data.frame()

    sample.treatment <- as.factor(sample.info$group)
    m_logInts <- df %>% select(sample.info$sample) %>% as.data.frame()

    hs_m_ints_eig1 <- ProteoMM::eig_norm1(m=m_logInts,treatment=sample.treatment,prot.info=m_prot.info)
    if(hs_m_ints_eig1$h.c > 0) {
        hs_m_ints_norm <- ProteoMM::eig_norm2(rv=hs_m_ints_eig1)
        output <- as_tibble(hs_m_ints_norm$norm_m, rownames = "FeatureID")
    } else {
        message("WARNING: No bias trends detected!")
        output <- df
    }
    return(output)
}

#' Median normalisation
#'
#' Median normalise the dataframe
#' @param df the tibble containing only FeatureIDs and sample abundance
#' @param mean the average value to normalise the dataframe to. Default as NA
#' which normalise the mean of the dataframe to the mean of the medians
#' @return a Median normalised tibble dataframe
#' @export
normMedian <- function(df, mean = NA) {
    rowMeans <- df %>% summarise_if(is.numeric, median, na.rm = TRUE) %>% pivot_longer(cols = everything(), values_to = "median")
    if (is.na(mean)) {
        overallMean <- rowMeans %>% summarise(mean = mean(median)) %>% .$mean
    } else {
        overallMean <- mean
    }
    median.df <- df %>% pivot_longer(-FeatureID) %>% inner_join(rowMeans, by = "name") %>%
                        mutate(value = value - median + overallMean) %>% select(-median) %>% pivot_wider()
    return(median.df)
}

#' Limma differential analysis
#'
#' Perform limma differential analysis
#' @param df the tibble containing only FeatureIDs and sample abundance
#' @param design.matrix limma design matrix. See limma
#' @param cont.matrix limma cont matrix. See limma
#' @return The log-fold change and p value of each feature per comparison
#' @export
runLimma <- function(df, design.matrix, cont.matrix) {
    limma.fit <- df %>% column_to_rownames(var = "FeatureID") %>% as.matrix() %>%
      limma::lmFit(design.matrix) %>%
      limma::contrasts.fit(cont.matrix) %>%
      limma::eBayes()

    comparisons <- colnames(cont.matrix)
    limma.tb.list <- list()
    for(i in 1:length(comparisons)) {
        limma.tb <- limma.fit %>% limma::topTable(coef = comparisons[i], adjust.method = "BH", number = Inf) %>% as_tibble(rownames = "FeatureID") %>% mutate(coef = comparisons[i])
        limma.tb.list[[i]] <- limma.tb
    }
    limma.summary <- bind_rows(limma.tb.list) %>% select(FeatureID, coef, logFC, adj.P.Val) %>%
            pivot_longer(c(-FeatureID, -coef)) %>% mutate(colname = paste(name, coef, sep = ".")) %>%
            select(FeatureID, colname, value) %>% pivot_wider(names_from = "colname")
    return(limma.summary)
}
