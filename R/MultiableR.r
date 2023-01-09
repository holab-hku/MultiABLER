library(tidyverse)
library(limma)
library(ProteoMM)
library(xcms)
library(CAMERA)
library(metid)

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
