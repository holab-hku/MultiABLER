% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MultiableR.r
\name{read_LipidFinder_csv}
\alias{read_LipidFinder_csv}
\title{Read csv files from LipidFinder output}
\usage{
read_LipidFinder_csv(file.path, omic = "Mx")
}
\arguments{
\item{file.path}{The file location}

\item{omic}{omic the file belongs to: accept Mx / Lx}
}
\value{
a tibble dataframe
}
\description{
Convert LipidFinder output into a dataframe
The function also automatically process the sample name
File format: Please make sure the column "Polarity" is included
and a column named "FeatureID" for indexing the column
}
