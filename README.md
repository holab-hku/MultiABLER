# MultiABLER: An R package for integrative analysis of metabolomic and lipidomic data from MultiABLE

MultiABLER is a set of R functions forms a seamless workflow that supports integrative processing and analysis of lipidomic and metabolomic data generated by the Multi-ABLE method ([Talib et al., 2019](https://pubs.acs.org/doi/10.1021/acs.analchem.9b01842)). To run MultiABLER, install the following packges in R and run the funcions in MultiABLER.r. A detailed tutorial on how to run MultiABLER can be found in [here](
https://htmlpreview.github.io/?https://github.com/holab-hku/MultiABLER/blob/main/tutorials/tutorial.html).

## Installation

To install via github, run the following code in R:

`devtools::install_github("holab-hku/MultiABLER", dependencies = TRUE)`

## Installation issues

In case of installation issues, try install the following dependencies and run the above code:

* [ProteoMM](https://www.bioconductor.org/packages/release/bioc/html/ProteoMM.html)
* [limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html)
* [Tidyverse](https://www.tidyverse.org)

## License

MIT
