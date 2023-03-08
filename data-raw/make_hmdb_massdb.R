library(massdatabase)
library(tidyr)
library(dplyr)
library(utils)
library(usethis)

data.dir <- "data-raw"

### hmdb database version 5.0, all metabolites xml
### url: https://hmdb.ca/downloads
### only run for first time
# hmdb.file <- "hmdb_metabolites.zip"
# utils::unzip(file.path(data.dir, hmdb.file), exdir = data.dir)
# hmdb.xml <- massdatabase::read_xml_data_hmdb("hmdb_metabolites.xml", path = data.dir)
# hmdb.massdb <- massdatabase::convert_hmdb2metid(hmdb.xml)

### download and read lipidmaps database
### only run for first time
# massdatabase::download_lipidmaps_lipid(path = data.dir)
# utils::unzip(file.path(data.dir, "LMSD.sdf.zip"), exdir = data.dir)
# lipidmaps.sdf <- massdatabase::read_sdf_data_lipidmaps(file.path(data.dir, "structures.sdf"), path = ".")
# lipidmaps.massdb <- massdatabase::convert_lipidmaps2metid(lipidmaps.sdf)

# save raw db
# save(hmdb.xml, lipidmaps.sdf, file = file.path(data.dir, "raw_db.Rda"))
# save(hmdb.massdb, lipidmaps.massdb, file = file.path(data.dir, "read_massdb.Rda"))

load(file.path(data.dir, "read_massdb.Rda"))

hmdb.tb <- tibble::as_tibble(hmdb.massdb@spectra.info)
lipm.tb <- tibble::as_tibble(lipidmaps.massdb@spectra.info)

## select HMDB entry with HMDB.ID/KEGG.ID/PUBCHEM.ID/INCHIKEY.ID/CHEBI.ID in lipidmap
hmdb.with.lipm.hmdb <- hmdb.tb %>% dplyr::select(Lab.ID, HMDB.ID) %>%
  dplyr::inner_join(lipm.tb %>% dplyr::select(HMDB.ID, LIPM.ID = Lab.ID), 
                    by = "HMDB.ID", multiple = "all")
hmdb.with.lipm.kegg <- hmdb.tb %>% dplyr::select(Lab.ID, KEGG.ID) %>%
  dplyr::inner_join(lipm.tb %>% dplyr::filter(!is.na(KEGG.ID)) %>% dplyr::select(KEGG.ID, LIPM.ID = Lab.ID), 
                    by = "KEGG.ID", multiple = "all")
hmdb.with.lipm.pubchem <- hmdb.tb %>% dplyr::select(Lab.ID, PUBCHEM.ID) %>%
  dplyr::inner_join(lipm.tb %>% dplyr::filter(!is.na(PUBCHEM.ID)) %>% dplyr::select(PUBCHEM.ID, LIPM.ID = Lab.ID), 
                    by = "PUBCHEM.ID", multiple = "all")
hmdb.with.lipm.inchi <- hmdb.tb %>% dplyr::select(Lab.ID, INCHIKEY.ID) %>%
  dplyr::inner_join(lipm.tb %>% dplyr::filter(!is.na(INCHIKEY.ID)) %>% dplyr::select(INCHIKEY.ID, LIPM.ID = Lab.ID), 
                    by = "INCHIKEY.ID", multiple = "all")
hmdb.with.lipm.chebi <- hmdb.tb %>% dplyr::select(Lab.ID, CHEBI.ID) %>%
  dplyr::inner_join(lipm.tb %>% dplyr::filter(!is.na(CHEBI.ID)) %>% dplyr::select(CHEBI.ID, LIPM.ID = Lab.ID), 
                    by = "CHEBI.ID", multiple = "all")

hmdb.in_lipm <- dplyr::bind_rows(hmdb.with.lipm.hmdb,
                                 hmdb.with.lipm.kegg,
                                 hmdb.with.lipm.pubchem,
                                 hmdb.with.lipm.inchi,
                                 hmdb.with.lipm.chebi) %>%
  dplyr::select(Lab.ID, LIPM.ID) %>% dplyr::group_by(Lab.ID) %>%
  dplyr::distinct() %>%
  dplyr::summarise(LIPM.ID = paste(LIPM.ID, collapse = ";"))

# filter and add lipidmap ID to HMDB as lipids
hmdb.lipids.tb <- hmdb.tb %>% dplyr::inner_join(hmdb.in_lipm, by = "Lab.ID")

# filter non lipids as mets
hmdb.mets.tb <- hmdb.tb %>% dplyr::filter(!(Lab.ID %in% hmdb.lipids.tb$Lab.ID))

hmdb.lipids.massdb <- hmdb.massdb
hmdb.lipids.massdb@spectra.info <- as.data.frame(hmdb.lipids.tb)
hmdb.lipids.massdb@database.info$version <- Sys.Date()
hmdb.lipids.massdb@database.info$Source <- "HMDB with LIPIDMAPS"

hmdb.mets.massdb <- hmdb.massdb
hmdb.mets.massdb@spectra.info <- as.data.frame(hmdb.mets.tb)
hmdb.mets.massdb@database.info$version <- Sys.Date()
hmdb.mets.massdb@database.info$Source <- "HMDB without LIPIDMAPS"

usethis::use_data(hmdb.lipids.massdb, compress = "xz", overwrite = T)
usethis::use_data(hmdb.mets.massdb, compress = "xz", overwrite = T)

