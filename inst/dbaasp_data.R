# peptides-complete from https://dbaasp.org/download-dataset?source=peptides
# accessed on 22-04-2022
#' This function pulls the data related to the aa's distribution
#'
#' @import dplyr
#' @importFrom utils read.csv write.csv
#'
#' @details
#' peptides-complete from \url{https://dbaasp.org/download-dataset?source=peptides}
#' accessed on 22-04-2022
#'

get_dbaasp_data <- function() {
    dbaasp <- read.csv("peptides-complete.csv")
    dbaasp %>%
        dplyr::mutate(AMP = grepl("Gram+", x = TARGET.GROUP) |
                          grepl("Gram-", x = TARGET.GROUP),
                      ACP = grepl("Cancer", x = TARGET.GROUP)) %>%
        dplyr::filter(ACP | AMP) %>%
        dplyr::select(ID, SEQUENCE, ACP, AMP) %>%
        unique %>%
        write.csv("./data/dbaasp_amp_acp.csv", row.names = FALSE)
}
