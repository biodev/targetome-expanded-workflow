#' Create a UniProt ID Mapping Table
#'
#' This function reads and processes the UniProt ID mapping file
#' (`HUMAN_9606_idmapping.dat.gz`) to produce a cleaned mapping table
#' between UniProt accession IDs and other identifier systems.
#' It standardizes the identifiers by removing version numbers and
#' isoform suffixes, and adds self-mapping entries for all UniProt IDs.
#'
#' @param unip.file Character string giving the path to the UniProt
#' ID mapping file, typically `"HUMAN_9606_idmapping.dat.gz"`.
#' This file should be the human ID mapping dataset obtained from UniProt,
#' containing tab-delimited fields of UniProt accessions and corresponding IDs.
#'
#' @return
#' A `data.table` containing the cleaned mapping data with columns:
#' \describe{
#'   \item{V1}{Cleaned UniProt accession IDs.}
#'   \item{V2}{Source database name or `"self"` for self-mappings.}
#'   \item{V3}{Cleaned mapped identifiers.}
#' }
form.uniprot.map <- function(unip.file) {
  unip.map <- fread(unip.file, header = F)

  unip.map[, V3 := sub("\\.\\d+$", "", V3)]

  unip.map[, V1 := sub("\\-\\d+$", "", V1)]

  # add V1 values to V3 for good measure

  self.map <- unique(unip.map[, .(V1 = V1, V2 = "self", V3 = V1)])

  unip.map <- rbind(unip.map, self.map)

  unip.map <- unip.map[!duplicated(paste(V1, V3))]

  unip.map
}


#' Capitalizes words in a sentence
#'
#' Adapted from the Examples in \code{toupper}
#'
capwords <- function(s, strict = FALSE) {
  cap <- function(s) {
    paste(toupper(substring(s, 1, 1)),
      {
        s <- substring(s, 2)
        if (strict) tolower(s) else s
      },
      sep = "",
      collapse = " "
    )
  }
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

#' Safely Apply Title Case Transformation
#'
#' A wrapper around \code{capwords} that gracefully handles encoding issues
#' and malformed strings. When a warning is thrown (e.g., "unable to translate")
#' or the result contains "NANA" due to encoding failures, the original string
#' is returned instead.  Note this function was written by Anthropic Claude Opus 4.5.
#'
#' @param s A character vector of strings to transform.
#'
#' @return A character vector the same length as \code{s}, with title case
#'   applied where possible. Strings that trigger warnings or produce invalid
#'   output are returned unchanged.
#'
safe_capwords <- function(s) {
  vapply(s, function(x) {
    tryCatch(
      {
        res <- capwords(x, strict = TRUE)
        if (is.na(res) || grepl("NANA", res, fixed = TRUE)) x else res
      },
      warning = function(w) x,
      error = function(e) x
    )
  }, character(1), USE.NAMES = FALSE)
}
