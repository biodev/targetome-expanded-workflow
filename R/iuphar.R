#' Form IUPHAR/GtoPdb interactions into a compatible database
#'
#' @param iuphar.inter.file
#'   Path to the IUPHAR/GtoPdb interactions file, usually the
#'   \code{"interactions.tsv"} file obtained from
#'   \url{https://www.guidetopharmacology.org/DATA/interactions.tsv}.
#'
#' @param iuphar.lig.file
#'   Path to the IUPHAR/GtoPdb ligand table, usually the
#'   \code{"ligands.tsv"} file obtained from
#'   \url{https://www.guidetopharmacology.org/DATA/ligands.tsv}. T
#'
#' @param unip.map
#'   A \code{data.table} describing a mapping between UniProt IDs and
#'   canonical UniProt identifiers.  It should contain
#'   at least two columns (typically named \code{V1} and \code{V3} when
#'   read), where:
#'   \itemize{
#'     \item \code{V3} holds UniProt IDs present in the IUPHAR data.
#'     \item \code{V1} holds the canonical UniProt IDs to be used in the
#'           final \code{targetome}.
#'   }
#'
#' @details
#' The function performs the following main steps:
#' \enumerate{
#'
#'   \item Read the IUPHAR interaction and ligand table and merge
#'   \item Clean affinities and identifiers:
#'         \itemize{
#'           \item Drop rows with empty \code{uniprot_id} or
#'                 \code{pubmed_id}.
#'         }
#'   \item Handle multiple UniProt IDs or PubMed IDs delimited by
#'         \code{"|"} by converting to separate rows.
#'   \item Harmonize UniProt IDs using \code{unip.map}:
#'   \item Remove exact duplicate assays
#' }
#'
#' @return A `data.table` containing the IUPHAR database
#'
form.iuphar.db <- function(iuphar.inter.file, iuphar.lig.file, unip.map) {
  iuphar <- fread(iuphar.inter.file)

  human.iuphar <- iuphar[
    `Target Species` == "Human" & is.na(`Ligand PubChem SID`) == F,
    .(`Ligand ID`,
      action = Type,
      uniprot_id = `Target UniProt ID`,
      pubmed_id = `PubMed ID`,
      assay_relation = `Original Affinity Relation`,
      assay_type = `Original Affinity Units`,
      assay_value = ifelse(is.na(`Original Affinity Median nm`),
        mapply(function(x, y) mean(c(x, y), na.rm = T), `Original Affinity Low nm`, `Original Affinity High nm`),
        `Original Affinity Median nm`
      )
    )
  ]

  # add in CID annotation

  ligs <- fread(iuphar.lig.file)

  human.iuphar <- merge(ligs[, .(`Ligand ID`, pubchem_cid = `PubChem CID`)], human.iuphar, by = "Ligand ID", all.x = F, all.y = T)

  human.iuphar <- human.iuphar[is.na(pubchem_cid) == F]

  human.iuphar[is.nan(assay_value), assay_value := NA_real_]

  human.iuphar <- human.iuphar[uniprot_id != "" & pubmed_id != ""]

  # need to fix uniprot_id and pubmed delimited by |

  hi.tbl <- tidyr::separate_longer_delim(tibble::as_tibble(human.iuphar), cols = c("uniprot_id"), delim = "|")

  hi.tbl <- tidyr::separate_longer_delim(hi.tbl, cols = c("pubmed_id"), delim = "|")

  hi.dt <- as.data.table(hi.tbl)

  # reconcile with uniprot map

  hi.dt <- merge(hi.dt, unip.map[, .(new_uniprot_id = V1, uniprot_id = V3)], by = "uniprot_id", all.x = T, all.y = F)

  # The remaining NAs appear to be non-human...
  # hi.dt[is.na(new_uniprot_id),.N,by=uniprot_id]
  #   uniprot_id N
  # 1:   P56856-2 1

  # manually add in...

  hi.dt[uniprot_id == "P56856-2", new_uniprot_id := "P56856"]

  stopifnot(hi.dt[uniprot_id != new_uniprot_id, .N] == 1)

  hi.dt[, uniprot_id := new_uniprot_id]

  # remove duplicates

  hi.dt[, assay_key := paste(pubchem_cid, pubmed_id, uniprot_id, assay_type, assay_relation, as.character(round(assay_value, digits = 7)))]

  hi.dt <- hi.dt[!duplicated(assay_key)]

  hi.dt
}

#' Merge the IUPHAR database into Targetome
#'
#' @param targetome
#'   A \code{data.table} containing Targetome, usually the
#'   return value of \code{\link{add.bindingdb.to.targetome}}. It is expected
#'   to have a harmonized schema containing, at minimum, the columns:
#'   \itemize{
#'     \item \code{database}: Data source label
#'           (e.g. \code{"pubchem_bioassay"}, \code{"sorger_sms"}).
#'     \item \code{pubchem_cid}: PubChem compound identifier.
#'     \item \code{pubmed_id}: PubMed identifier associated with the assay.
#'     \item \code{uniprot_id}: UniProt target identifier.
#'     \item \code{assay_type}: Activity type (e.g. Ki, IC50, etc.).
#'     \item \code{assay_relation}: Qualifier for the value
#'           (e.g. \code{"="}, \code{"<"}, \code{">"}).
#'     \item \code{assay_value}: Numeric activity value.
#'   }
#'
#' @param hi.dt A \code{data.table} containing the IUPHAR database with
#'   at least columns for:
#'   \itemize{
#'     \item \code{pubchem_cid}: PubChem compound identifier.
#'     \item \code{pubmed_id}: PubMed identifier associated with the assay.
#'     \item \code{uniprot_id}: UniProt target identifier.
#'     \item \code{assay_type}: Activity type (e.g. Ki, IC50, etc.).
#'     \item \code{assay_relation}: Qualifier for the value
#'           (e.g. \code{"="}, \code{"<"}, \code{">"}).
#'     \item \code{assay_value}: Numeric activity value.
#'     }
#'
#' @return
#'   A \code{data.table} equal to the input \code{targetome} with
#'   additional IUPHAR-derived data appended. The resulting
#'   table preserves the existing schema of \code{targetome}.
#'
add.iuphar.to.targetome <- function(targetome, hi.dt) {
  hi.dt <- merge(hi.dt, unique(targetome[, .(pubchem_cid, uniprot_id, pubmed_id, in_tome = T)]), by = c("pubchem_cid", "uniprot_id", "pubmed_id"), all.x = T, all.y = F)

  hi.dt$database <- "iuphar"

  # The existing annotation will be carried along with iuphar

  targetome <- rbind(
    targetome,
    hi.dt[(is.na(in_tome) == T) & (is.na(assay_value) == F), names(targetome), with = F]
  )

  targetome
}
