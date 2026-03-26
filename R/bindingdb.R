#' Form BindingDB into a compatible database
#'
#' @param unip.map
#'   A \code{data.table} describing a mapping between UniProt IDs and
#'   another protein identifier space, as in
#'   \code{\link{add.sorger.to.bioassay}}. The function expects at least
#'   two columns (typically named \code{V1} and \code{V3} when read),
#'   where:
#'   \itemize{
#'     \item \code{V3} contains UniProt IDs present in BindingDB.
#'     \item \code{V1} contains the "canonical" UniProt IDs to be used
#'           in the final targetome.
#'   }
#'
#' @param bdb.file
#'   Path to the BindingDB full export file, such as
#'   \code{"BindingDB_All_202311.tsv"}, downloadable from the BindingDB
#'   database (\url{https://www.bindingdb.org/bind/downloads/}).
#'
#' @details
#' The function performs the following major steps:
#' \enumerate{
#'   \item Read the BindingDB export via \code{readr::read_delim}.
#'   \item Restrict to human targets
#'   \item Split the data into:
#'         \itemize{
#'           \item \emph{Single-chain targets}: rows with
#'                 \code{Number of Protein Chains in Target} equal to 1.
#'                 For these, derive \code{uniprot_id} by preferring the
#'                 SwissProt ID and falling back to the TrEMBL ID when
#'                 SwissProt is missing.
#'           \item \emph{Multichain targets}: rows with
#'                 \code{Number of Protein Chains in Target} greater than 1.
#'                 For these, the multiple "Primary ID" columns are
#'                 gathered into long format, parsed to extract the ID
#'                 type (SwissProt vs. TrEMBL) and chain index, then
#'                 reshaped back to wide format to select one UniProt ID
#'                 per chain (preferring SwissProt when available).
#'         }
#'   \item Combine single-chain and multichain subsets into a unified
#'         BindingDB table with one row per compound–target–reference,
#'         selecting only:
#'         \code{pubchem_cid}, \code{uniprot_id}, \code{pubmed_id}
#'         and the affinity columns (\code{Ki (nM)}, \code{IC50 (nM)},
#'         \code{Kd (nM)}, \code{EC50 (nM)}).
#'   \item Reshape the affinity columns into long form using
#'         \code{tidyr::pivot_longer}, yielding one row per
#'         compound–target–reference–assay:
#'         \itemize{
#'           \item \code{assay_type}: Derived from the column name
#'                 (e.g. \code{"Ki"}, \code{"IC50"}).
#'           \item \code{value}: The original character value from the
#'                 affinity column.
#'         }
#'   \item Parse \code{value} to extract:
#'         \itemize{
#'           \item \code{assay_relation}: Comparison operator (e.g.
#'                 \code{"<"}, \code{">"}, or \code{"="}), inferred via
#'                 regular expressions and normalized by removing
#'                 whitespace.
#'           \item \code{assay_value}: Numeric portion of the value,
#'                 converted to \code{numeric}.
#'         }
#'         Only rows with non-missing \code{pubmed_id},
#'         \code{assay_relation}, and \code{assay_value} are kept.
#'   \item Some rows contain multiple UniProt IDs in a single field
#'         (e.g. \code{"Q96CA5 P98170"}). These are handled by
#'         \code{tidyr::separate_longer_delim}, which splits
#'         \code{uniprot_id} on whitespace into separate rows, yielding
#'         one UniProt ID per row.
#'   \item Convert the resulting data to a \code{data.table} and
#'         restrict to rows with non-missing \code{pubchem_cid}.
#'   \item Harmonize UniProt IDs by merging with \code{unip.map}
#'   \item Remove exact duplicate assays by building an \emph{assay key}
#'         from \code{pubchem_cid}, \code{pubmed_id}, \code{uniprot_id},
#'         \code{assay_type}, \code{assay_relation}, and a rounded
#'         (7-decimal) version of \code{assay_value}, and keeping only
#'         the first occurrence of each key.
#' }
#'
#' @return A `data.table` containing the BindingDB database
#'
form.bindingdb.db <- function(unip.map, bdb.file) {
  # note parsing warnings are expected here, *seem* to be ok
  bdb <- readr::read_delim(bdb.file)

  human.bdb <- filter(bdb, `Target Source Organism According to Curator or DataSource` == "Homo sapiens")

  # first reconcile the numerous protein IDs

  human.single <- filter(human.bdb, `Number of Protein Chains in Target (>1 implies a multichain complex)` == 1)

  human.single <- mutate(human.single, uniprot_id = ifelse(is.na(`UniProt (SwissProt) Primary ID of Target Chain 1`), `UniProt (TrEMBL) Primary ID of Target Chain 1`, `UniProt (SwissProt) Primary ID of Target Chain 1`))

  human.multi <- filter(human.bdb, `Number of Protein Chains in Target (>1 implies a multichain complex)` > 1)

  human.bdb.prot <- pivot_longer(select(human.multi, `BindingDB Reactant_set_id`, `Target Name`, contains("Primary ID")), cols = contains("Primary ID"))

  human.bdb.prot <- mutate(human.bdb.prot,
    id_type = str_match(name, "\\((.+)\\)")[, 2],
    index = str_match(name, "...(\\d+)")[, 2]
  )

  human.bdb.wide <- pivot_wider(human.bdb.prot, id_cols = c(`BindingDB Reactant_set_id`, `Target Name`, index), names_from = id_type, values_from = value)

  human.bdb.wide <- filter(human.bdb.wide, (is.na(SwissProt) == F) | (is.na(TrEMBL) == F))

  human.bdb.wide <- mutate(human.bdb.wide, uniprot_id = ifelse(is.na(SwissProt), TrEMBL, SwissProt))

  human.multi <- inner_join(human.multi, select(human.bdb.wide, `BindingDB Reactant_set_id`, uniprot_id), by = "BindingDB Reactant_set_id")

  comb.bdb <- bind_rows(human.single, human.multi[, names(human.single)])
  
  comb.bdb <- select(comb.bdb, pubchem_cid = `PubChem CID`, uniprot_id, pubmed_id = PMID, `Ki (nM)`, `IC50 (nM)`, `Kd (nM)`, `EC50 (nM)`)

  # fix issue with pmid being converted to 1e7 and then into a character
  
  comb.bdb <- mutate(comb.bdb, pubmed_id = as.character(as.integer(pubmed_id)))
  
  # next reconcile the assay values

  comb.bdb.assay <- pivot_longer(data = comb.bdb, cols = `Ki (nM)`:`EC50 (nM)`, names_to = "assay_type")

  # logic to pull out > < etc otherwise = and convert value to numeric

  comb.bdb.assay <- mutate(comb.bdb.assay,
    assay_relation = str_match(value, "^(\\D+)[\\d\\.]+$")[, 2],
    assay_relation = sub(" ", "", ifelse(assay_relation == " ", "=", assay_relation)),
    assay_value = as.numeric(str_match(value, "^\\D+([\\d\\.]+)$")[, 2])
  )

  comb.bdb.assay <- mutate(comb.bdb.assay, assay_type = sub(" (nM)", "", assay_type, fixed = T))

  comb.bdb.assay <- filter(comb.bdb.assay, is.na(pubmed_id) == F)

  # keep one assay value where available

  full.triple <- unique(select(comb.bdb.assay, pubchem_cid, uniprot_id, pubmed_id))

  non.na.values <- filter(comb.bdb.assay, is.na(assay_relation) == F & is.na(assay_value) == F)

  full.bdb <- left_join(full.triple, non.na.values, by = c("pubchem_cid", "uniprot_id", "pubmed_id"))

  # need to fix these: 'Q96CA5 P98170'

  multi.unips <- full.bdb %>% filter(grepl(" ", uniprot_id))

  fixed.unips <- separate_longer_delim(multi.unips, cols = uniprot_id, delim = " ")

  full.bdb <- filter(full.bdb, uniprot_id %in% multi.unips$uniprot_id == F)

  full.bdb <- bind_rows(
    full.bdb,
    fixed.unips
  )

  bdb.dt <- as.data.table(full.bdb)

  bdb.dt <- bdb.dt[is.na(pubchem_cid) == F]

  # check against uniprot ref

  bdb.dt <- merge(bdb.dt, unip.map[, .(new_uniprot_id = V1, uniprot_id = V3)], by = "uniprot_id", all.x = T, all.y = F)

  # The remaining NAs appear to be non-human...
  # bdb.dt[is.na(new_uniprot_id),.N,by=uniprot_id]

  bdb.dt <- bdb.dt[is.na(new_uniprot_id) == F]

  stopifnot(bdb.dt[uniprot_id != new_uniprot_id, .N] == 0)

  # so can just keep uniprot_id for the rest

  # remove duplicates

  bdb.dt[, assay_key := paste(pubchem_cid, pubmed_id, uniprot_id, assay_type, assay_relation, as.character(round(assay_value, digits = 7)))]

  bdb.dt <- bdb.dt[!duplicated(assay_key)]

  bdb.dt[, assay_key := NULL]

  bdb.dt
}

#' Merge the BindingDB database into Targetome
#'
#' @param targetome
#'   A \code{data.table} containing Targetome, usually the
#'   return value of \code{\link{add.sorger.to.targetome}}. It is expected
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
#' @param bdb.dt A \code{data.table} containing the BindingDB database containing
#'   containing at least columns for:
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
#' A \code{data.table} equal to the input \code{targetome} with
#' additional BindingDB-derived data appended. The resulting
#' table preserves the existing schema of \code{targetome}.
#'
add.bindingdb.to.targetome <- function(targetome, bdb.dt) {
  bdb.dt <- merge(bdb.dt, targetome[, .(pubchem_cid, uniprot_id, pubmed_id, in_tome = T)], by = c("pubchem_cid", "uniprot_id", "pubmed_id"), all.x = T, all.y = F)

  # from previous work, remove data from reference 18183025
  ## e.g. Imatinib and ABL1 has a value > 10,000 but doesn't seem to be represented in resource
  
  bdb.dt <- bdb.dt[pubmed_id != "18183025"]

  bdb.dt$database <- "bindingdb"

  targetome <- rbind(
    targetome,
    bdb.dt[(is.na(in_tome) == T) & (is.na(assay_value) == F), names(targetome), with = F]
  )

  targetome
}
