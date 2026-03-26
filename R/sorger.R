#' Map Sorger compounds to unique PubChem CIDs
#'
#' Given Sorger LSP compound annotation files and PubChem CID mapping
#' tables, construct a one-to-one mapping between Sorger compounds
#' (identified by InChI) and PubChem compound identifiers (CIDs),
#' restricted to compounds that appear in the Sorger biochemical
#' activity table. The function resolves multiple CID matches per InChI
#' using a series of heuristics based on preferred names and PubChem
#' parent/child relationships.
#'
#' @param sorger.comp.dict.file
#'   Path to the Sorger compound dictionary file
#'   (e.g. `"sorger/lsp_compound_dictionary.csv.gz"`). This file is
#'   expected to be a comma-separated table readable by
#'   \code{data.table::fread} and to contain at least the columns
#'   \code{inchi}, \code{lspci_id}, and \code{pref_name}.
#'
#' @param sorger.biochem.file
#'   Path to the Sorger biochemical activity file
#'   (e.g. `"sorger/lsp_biochem.csv.gz"`). This file is expected to
#'   contain biochemical assay information with a column
#'   \code{lspci_id} that links back to the compound dictionary.
#'
#' @param pubchem.cid.inchi.file
#'   Path to the PubChem CID–InChI mapping file
#'   (e.g. `"pubchem/CID-InChI-Key.gz"`). This file is read with
#'   \code{data.table::fread} and is expected to have at least two
#'   columns, where the first column (\code{V1}) is the PubChem CID and
#'   the second column (\code{V2}) contains the InChI
#'   string used to merge with the Sorger compound dictionary.
#'
#' @param pubchem.cid.title.file
#'   Path to the PubChem CID–title mapping file
#'   (e.g. `"pubchem/CID-Title.gz"`). This file is read with
#'   \code{data.table::fread} (with \code{quote = ""}) and is expected
#'   to have at least two columns, where the first column (\code{V1})
#'   is the PubChem CID and the second column (\code{V2}) is a
#'   human-readable PubChem compound name or description.
#'
#' @param pubchem.cid.parent.file
#'   Path to the PubChem parent–child CID mapping file
#'   (e.g. `"pubchem/CID-Parent.gz"`). This file is read with
#'   \code{data.table::fread} and is expected to have at least two
#'   columns, where the first column (\code{V1}) is the  CID and
#'   the second column (\code{V2}) is the parent CID.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads the Sorger compound dictionary and PubChem CID–InChI
#'         mapping, and merges them by InChI to assign PubChem CIDs
#'         to Sorger compounds.
#'   \item Reads the PubChem CID–title file and merges it in to add a
#'         \code{pubchem_desc} field.
#'   \item Filters the merged table to rows with non-missing CIDs and
#'         restricts to compounds that are present in the Sorger
#'         biochemical activity table.
#'   \item For each InChI with multiple matched CIDs, resolves
#'         duplicates using the following prioritized heuristics:
#'         \enumerate{
#'           \item If exactly one CID has a PubChem title matching the
#'                 Sorger \code{pref_name} (case-insensitive), keep
#'                 that mapping.
#'           \item Otherwise, if parent–child CID relationships are
#'                 available, prefer the CID that has a unique parent
#'                 according to \code{pubchem.cid.parent.file}.
#'         }
#'   \item Asserts that each InChI maps to exactly one CID in the final
#'         table (via \code{stopifnot}), and returns the resulting
#'         data.
#'
#' @return
#' A \code{data.table} containing a one-to-one mapping between Sorger
#' compounds and PubChem CIDs for all compounds that appear in the
#' biochemical activity file. The returned table includes at least:
#' \itemize{
#'   \item \code{lspci_id}: Sorger compound identifier.
#'   \item \code{inchi}: InChI string used for matching.
#'   \item \code{CID}: resolved PubChem compound identifier.
#'   \item \code{pref_name}: Sorger preferred compound name.
#'   \item \code{pubchem_desc}: PubChem compound title/description,
#'         when available.
#'   \item Any additional columns carried through from the input
#'         Sorger tables.
#' }
make.sorger.pubchem.map <- function(sorger.comp.dict.file, sorger.biochem.file,
                                    pubchem.cid.inchi.file, pubchem.cid.title.file,
                                    pubchem.cid.parent.file) {
  dir.create("tmp")

  sorger.drugs <- fread(sorger.comp.dict.file, sep = ",", tmpdir = "tmp")

  inchi.to.cid <- fread(pubchem.cid.inchi.file, tmpdir = "tmp")

  sorger.cid <- merge(sorger.drugs, inchi.to.cid, by.x = "inchi", by.y = "V2", all.x = T, all.y = F)

  sorger.cid[, CID := V1]

  pubchem.title <- fread(pubchem.cid.title.file, tmpdir = "tmp", quote = "")

  # merge

  sorger.cid <- merge(sorger.cid, pubchem.title[, .(CID = V1, pubchem_desc = V2)], by = "CID", all.x = T, all.y = F)

  sorger.cid[, `:=`(V1 = NULL, V3 = NULL)]

  # limit these to those with matching CIDs and merge with biochem below

  sorger.cid <- sorger.cid[is.na(CID) == F]

  sorger.biochem <- fread(sorger.biochem.file, tmpdir = "tmp")

  # try to resolve duplicate inchi/CID entries

  ## first limit to only those CIDs covered in the biochem entries

  sorger.cid.in.biochem <- sorger.cid[lspci_id %in% sorger.biochem$lspci_id]

  sorger.cid.in.biochem[, inchi_count := .N, by = inchi]

  duplicate.inchi <- sorger.cid.in.biochem[inchi_count > 1]

  ## find those that match the `pref_name` of sorger

  best.matches <- duplicate.inchi[tolower(pref_name) == tolower(pubchem_desc)]

  stopifnot(best.matches[, .N, by = inchi][, all(N == 1)])

  ## now remove from the duplicate set

  duplicate.inchi <- duplicate.inchi[inchi %in% best.matches$inchi == F]

  ## limit to only parents

  cid.parents <- fread(pubchem.cid.parent.file, tmpdir = "tmp")

  norm.to.parents <- cid.parents[V1 %in% duplicate.inchi$CID]

  parent.resolved <- duplicate.inchi[CID %in% norm.to.parents$V2]

  parent.resolved[, inchi_count := .N, by = inchi]

  parent.resolved <- parent.resolved[inchi_count == 1]

  sorger.cid.in.biochem <- sorger.cid.in.biochem[inchi_count == 1]

  sorger.cid.in.biochem <- rbind(
    sorger.cid.in.biochem,
    best.matches[, names(sorger.cid.in.biochem), with = F],
    parent.resolved[, names(sorger.cid.in.biochem), with = F]
  )

  sorger.cid.in.biochem[, inchi_count := .N, by = inchi]

  stopifnot(sorger.cid.in.biochem[, all(inchi_count == 1)])

  sorger.cid.in.biochem
}


#' Form Sorger SMS biochemical data into a compatible database
#'
#' Combine biochemical assay data from the Sorger Small Molecule
#' Suite (SMS) with the Sorger–PubChem CID mapping
#' (from the output of \code{\link{make.sorger.pubchem.map}}),
#' Sorger assay/target/reference tables, and a UniProt ID mapping.
#'
#' @param sorger.cid.in.biochem
#'   A \code{data.table} mapping Sorger compounds to PubChem CIDs,
#'   typically the output of \code{\link{make.sorger.pubchem.map}}.
#'
#' @param sorger.biochem.file
#'   Path to the Sorger biochemical activity file
#'   (e.g. `"sorger/lsp_biochem.csv.gz"`).
#'
#' @param sorger.target.file
#'   Path to the Sorger target mapping file
#'   (e.g. `"sorger/lsp_target_mapping.csv.gz"`).
#'
#' @param sorger.ref.file
#'   Path to the Sorger reference mapping file
#'   (e.g. `"sorger/lsp_references.csv.gz"`).
#'
#' @param unip.map
#'   A \code{data.table} describing a mapping between UniProt IDs and
#'   another protein identifier space. The function expects at least
#'   two columns (named \code{V1} and \code{V3} when read) such that
#'   \code{V3} corresponds to the UniProt IDs present in
#'   \code{sorger.target.file} and \code{V1} corresponds to the
#'   "canonical" UniProt IDs used in the targetome.
#'
#' @details
#' The function performs the following main steps:
#' \enumerate{
#'   \item Read the Sorger biochemical activity table, filter out rows
#'         with missing \code{value_unit}, and select the core assay
#'         columns. Merge with \code{sorger.cid.in.biochem} to attach
#'         PubChem CIDs for each Sorger compound.
#'   \item Read the Sorger target mapping, merge in UniProt IDs via
#'         \code{lspci_target_id}, and thus link each assay to a
#'         UniProt target.
#'   \item Read the Sorger references table, filter to
#'         \code{reference_type == "pubmed_id"}, and merge in the
#'         corresponding PubMed IDs.
#'   \item Harmonize UniProt IDs by merging with \code{unip.map}, then
#'         drop rows without a mapped \code{new_uniprot_id}. A
#'         \code{stopifnot} check ensures that original and mapped
#'         UniProt IDs agree (\code{uniprot_id == new_uniprot_id})
#'         for the retained rows, after which the original
#'         \code{uniprot_id} is used.
#'   \item Construct an \emph{assay key} for each row as a concatenation
#'         of \code{pubchem_cid}, \code{pubmed_id}, \code{uniprot_id},
#'         \code{assay_type}, \code{assay_relation}, and a rounded
#'         version of \code{assay_value} (7 decimal places), and remove
#'         exact duplicates based on this key.
#' }
#'
#' @return A `data.table` containing the Sorger SMS database
#'
form.sorger.db <- function(sorger.cid.in.biochem, sorger.biochem.file,
                           sorger.target.file, sorger.ref.file,
                           unip.map) {
  if (dir.exists("tmp")) {
    unlink("tmp", recursive = T)
  }

  sorger.biochem <- fread(sorger.biochem.file, tmpdir = "tmp")

  sorger.biochem <- merge(
    sorger.biochem[is.na(value_unit) == F, .(lspci_id, lspci_target_id, reference_id,
      assay_value = value, assay_type = value_type,
      assay_relation = value_relation
    )],
    sorger.cid.in.biochem[, .(lspci_id, pubchem_cid = CID)],
    by = "lspci_id"
  )

  target.map <- fread(sorger.target.file)

  sorger.biochem <- merge(sorger.biochem, unique(target.map[is.na(uniprot_id) == F, .(lspci_target_id, uniprot_id)]), by = "lspci_target_id")

  ref.map <- fread(sorger.ref.file)

  ref.map <- ref.map[reference_type == "pubmed_id"]

  sorger.biochem <- merge(sorger.biochem, ref.map[, .(reference_id, pubmed_id = reference_value)], by = "reference_id")

  sorger.biochem <- merge(sorger.biochem, unip.map[, .(uniprot_id = V3, new_uniprot_id = V1)], by = "uniprot_id", all.x = T, all.y = F)

  # can remove the nas as they appear to be mouse and rat references

  sorger.biochem <- sorger.biochem[is.na(new_uniprot_id) == F]

  stopifnot(sorger.biochem[uniprot_id != new_uniprot_id, .N] == 0)

  # can simply keep uniprot_id

  # remove duplicates

  sorger.biochem[, assay_key := paste(pubchem_cid, pubmed_id, uniprot_id, assay_type, assay_relation, as.character(round(assay_value, digits = 7)))]

  sorger.biochem <- sorger.biochem[!duplicated(assay_key)]

  sorger.biochem
}

#' Merge Sorger SMS biochemical data with the Targetome database
#'
#' @param targetome
#'   A \code{data.table} of existing bioassay target information (the
#'   "targetome"), generally derived from PubChem BioAssay. Must
#'   contain at least the columns:
#'   \itemize{
#'     \item \code{CID}: PubChem compound identifier.
#'     \item \code{uniprot_id}: UniProt target identifier.
#'     \item \code{PMID}: PubMed identifier for the assay.
#'     \item \code{Activity Name}: Assay type.
#'     \item \code{Activity Qualifier}: Relation symbol for the
#'           activity value (e.g. \code{"="}, \code{"<"},  \code{">"}).
#'     \item \code{Activity Value}: Numeric activity value.
#'   }
#'
#' @param sorger.biochem A \code{data.table} containg the Sorger SMS database.
#'   containing at least columns for:
#'   \itemize{
#'     \item \code{pubchem_cid}: PubChem compound identifier.
#'     \item \code{uniprot_id}: UniProt target identifier.
#'     \item \code{pubmed_id}: PubMed identifier for the assay.
#'     \item \code{assay_type}: Assay type.
#'     \item \code{assay_relation}: Relation symbol for the activity value.
#'     \item \code{assay_value}: Numeric activity value.
#'
#' @return A \code{data.table} combining the original \code{targetome} entries
#' (labeled with \code{database = "pubchem_bioassay"}) and new,
#' non-duplicated Sorger SMS biochemical assays
#' (\code{database = "sorger_sms"})
#'
add.sorger.to.targetome <- function(targetome, sorger.biochem) {
  sorger.biochem <- merge(sorger.biochem, targetome[, .(pubchem_cid = CID, uniprot_id, pubmed_id = as.character(PMID), in_tome = T)], by = c("pubchem_cid", "uniprot_id", "pubmed_id"), all.x = T, all.y = F)

  sorger.biochem$database <- "sorger_sms"

  targetome <- rbind(
    targetome[, .(database = "pubchem_bioassay", pubchem_cid = CID, pubmed_id = as.character(PMID), uniprot_id, assay_type = `Activity Name`, assay_relation = `Activity Qualifier`, assay_value = `Activity Value`)],
    sorger.biochem[(is.na(in_tome) == T) & (is.na(assay_value) == F), .(database, pubchem_cid, pubmed_id, uniprot_id, assay_type, assay_relation, assay_value)]
  )

  targetome
}
