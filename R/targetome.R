#' Construct an enriched drug–targetome combining assay data, classification & annotation
#'
#' Build a “first-pass” enriched targetome for drugs and their targets by
#' combining a pre-existing targetome (Drugs <-> UniProt targets <-> Pubmed IDs <-> Assays)
#' with multiple annotation sources covering classification, synonyms, synonyms,
#' parent-child compound groups, MeSH terms, DrugBank metadata, IUPHAR ligand
#' metadata, PubChem titles, and ATC codes.
#'
#' @param targetome
#'   A \code{data.table} representing the unified compound–target–assay
#'   dataset.  Required to contain at least the columns:
#'   \itemize{
#'     \item \code{pubchem_cid} — PubChem compound identifier.
#'     \item \code{uniprot_id} — UniProt target identifier.
#'     \item \code{pubmed_id} — PubMed reference identifier.
#'     \item \code{assay_value} — Numeric assay/affinity value.
#'     \item \code{database}, \code{assay_type}, \code{assay_relation} — assay metadata.
#'   }
#'
#' @param cid.sum
#'   A \code{data.table} (or tibble) produced by \code{\link{combine.mesh.pubchem}},
#'   summarizing PubChem CIDs into MeSH terms and MeSH pharmacological-action classes.
#'   Must contain at least:
#'   \itemize{
#'     \item \code{CID} — numeric PubChem CID.
#'     \item \code{mesh} — pipe-delimited character string of MeSH terms.
#'     \item \code{pharm_mesh} — pipe-delimited character string of MeSH pharmacological actions.
#'   }
#'
#' @param pubchem.syn.file
#'   Path to the PubChem synonym file (e.g. \code{"CID-Synonym-filtered.gz"})
#'   from the PubChem Compound/Extras FTP directory. This file provides
#'   synonyms for each CID and is used to annotate drug names and alternate names.
#'
#' @param cid.par.file
#'   Path to the PubChem CID-Parent file (e.g. \code{"CID-Parent.gz"})
#'   from the PubChem Compound/Extras FTP directory. Used to map CIDs to
#'   their “parent” CID (e.g. canonical compound forms) and collapse assay
#'   values at the parent level.
#'
#' @param cid.title.file
#'   Path to the PubChem CID title file (e.g. \code{"CID-Title.gz"})
#'   from the PubChem Compound/Extras FTP directory. Provides a primary “drug_name”
#'   for each CID.
#'
#' @param pubchem.atc.file
#'   A CSV file derived from the WHO ATC classification system applied to
#'   PubChem CIDs (e.g. mapping each CID to ATC codes).
#'
#' @details
#' The function proceeds in major stages:
#' \enumerate{
#'   \item Read all the files.
#'   \item Merge the CID-parent mapping with the input \code{targetome} to convert each row from \code{pubchem_cid} to its \code{parent_cid} and then collapse assay values:
#'         \itemize{
#'           \item Group by \code{parent_cid, uniprot_id, pubmed_id, database, assay_type, assay_relation} and summarize by the median of unique assay values
#'         }
#'   \item Extract unique drug CIDs and annotate them by:
#'         \item Merge with CID titles to get a \code{drug_name}.
#'         \item Merge with DrugBank metadata (filtered to non-missing \code{pubchem_cid}) to add \code{drugbank_name, drugbank_status, drugbank_type}.
#'         \item Process the IUPHAR ligand table by pivoting longer on \code{Approved, Withdrawn} columns, converting to a data.table, summarizing by
#'          PubChem CID concatenating where necessary (e.g. “approved/withdrawn”).
#'         \item Merge the processed IUPHAR ligand metadata into the drug annotation table.
#'         \item Merge PubChem ATC classification.
#'         \item Merge MeSH annotations for the same CIDs to annotate \code{mesh} and \code{pharm_mesh}.
#'         \item Read PubChem synonyms and collapse for relevant CIDs into a \code{synonyms} string (pipe-delimited) and merge into the drug annotation table.
#' }
#'   \item Return a list with:
#'         \item \code{drugs} — the drug annotation table described above.
#'         \item \code{targets} — the collapsed assay table, with one row per
#'           parent CID–UniProt–PubMed–database-assay_type–assay_relation combination and median assay_value.
#' }
#'
#' @return
#' A list with two elements:
#' \itemize{
#'   \item \code{drugs}: A \code{data.table} of unique drugs annotated with:
#'         \code{pubchem_cid}, \code{drug_name}, \code{drugbank_name}, \code{drugbank_status},
#'         \code{drugbank_type}, \code{iuphar_type}, \code{iuphar_status}, \code{atc},
#'         \code{mesh}, \code{pharm_mesh}, \code{synonyms}.
#'   \item \code{targets}: A \code{data.table} of unique parent-CID–target–reference assay rows with columns:
#'         \code{pubchem_cid} (parent), \code{uniprot_id}, \code{pubmed_id},
#'         \code{database}, \code{assay_type}, \code{assay_relation}, \code{assay_value}.
#' }
#'
form.initial.targetome <- function(targetome, cid.sum, pubchem.syn.file, cid.par.file, cid.title.file, pubchem.atc.file) {
  cid.par <- fread(cid.par.file)

  # this is no longer a param
  # ligs <- fread(iuphar.lig.file)

  cid.titles <- fread(cid.title.file, quote = "")

  # this is no longer a param
  # drb <- fread(drugbank.drug.file)

  pubchem.atc <- fread(pubchem.atc.file)

  # removing those without a parent for now, but can add them back in later (n=9672 entries, 25,668 interactions)

  tome.m <- merge(cid.par[, .(pubchem_cid = V1, parent_cid = V2)], targetome, by = "pubchem_cid", all = F)

  .median.of.unique <- function(x) {
    x.char <- as.character(round(x, digits = 7))

    x <- x[!duplicated(x.char)]

    median(x, na.rm = T)
  }

  # summarize by parent

  uniq.tome <- tome.m[, .(assay_value = .median.of.unique(assay_value)), by = .(pubchem_cid = parent_cid, uniprot_id, pubmed_id, database, assay_type, assay_relation)]

  tome.drugs <- unique(uniq.tome[, .(pubchem_cid)])

  # add in names (titles)

  tome.drugs <- merge(tome.drugs, cid.titles[, .(pubchem_cid = V1, drug_name = V2)], by = "pubchem_cid")

  # add in categorizations, first ATC

  ## annotation column appears to be only ATC, delimited by '|'

  tome.drugs <- merge(tome.drugs, pubchem.atc[, .(pubchem_cid = cid, atc = annotation)], all.x = T, all.y = F)

  ## add in mesh terms as well

  tome.drugs <- merge(tome.drugs, cid.sum, by.x = "pubchem_cid", by.y = "CID", all.x = T, all.y = F)

  ## finally add in synonyms
  
  syns <- fread(pubchem.syn.file, quote = "")
  
  syns.par <- merge(syns[,.(pubchem_cid = V1, synonym = V2)], cid.par[, .(pubchem_cid = V1, parent_cid = V2)], by="pubchem_cid")

  uniq.syns <- syns.par[,.(synonyms = paste(synonym, collapse = "|")), by = .(pubchem_cid = parent_cid)]

  tome.drugs.m <- merge(tome.drugs, uniq.syns, by = "pubchem_cid", all.x = T, all.y = F)

  stopifnot(tome.drugs[,.N] == tome.drugs.m[,.N])
  
  list(drugs = tome.drugs.m, targets = uniq.tome)
}

#' Add InChI / InChIKey annotations to the targetome based on PubChem IDs
#'
#' Extend an existing drug–targetome by joining PubChem InChI and
#' InChIKey strings to each compound (PubChem CID) on the drug side. The
#' function uses the PubChem “CID – InChIKey” flat file to map each
#' PubChem CID to its full InChI string and standard InChIKey, and then
#' merges this data into the \code{drugs} table in the targetome workinglist. The
#' \code{targets} part of the targetome is retained unmodified.
#'
#' @param targetome
#'   A \code{list} with at least two elements:
#'   \itemize{
#'     \item \code{drugs}: A \code{data.table} of drug annotation, containing at least column \code{pubchem_cid}.
#'     \item \code{targets}: A \code{data.table} of compound–target–reference assay lines,
#'     typically carrying \code{pubchem_cid}, \code{uniprot_id}, etc.
#'   }
#'
#' @param inchi.file
#'   Path to the PubChem “CID–InChIKey” mapping file
#'   (e.g. \code{"CID-InChI-Key.gz"}) from the PubChem FTP site
#'   (\url{https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/}). This file
#'   is expected to provide at least three columns: the first column (CID),
#'   a second column (InChI string),and a third column (InChIKey string).
#'   The function renames these to \code{pubchem_cid}, \code{inchi}, and
#'   \code{inchi_key} respectively.
#'
#'
#' @return
#' A \code{list} with two elements:
#' \itemize{
#'   \item \code{drugs}: A \code{data.table} annotated with columns \code{pubchem_cid},
#'         \code{inchi}, \code{inchi_key}, and all original drug‐annotation columns.
#'   \item \code{targets}: The original \code{targets} table from the input,
#'         unchanged.
#' }
#'
add.inchi <- function(targetome, inchi.file) {
  if (dir.exists("tmp") == F) {
    dir.create("tmp")
  }

  cid2inchi <- fread(inchi.file, tmpdir = "tmp")

  tome.drugs <- merge(cid2inchi[, .(pubchem_cid = V1, inchi = V2, inchi_key = V3)], targetome$drugs, by = "pubchem_cid", all.x = F, all.y = T)

  list(drugs = tome.drugs, targets = targetome$targets)
}

#' Add annotations from ChEMBL to the targetome
#'
#' Extend an existing targetome list of \code{data.table}'s by merging in ChEMBL molecule-level metadata
#' (e.g., ChEMBL ID, clinical phase, therapeutic flag, molecule type, withdrawn flag, chemical probe, natural product status)
#' based on InChIKey matching.
#'
#' @param targetome
#'   A \code{list} with at least two elements:
#'   \itemize{
#'     \item \code{drugs}: A \code{data.table} of drug annotation, containing at minimum \code{inchi_key}, plus \code{pubchem_cid} and other identifiers.
#'     \item \code{targets}: A \code{data.table} of compound–target–reference-assay rows (unchanged by this function).
#'   }
#'
#' @param chembl.file
#'   Path to the ChEMBL SQLite database file (e.g., downloaded from
#'   \url{https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_sqlite.tar.gz}).
#'   The schema includes tables \code{molecule_dictionary} and \code{compound_structures} (among others)
#'   in which the columns \code{standard_inchi_key}, \code{pref_name}, \code{chembl_id}, \code{therapeutic_flag},
#'   \code{molecule_type}, \code{withdrawn_flag}, \code{chemical_probe}, \code{natural_product}, and \code{max_phase}
#'   are available.
#'
#' @details
#' The function executes the following steps:
#' \enumerate{
#'   \item Connect to the ChEMBL SQLite database.
#'   \item Query and join the tables to assemble the drug dictionary metadata.
#'   \item Use a lookup table converting \code{max_phase} values to descriptive \code{clinical_phase} strings (for example, \code{4} → “Approved”, etc.).
#'   \item Merge the drug dictionary to the \code{targetome$drugs} data table by the \code{inchi_key} column (matching ChEMBL’s \code{standard_inchi_key}).
#'   \item Retain only drugs present in \code{targetome$drugs} (all rows preserved except those without matching metadata).
#'   \item Return the updated targetome as a list
#'   }
#'
#' @return
#' A list with two elements:
#' \itemize{
#'   \item \code{drugs}: A \code{data.table} of drug entities annotated with ChEMBL metadata merged onto the existing drug list from the targetome.
#'   \item \code{targets}: The original \code{targetome$targets} data table (unchanged).
#' }
#'
add.chembl.annots <- function(targetome, chembl.file) {
  con <- dbConnect(SQLite(), chembl.file)

  drug.dict <- data.table(dbGetQuery(con, "select * from molecule_dictionary join compound_structures using (MOLREGNO);"))

  dbDisconnect(con)

  # note parsed from https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html

  phase.dt <- data.table(
    max_phase = c(NA_real_, -1, .5, 1, 2, 3, 4),
    clinical_phase = c(
      "preclinical compounds with bioactivity data",
      "Clinical Phase unknown for drug or clinical candidate drug",
      "Early Phase 1 Clinical Trials",
      "Phase 1 Clinical Trials",
      "Phase 2 Clinical Trials",
      "Phase 3 Clinical Trials",
      "Approved"
    )
  )

  drug.dict <- merge(drug.dict, phase.dt, by = "max_phase")

  tome.drugs <- merge(drug.dict[, .(inchi_key = standard_inchi_key, pref_name, clinical_phase, chembl_id, therapeutic_flag, molecule_type, withdrawn_flag, chemical_probe, natural_product)], targetome$drugs, by = "inchi_key", all.x = F, all.y = T)

  list(drugs = tome.drugs, targets = targetome$targets)
}

#' Extract a synonym lookup table from the targetome
#'
#' Build a separate synonym table for drug entities in an existing
#' \code{targetome} \code{data.table} by expanding the pipe-delimited
#' \code{synonyms} field. The resulting synonym table is attached back into
#' the \code{targetome} list as a component, and the original \code{synonyms}
#' column in \code{targetome$drugs} is dropped.
#'
#' @param targetome
#'   A \code{list} (or object) that includes at least:
#'   \itemize{
#'     \item \code{drugs}: A \code{data.table} of drug entities, containing at minimum
#'           columns \code{pubchem_cid}, \code{inchi_key}, and \code{synonyms}.
#'     \item \code{targets}: A \code{data.table} of compound–target–reference-assay rows
#'           (this component is not modified by this function).
#'   }
#'
#' @details
#' The function carries out the following steps:
#' \enumerate{
#'   \item Convert \code{targetome$drugs} to a tibble for convenience.
#'   \item Expand the existing \code{synonyms} column (which may contain multiple
#'         synonyms delimited by “\(|\)”), producing one row per synonym.
#'   \item Add in entries for the standard drug names
#'   \item Add a column \code{lower_name}, which is the lowercase version of each
#'         \code{synonyms} value.
#'   \item Remove rows with missing \code{synonyms} and drop duplicates.
#'   \item Convert back to a \code{data.table}, assign it to \code{targetome$synonyms},
#'         and remove the original \code{synonyms} column from \code{targetome$drugs}.
#' }
#'
#' @return
#' A modified version of the input \code{targetome} (same structure: list with
#' \code{drugs} and \code{targets}), where:
#' \itemize{
#'   \item \code{targetome$drugs} no longer includes the original \code{synonyms} field.
#'   \item A new component \code{targetome$synonyms} (a \code{data.table}) is added,
#'         containing the columns:
#'         \code{pubchem_cid}, \code{inchi_key}, \code{synonyms}, \code{lower_name}.
#' }
#'
add.synonym.table <- function(targetome) {
  
  drug.tbl <- as_tibble(targetome$drugs)

  syn.tab <- separate_longer_delim(drug.tbl, cols = synonyms, delim = "|")

  # also add drug_name and pref_name to synonym table

  pref.tab <- filter(drug.tbl, is.na(pref_name) == F) %>% select(pubchem_cid, inchi_key, synonyms = pref_name)

  name.tab <- filter(drug.tbl, is.na(drug_name) == F) %>% select(pubchem_cid, inchi_key, synonyms = drug_name)
  
  all.syns <- bind_rows(
    select(syn.tab, pubchem_cid, inchi_key, synonyms),
    pref.tab,
    name.tab
  )

  all.syns <- mutate(all.syns, lower_name = sapply(synonyms, function(x) tryCatch(tolower(x), error = function(x) NA_character_)))

  all.syns <- filter(all.syns, is.na(synonyms) == F)
  
  all.syns <- unique(all.syns)

  targetome$synonyms <- as.data.table(all.syns)

  targetome$drugs$synonyms <- NULL

  targetome
}

#' Map drugs in the targetome to NCI Thesaurus concepts and populate annotations
#'
#' Integrate controlled-vocabulary annotations from the NCI Thesaurus (NCIt) into the targetome.
#' The function uses the NCIt flat file (Thesaurus.FLAT) to extract concept codes, synonyms, parent–child relationships, and display terms,
#' then links drug names/synonyms in the targetome to NCIt codes and traces parent paths to identify broader ontology placements (e.g., “Drug, Food, Chemical or Biomedical Material”).
#'
#' @param targetome
#'   A \code{list} containing at least two elements:
#'   \itemize{
#'     \item \code{drugs}: A \code{data.table} of drug entities, containing at minimum \code{pubchem_cid}, \code{inchi_key}, \code{synonyms}, and \code{lower_name}.
#'     \item \code{targets}: A \code{data.table} of compound–target–reference-assay rows
#'   }
#'   This object is typically produced by \code{\link{add.synonym.table}}.
#'
#' @param nci.file
#'   Path to the NCIt flat-file (e.g. \code{"Thesaurus.FLAT"} or included in \code{"Thesaurus.FLAT.zip"}) downloaded from the EVS FTP site
#'   (\url{https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus.FLAT.zip}) which contains NCIt concepts in a tab-delimited format.
#'
#' @details
#' The function executes the following steps:
#' \enumerate{
#'   \item Read the NCIt flat file.
#'   \item Rename columns to:
#'         \code{code}, \code{concept}, \code{parents}, \code{synonyms}, \code{definition},
#'         \code{display name}, \code{concept status}, \code{semantic type}, \code{concept in subset}.
#'   \item Compute \code{first_syn} (first synonym from \code{synonyms} delimited by “|”), and \code{use_name} as the display name if not empty, otherwise the first synonym.
#'   \item Expand the \code{parents} field to produce a two-column table \code{(child=code, parent=parents)} (converted to a \code{data.table}).
#'   \item Build a directed graph using \code{igraph::graph_from_data_frame()} where vertices are NCIt codes and edges are parent→child.
#'   \item Expand the \code{synonyms} field similarly to create a \code{syn.map} table of \code{(code, synonyms)}; add \code{lower_name = tolower(synonyms)}.
#'   \item Merge the \code{syn.map} with the targetome synonyms table (\code{targetome$synonyms$lower_name}) by \code{lower_name} to identify drugs with matched NCIt codes.
#'   \item For each matched NCIt \code{code}, compute *all shortest paths* in the graph from that code up to root(s). Then select a maximal set of distinct longest paths, convert each path into a string of “code – use_name” entries joined by “ > ”, discard the last element that corresponds to the drug itself (i.e., \code{code == x}). Glue multiple paths with “|”.
#'   \item Keep only those NCIt paths that stem from the broad code “C1908 – Drug, Food, Chemical or Biomedical Material” (via \code{grepl("1908", nci)==TRUE}).
#'   \item Build \code{cid.to.nci} as a table of \code{(pubchem_cid, nci)} where \code{nci} is the pipe-delimited set of path-strings, then merge this into \code{targetome$drugs} by \code{pubchem_cid}.
#'   \item Return the modified \code{targetome} \code{list}, with the \code{drugs} component now including a new column \code{nci}.
#' }
#'
#' @return
#' A modified version of the input \code{targetome} (same structure) in which:
#' \itemize{
#'   \item \code{targetome$drugs} now includes a new column \code{nci} (character) with pipe-delimited NCIt concept path annotations for each drug (or \code{NA} if none matched).
#'   \item \code{targetome$targets} remains unchanged.
#'   \item \code{targetome$synonyms} remains unchanged.
#' }
#'
add.nci.annots <- function(targetome, nci.file) {
  nci <- fread(nci.file, quote = "")

  names(nci) <- c("code", "concept", "parents", "synonyms", "definition", "display name", "concept status", "semantic type", "concept in subset")

  nci[, first_syn := tstrsplit(synonyms, "\\|", keep = 1)]

  nci[, use_name := ifelse(`display name` == "", first_syn, `display name`)]

  nci.pars <- as.data.table(tidyr::separate_longer_delim(data = dplyr::select(tibble::as_tibble(nci), code, parents), col = parents, delim = "|"))

  nci.graph <- graph_from_data_frame(nci.pars[, .(parents, code)])

  syn.map <- data.table(tidyr::separate_longer_delim(dplyr::select(tibble::as_tibble(nci), code, synonyms), cols = synonyms, delim = "|"))

  syn.map[, lower_name := tolower(synonyms)]

  found.drugs <- merge(unique(syn.map[, .(code, lower_name)]), unique(targetome$synonyms[, .(pubchem_cid, lower_name)]), by = "lower_name")

  drug.nci.terms <- rbindlist(lapply(unique(found.drugs$code), function(x) {
    tmp.sp <- all_shortest_paths(nci.graph, from = V(nci.graph)[name == x], mode = "in")

    # for each of the paths, iteratively find the longest and remove the subpaths
    ## leaving the longest distinct paths

    remaining.paths <- tmp.sp$res

    kept.paths <- list()

    while (length(remaining.paths) > 0) {
      longest.path <- which.max(lengths(remaining.paths))

      is.part.longest <- sapply(remaining.paths, function(y) {
        all(y %in% remaining.paths[[longest.path]])
      })

      kept.paths <- append(kept.paths, remaining.paths[longest.path])

      remaining.paths <- remaining.paths[is.part.longest == F]
    }

    nci.paths <- sapply(kept.paths, function(y) {
      y.dt <- data.table(code = y$name, index = seq_along(y$name))

      found.nci <- merge(nci[, .(code, use_name)], y.dt, by = "code", all = F)

      # don't link to the drug itself as that interferes with aggregating across synonyms
      found.nci <- found.nci[code != x, ]

      found.nci[order(index, decreasing = T), paste(paste(code, use_name, sep = " - "), collapse = " > ")]
    })

    data.table(code = x, nci = paste(nci.paths, collapse = "|"))
  }))

  # filter to those stemming from 'C1908 - Drug, Food, Chemical or Biomedical Material'

  drug.nci.terms <- drug.nci.terms[grepl("1908", nci) == T]

  found.drugs <- merge(found.drugs, drug.nci.terms, by = "code", all = F)

  cid.to.nci <- found.drugs[is.na(nci) == F, .(nci = paste(unique(nci), collapse = "|")), by = pubchem_cid]

  targetome$drugs <- merge(targetome$drugs, cid.to.nci, by = "pubchem_cid", all.x = T, all.y = F)

  targetome
}

#' Retrieve additional Uniprot information including annotation score
#'
#' Note this was written by ChatGPT 5.
#'
.uniprot.gene.and.score <- function(ids) {
  stopifnot(length(ids) >= 1)
  # Build a query like: (accession:P05067 OR accession:Q9Y2T1)
  q <- paste0("(", paste(sprintf("accession:%s", ids), collapse = " OR "), ")")

  url <- "https://rest.uniprot.org/uniprotkb/search"
  # Ask UniProt to only return the fields we need, as TSV
  params <- list(
    query  = q,
    fields = "accession,gene_primary,annotation_score",
    format = "tsv",
    size   = length(ids) # enough to return all rows in one go
  )

  resp <- httr::GET(url, query = params, httr::user_agent("R-UniProt-demo/1.0"))
  httr::stop_for_status(resp)

  # Parse the TSV to a data frame
  out <- readr::read_tsv(httr::content(resp, as = "raw"), show_col_types = FALSE)
  # Make columns easier to work with
  names(out) <- c("accession", "gene_name", "annotation_score")
  out
}

#' Add gene symbol information to targetome
#'
#' The targetome identifies targets using Uniprot identifiers.  Convert these
#' to gene symbols to simplify use in downstream applications.  As part of this
#' attempt to deal with the multi-mapping and ambiguity issues that typically
#' arise as part of translation from one type of identifier to another.
#' Internal function only to be used in \code{postprocess.targetome}.
#'
.add.symbols.targets <- function(tome, uniprot.file, entrez.file) {
  entrez <- fread(entrez.file)

  uniprot <- fread(uniprot.file, header = F)

  unip.gene.ids <- uniprot[V2 == "GeneID"]

  tome.gene.ids <- unip.gene.ids[V1 %in% unique(tome$uniprot_id)]

  # fix the multi-matching

  mm.genes <- tome.gene.ids[, .N, by = V3][N > 1]

  if (mm.genes[, .N] > 0) {
    mm.genes.prots <- tome.gene.ids[V3 %in% mm.genes$V3]

    mm.genes.annot <- .uniprot.gene.and.score(unique(mm.genes.prots$V1))

    mm.genes.prots <- merge(mm.genes.prots, mm.genes.annot, by.x = "V1", by.y = "accession")

    mm.genes.prots <- merge(mm.genes.prots, entrez[, .(V3 = as.character(GeneID), Symbol)], by = "V3")

    rm.gene <- mm.genes.prots[(gene_name == Symbol & annotation_score == 5) == F]

    tome.gene.ids <- merge(tome.gene.ids, rm.gene[, .(V1, V3, rm_gene = T)], by = c("V1", "V3"), all.x = T)
  } else {
    tome.gene.ids[, rm_gene := NA]
  }

  mm.prots <- tome.gene.ids[, .N, by = V1][N > 1]

  if (mm.prots[, .N] > 0) {
    mm.prots.genes <- tome.gene.ids[V1 %in% mm.prots$V1]

    mm.prots.annot <- .uniprot.gene.and.score(unique(mm.prots.genes$V1))

    mm.prots.genes <- merge(mm.prots.genes, mm.prots.annot, by.x = "V1", by.y = "accession", all.x = T)

    mm.prots.genes <- merge(mm.prots.genes, entrez[, .(V3 = as.character(GeneID), Symbol)], by = "V3", all.x = T)

    rm.prots <- mm.prots.genes[(Symbol == gene_name | grepl(";", gene_name)) == F]

    tome.gene.ids <- merge(tome.gene.ids, rm.prots[, .(V1, V3, rm_prots = T)], by = c("V1", "V3"), all.x = T)
  }

  tome.gene.ids <- tome.gene.ids[is.na(rm_gene) & is.na(rm_prots)]

  # add in symbols

  tome.gene.ids <- merge(tome.gene.ids, entrez[, .(V3 = as.character(GeneID), Symbol)], by = "V3")
  
  #stopifnot(tome.gene.ids[, .N, by = Symbol][N > 1, .N] == 0)

  tome.m <- merge(tome.gene.ids[, .(uniprot_id = V1, gene_id = V3, symbol = Symbol)], tome, by = "uniprot_id")

  tome.m
}


#' Post-process targetome
#'
#' Carries out the following post-processing steps:
#'   \itemize{
#'     \item Limit the targetome targets to only 'IC50', 'Kd', 'Ki' and 'EC50' types.
#'     \item Ensure the drugs and synonyms are filtered appropriately so that the list of drugs is consistent
#'     \item Fix irregularities in assay value.
#'     \item Add in gene symbols
#'     \item Add in drug information including inchi-keys where possible.
#'
#' @param tome.list
#'   A \code{list} containing three elements:
#'   \itemize{
#'     \item \code{drugs}: A \code{data.table} of drug entities.
#'     \item \code{targets}: A \code{data.table} of compound–target–reference-assay rows.
#'     \item \code{synonyms}: A \code{data.table} of drugs and synonyms.
#'   }
#'
#' @param uniprot A `data.table` produced by [form.uniprot.map()],
#' representing the mapping between UniProt accession IDs and other database
#' identifiers (e.g., Ensembl, GeneID). It must contain at least columns `V1`
#' (UniProt ID), `V2` (identifier type) and `V3` (mapped identifier, such as a protein accession).
#'
#' @param entrez.file An Entrez gene annotation file as downloaded from
#' \url{https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz}.
#'
#' @return
#' A modified version of the input \code{targetome} (same structure) specifically:
#' \itemize{
#'   \item \code{targetome$targets}
#'   \item \code{targetome$drugs}
#'   \item \code{targetome$synonyms}
#' }
#'
postprocess.targetome <- function(tome.list, uniprot.file, entrez.file) {
  # limit to desired assay types
  cur.targs <- tome.list$targets[assay_type %in% c("IC50", "Kd", "Ki", "EC50")]
  
  # add symbols to the targets
  cur.targs <- .add.symbols.targets(cur.targs, uniprot.file, entrez.file)

  limited.drugs <- tome.list$drugs[pubchem_cid %in% cur.targs$pubchem_cid]

  limited.targs <- cur.targs[pubchem_cid %in% limited.drugs$pubchem_cid]

  limited.syns <- tome.list$synonyms[pubchem_cid %in% limited.drugs$pubchem_cid]

  # Convert to nanomolar from micromolar for pubchem
  limited.targs[database == "pubchem_bioassay", assay_value := assay_value * 1000]

  # Adjust drug names
  limited.drugs[, pref_name2 := ifelse(is.na(pref_name), mesh, pref_name)]

  limited.drugs[, pref_name2 := ifelse(is.na(pref_name2), drug_name, pref_name2)]

  limited.drugs[, pref_name2 := safe_capwords(pref_name2)]

  limited.drugs[, `:=`(drug_name = NULL, pref_name = NULL)]

  limited.drugs[, drug_name := pref_name2]

  # Add additional drug info including inchi-key annotation to the targets
  limited.targs <- merge(limited.drugs[, .(pubchem_cid, inchi_key, drug_name, clinical_phase)],
    limited.targs,
    by = "pubchem_cid"
  )
  
  limited.targs[,evidence_level:="III"]

  list(drugs = limited.drugs, targets = limited.targs, synonyms = limited.syns)
}


#' Format and export the full drug set from targetome
#'
#' Given a comprehensive targetome (typically enriched with synonyms, NCI annotations, etc.),
#' prepare and save the drug annotation table and associated synonyms and target tables.
#'
#' @param tome.list
#'   A \code{list} containing at least three components:
#'   \itemize{
#'     \item \code{drugs}: A \code{data.table} of drug entities annotated with multiple metadata fields.
#'     \item \code{synonyms}: A \code{data.table} of drug synonyms (see \code{\link{add.synonym.table}}).
#'     \item \code{targets}: A \code{data.table} of compound–target–reference-assay rows.
#'   }
#'
#' @param tkg.preds A `data.table` of the Targetome Knowledge Graph predictions containing
#' columns for at least 'inchikey', 'uniprot_id', 'gene_id' and 'symbol'.
#'
#' @return
#' A character scalar giving the file path of the saved targets file,
#' and as a side-effect files containing drug and synonym information are generated
#' as well.
#'
write.full.drugs <- function(tome.list, tkg.preds, drugbank=NULL) {
  
  out.dir <- "outputs"

  limited.drugs <- tome.list$drugs[, .(
    pubchem_cid, inchi_key, clinical_phase, drug_name, chembl_id, therapeutic_flag,
    molecule_type, withdrawn_flag, chemical_probe, natural_product,
    atc, mesh, pharm_mesh, nci, inchi
  )]

  syns <- tome.list$synonyms

  # for the preds, get the requisite drug information from targetome expanded

  tkg.w.drug.info <- merge(
    limited.drugs[, .(pubchem_cid, inchi_key, drug_name, clinical_phase)],
    tkg.preds[, .(inchi_key = inchikey, uniprot_id, gene_id, symbol)],
    by = c("inchi_key")
  )

  stopifnot(tkg.w.drug.info[, .N] == tkg.preds[, .N])

  # add in extra information

  tkg.w.drug.info[, `:=`(
    pubmed_id = NA_character_, database = NA_character_,
    assay_type = "Predicted", assay_relation = "<=",
    assay_value = 100,
    evidence_level = "0"
  )]
  
  targs <- rbind(
    tome.list$targets,
    tkg.w.drug.info[, names(tome.list$targets), with = F]
  )
  
  if (is.null(drugbank) == F){
    
    db.targs <- merge(
      as.data.table(drugbank$targets), 
      drugbank$drugs[,.(pubchem_cid, inchi_key, clinical_phase)],
      by="pubchem_cid"
    )
    
    db.targs[,`:=`(database="drugbank", assay_type=NA_character_, assay_relation=NA_character_, 
                   assay_value=NA_real_, evidence_level="II")]
    
    targs <- rbind(
      targs,
      db.targs[,names(targs),with=F]
    )
    
    drugbank$drugs[,nci:=NA_character_]
    
    missing.drugs <- setdiff(drugbank$drugs$pubchem_cid, limited.drugs$pubchem_cid)
    
    limited.drugs <- rbind(
      limited.drugs,
      drugbank$drugs[pubchem_cid %in% missing.drugs,names(limited.drugs),with=F]
    )
    
    syns <- rbind(
      syns,
      drugbank$syns[pubchem_cid %in% missing.drugs, names(syns), with=F]
    )
    
    suff <- "_private_"
    
  }else{
    suff <- "_"
  }

  if (dir.exists(out.dir) == F) {
    dir.create(out.dir)
  }

  cur.date <- format(Sys.Date(), "%m-%d-%y")

  targets.file <- file.path(out.dir, paste0("targetome_expanded", suff, cur.date, ".tsv.gz"))

  fwrite(targs, file = targets.file, sep = "\t", col.names = T, row.names = F, quote = F)

  # Files output as a side-effect

  fwrite(limited.drugs, file = file.path(out.dir, paste0("targetome_expanded_drugs", suff, cur.date, ".tsv.gz")), sep = "\t", col.names = T, row.names = F, quote = T)

  fwrite(syns, file = file.path(out.dir, paste0("targetome_expanded_syns", suff, cur.date, ".tsv.gz")), sep = "\t", col.names = T, row.names = F, quote = T)

  targets.file
}
