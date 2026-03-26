#' Preprocess Cancer Targetome Data
#'
#' This function preprocesses the Cancer Targetome dataset by standardizing drug names,
#' filtering for high-quality evidence, harmonizing protein targets to UniProt IDs,
#' and mapping drugs to PubChem compound IDs (CIDs).
#'
#' @param targetome.file Character string. Path to the Targetome flat file containing
#'   drug-target interaction data. Expected format is tab-delimited with columns including
#'   Drug_Found, EvidenceLevel_Assigned, Target_Species, Assay_Type, Target_UniProt,
#'   Reference, Assay_Relation, etc.
#'   Reference file: https://raw.githubusercontent.com/ablucher/The-Cancer-Targetome/refs/heads/beta-V2/results_V2beta/Targetome_FullEvidence_210618_All.txt
#' @param cur.syns A `tibble` containing at least columns for pubchem_cid, inchi_key, synonyms and lower_name which
#'   is a lower-case version of the 'synonyms' column
#' @param unip.map A `tibble` as derived from the `form.uniprot.map` function in `utils.R`
#'
#' @return A tibble/data frame containing the preprocessed targetome data with renamed columns:
#'   Reference -> pubmed_id, Assay_Type -> assay_type, Assay_Relation -> assay_relation
#'
#' @details
#' The function performs several data cleaning steps:
#' \itemize{
#'   \item Removes drug class prefixes (e.g., "Aurora Kinase Inhibitor MLN8054" -> "MLN8054")
#'   \item Corrects known misspellings (AT-101, Baiclein, Ruxolotinib, Vargetef)
#'   \item Strips chemical formulation suffixes (Hydrochloride, Sulfate, Citrate, etc.)
#'   \item Harmonized UniProt IDs mapped to canonical identifiers
#'   \item Maps drugs to Targetome-expanded version using a synonym frequency approach to resolve ambiguities
#' }
#'

preprocess.targetome <- function(targetome.file, tome.list, unip.map) {
  # Initial preprocessing

  tome.flat <- readr::read_delim(targetome.file, quote = "")

  # If a drug is specified as 'xxx inhibitor yyy' then remove the 'xxx inhibitor ' portion.
  # See below
  # filter(tome.flat, grepl("^.+Inhibitor ", Drug_Found)) %>%
  #    {unique(.$Drug_Found)}

  # e.g. "Aurora Kinase Inhibitor MLN8054"          "Akt Inhibitor MK2206"                     "CDK Inhibitor AT7519"

  tome.flat <- mutate(tome.flat, pp_drug = sub("^.+Inhibitor ", "", Drug_Found))

  # Fix a few mispellings, near synonym misses

  tome.flat <- mutate(tome.flat, pp_drug = ifelse(pp_drug == "AT-101", "AT101", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = ifelse(pp_drug == "Baiclein", "Baicalein", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = ifelse(pp_drug == "Ruxolotinib", "Ruxolitinib", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = ifelse(pp_drug == "Vargetef", "Vargatef", pp_drug))

  # Adjust for a few cases where there is a specific formulation specified (e.g. 'Procarbazine Hydrochloride' -> 'Procarbazine')

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Hydrochloride", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Sulfate", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Triacetate", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Citrate", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Acetate", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Phosphate", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Disodium", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = sub(" Olamine", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = ifelse(pp_drug == "Liposomal Doxorubicin ", "Doxorubicin", pp_drug))

  # Finally remove another prefix

  tome.flat <- mutate(tome.flat, pp_drug = sub("BH3 Mimetic ", "", pp_drug))

  tome.flat <- mutate(tome.flat, pp_drug = tolower(pp_drug))

  tome.lev3 <- filter(tome.flat, EvidenceLevel_Assigned == "III" & Target_Species == "Homo sapiens" & Assay_Type %in% c("IC50", "Kd", "Ki", "EC50"))

  # Harmonize drug names
  
  ## get rid of trailing parentheses in say 'Ixazomib Citrate  (MLN-9708)'
  
  cur.syns <- mutate(tibble::as_tibble(tome.list$synonyms),
                     lower_name = sub("\\s+[\\(\\[].+", "", lower_name)
                     # lower_name=sub(" ", "", lower_name)
  )
  cur.syns <- unique(cur.syns)
  
  ## for each synonym, count the frequency of each cid as a way of choosing the best one
  
  # summarize(cur.syns, n=n(), .by=lower_name) %>%
  #    filter(n > 1)
  
  syn.freq <- summarize(cur.syns, n = n(), .by = c("lower_name", "pubchem_cid")) %>%
    arrange(desc(n)) %>%
    filter(!duplicated(lower_name))
  
  tome.lev3.drugs <- tome.lev3 %>%
    select(pp_drug) %>%
    unique()
  
  # left join to be exported
  
  tome.lev3.drug.status <- left_join(
    tome.lev3.drugs,
    syn.freq,
    by=c("pp_drug"="lower_name")
  )
  
  # add in other drug information
  
  tome.lev3.drug.status <- left_join(
    tome.lev3.drug.status,
    select(tibble::as_tibble(tome.list$drugs), -nci, -atc),
    by="pubchem_cid"
  ) %>%
    rename(orig_drug=pp_drug, pubchem_cid_from_new=pubchem_cid)
  
  openxlsx::write.xlsx(list(tome_lev3_drugs=tome.lev3.drug.status), 
                       file="outputs/orig_targetome_drugs.xlsx")
  
  # inner join to be kept for downstream analyses
  
  tome.drug.map <- inner_join(
    tome.lev3.drugs,
    syn.freq,
    by = c("pp_drug" = "lower_name")
  )
  
  stopifnot(nrow(tome.lev3.drugs) - 6 == nrow(tome.drug.map))
  
  tome.lev3 <-
    inner_join(
      tome.lev3,
      tome.drug.map,
      by = "pp_drug"
    )
  
  # Harmonize uniprot

  tome.lev3.targs <- tome.lev3 %>%
    select(Target_UniProt) %>%
    unique()

  tome.lev3.targs <- inner_join(
    tome.lev3.targs,
    unip.map,
    by = c(Target_UniProt = "V3")
  )

  tome.lev3.targs %>%
    summarize(n = n(), .by = Target_UniProt) %>%
    filter(n > 1) %>%
    {
      stopifnot(nrow(.) == 0)
    }

  tome.lev3.targs %>%
    summarize(n = n(), .by = V1) %>%
    filter(n > 1) %>%
    {
      stopifnot(nrow(.) == 0)
    }

  tome.lev3.mapped <-
    inner_join(
      tome.lev3,
      select(tome.lev3.targs, Target_UniProt, uniprot_id = V1),
      by = c("Target_UniProt")
    )

  tome.lev3.mapped <- rename(tome.lev3.mapped,
    pubmed_id = Reference,
    assay_type = Assay_Type,
    assay_relation = Assay_Relation
  )

  tome.lev3.mapped
}

#' Subset the original Targetome to valid PubMed IDs
#'
#' The original Targetome contained reference to other entities besides PubMed IDs
#' such as information from Chembl and PubChem assays.  Here we query PubChem to
#' filter out any invalid entries.
#'
#' @param tome.lev3.mapped
#'  A `tibble` containing at least a `pubmed_id` column containing potential PubMed
#'  IDs to be queried.
#'
#'
#' @return
#'  A filtered `tibble` with only valid PubMed IDs
#'
subset.to.valid.pmids <- function(tome.lev3.mapped) {
  unique.refs <- unique(tome.lev3.mapped$pubmed_id)
  unique.refs <- unique.refs[is.na(unique.refs) == F]

  refs.val <- validate_pubmed_ids(unique.refs)

  valid.pmids <- refs.val[refs.val$is_valid == T, ]

  tome.lev3.mapped.val <- filter(tome.lev3.mapped, pubmed_id %in% valid.pmids$pmid)

  tome.lev3.mapped.val
}

#' Compare Original and Expanded Targetome Datasets
#'
#' This function compares two versions of targetome drug-target interaction data
#' (the original and  expanded dataset) to assess agreement, generating an Excel
#' workbook with agreement and details on missing interactions.
#'
#' @param cur.targs The expanded/current targetome dataset containing
#'   drug-target interactions with columns: pubchem_cid, uniprot_id, pubmed_id,
#'   assay_type, assay_relation, assay_value, and other metadata.
#' @param orig.targs The original targetome dataset with columns:
#'   pubchem_cid, uniprot_id, pubmed_id, assay_type, assay_relation, Assay_Value
#'   and other metadata. Multiple measurements for the same interaction will be
#'   summarized using the median assay value.
#'
#' @return The function does not explicitly return a value but writes an Excel file
#'   to "outputs/orig_exp_eval_summary.xlsx" containing three sheets:
#'   \itemize{
#'     \item \strong{agreement_orig}: Summary statistics showing agreement between datasets
#'           for all interactions and interactions <100nM
#'     \item \strong{missing_expanded}: Interactions present in original but missing
#'           in expanded dataset (expected to be empty based on assertions)
#'
compare.orig.w.expanded <- function(cur.targs, orig.targs) {
  cur.targs <- mutate(cur.targs,
    assay_value_key = as.character(round(assay_value, digits = 7)),
    pubmed_id = as.character(pubmed_id),
    tome_ext = T
  )


  summarized.tome <- summarize(orig.targs,
    assay_value = median(as.numeric(Assay_Value)),
    .by = c("pubchem_cid", "uniprot_id", "pubmed_id", "assay_type", "assay_relation")
  )

  summarized.tome <- mutate(summarized.tome,
    assay_value_key = as.character(round(assay_value, digits = 7)),
    tome_orig = T
  )

  # Are all interactions accounted for?

  unique.expanded <- unique(select(cur.targs, pubchem_cid, uniprot_id, tome_ext))
  unique.tome <- unique(select(summarized.tome, pubchem_cid, uniprot_id, tome_orig))

  all.inters <-
    right_join(
      unique.expanded,
      unique.tome,
      by = c("pubchem_cid", "uniprot_id")
    )

  all.both <- filter(all.inters, tome_orig == T & tome_ext == T) %>%
    nrow()

  all.total <- all.inters %>% nrow()

  # Are interactions <= 100nM accounted for?

  lt100.expanded <- unique(select(filter(cur.targs, assay_relation %in% c("=", "<", "<=", "<<", "~") & assay_value < 100), pubchem_cid, uniprot_id, tome_ext))
  lt100.tome <- unique(select(filter(summarized.tome, assay_relation %in% c("=", "<", "<=", "<<", "~") & assay_value < 100), pubchem_cid, uniprot_id, tome_orig))

  lt100.inters <- right_join(
    lt100.expanded,
    lt100.tome,
    by = c("pubchem_cid", "uniprot_id")
  )

  lt100.both <- filter(lt100.inters, tome_orig == T & tome_ext == T) %>%
    nrow()

  lt100.total <- lt100.inters %>% nrow()

  # Record the missing ones from both categorizations (lt100 should be a subset)

  missing.in.ext <- inner_join(
    summarized.tome,
    filter(all.inters, is.na(tome_ext)),
    by = c("pubchem_cid", "uniprot_id", "tome_orig")
  ) %>%
    filter(is.na(tome_ext))

  select(missing.in.ext, pubchem_cid, uniprot_id) %>%
    unique() %>%
    nrow() %>%
    {
      stopifnot(nrow(.) == 0)
    }

  new.inters <- left_join(
    unique.expanded,
    unique.tome,
    by = c("pubchem_cid", "uniprot_id")
  ) %>%
    filter(is.na(tome_orig))

  stopifnot((nrow(new.inters) + all.both) == nrow(unique.expanded))

  new.exp <- new.inters %>%
    nrow()

  all.exp.total <- unique.expanded %>% nrow()

  # Make a summary table of interactions

  inter.summary <- tribble(
    ~comparison, ~count, ~total,
    "All_orig", all.both, all.total,
    "lt100_orig", lt100.both, lt100.total,
    "New_exp", new.exp, all.exp.total
  ) %>%
    mutate(
      percentage = scales::percent(count / total, accuracy = 0.001)
    )
  
  # Compute summary of coverage
  
  cov.expanded <- summarize(unique.expanded, n_expanded=n(), .by=uniprot_id)
  
  cov.orig <- summarize(unique.tome, n_orig=n(), .by=uniprot_id)

  cov.summary <- full_join(
    cov.expanded,
    cov.orig,
    by=c("uniprot_id")
  )
  
  openxlsx::write.xlsx(
    list(
      agreement_orig = inter.summary,
      missing_expanded = missing.in.ext,
      target_coverage = cov.summary
    ),
    file = "outputs/orig_exp_eval_summary.xlsx"
  )

  "outputs/orig_exp_eval_summary.xlsx"
}


#' Compare Targetome-expanded with Other databases
#'
#' Evaluates the overlap between the Targetome-expanded DTIs and
#' other supplied DTI databases. Generates a
#' summary report showing the frequency of database overlap
#'
#' @param targetome Current Targetome-expanded interaction data
#'   that must contain at least the columns 'pubchem_cid', 'uniprot_id', 'pubmed_id', 
#'   'assay_type', 'assay_relation' and 'assay_value'.
#' @param cid.par.file PubChem parent file.
#' @param ... Databases to compare, must contain at least 'pubchem_cid' and 'uniprot_id'
#' columns.
#'
#' @return The file path to the generated Excel summary file
#'   ("outputs/targetome_db_eval.xlsx").  As a side effect a
#'   file with details is generated (outputs/targetome_db_representation.csv).
#'
compare.expanded.w.dbs <- function(targetome, cid.par.file, ...) {
  
  unique.tome <- mutate(targetome, lt100=assay_value < 100) %>%
   select(
     pubchem_cid, uniprot_id, 
     pubmed_id, assay_type,
     assay_relation, lt100) %>%
    unique()

  data.list <- list(...)

  data.list$pubchem <- rename(
    data.list$pubchem,
    pubchem_cid = "CID",
    pubmed_id = "PMID",
    assay_type = "Activity Name",
    assay_relation = "Activity Qualifier"
  ) %>%
    mutate(assay_value = `Activity Value` * 1000)

  cid.par <- readr::read_delim(cid.par.file)
  
  for (data.source in names(data.list)) {
    print(data.source)
    col.name <- "in_ds"
    names(col.name) <- data.source

    unique.dt <- mutate(data.list[[data.source]], 
                        lt100 = assay_value < 100) %>%
      select(
        pubchem_cid, uniprot_id, 
        pubmed_id, assay_type,
        assay_relation, lt100) %>%
      unique() %>%
      mutate(in_ds = T)
    
    if (class(unique.dt$pubmed_id) != "character"){
      unique.dt$pubmed_id <- as.character(unique.dt$pubmed_id)
    }

    unique.w.pars <- inner_join(
      select(cid.par, pubchem_cid = `1`, parent_cid = `2`),
      unique.dt,
      by = "pubchem_cid"
    ) %>%
      select(-pubchem_cid) %>%
      rename(pubchem_cid = "parent_cid") %>%
      unique()

    unique.tome <- left_join(
      unique.tome,
      unique.w.pars,
      by = c("pubchem_cid", "uniprot_id", "pubmed_id", "assay_type", "assay_relation", "lt100")
    )

    unique.tome <- mutate(
      unique.tome,
      in_ds = ifelse(is.na(in_ds), F, in_ds)
    ) %>%
      rename(all_of(col.name))
  }

  tome.repr.summary <- pivot_longer(unique.tome, cols = -c(1:6), names_to = "database") %>%
    summarize(
      counted = sum(value),
      total = n(),
      dbs = paste(database[value == T], collapse = ";"),
      .by = c("pubchem_cid", "uniprot_id", "pubmed_id", "assay_type", "assay_relation", "lt100")
    )

  readr::write_csv(tome.repr.summary, file = "outputs/targetome_db_representation.csv")

  by.db.counts <- summarize(tome.repr.summary, n = n(), .by = counted)
  by.dbs <- summarize(tome.repr.summary, n = n(), .by = dbs)

  openxlsx::write.xlsx(list(db_counts = by.db.counts, by_dbs = by.dbs), file = "outputs/targetome_db_eval.xlsx")

  "outputs/targetome_db_eval.xlsx"
}