library(targets)

source("R/utils.R")
source("R/pubchem.R")
source("R/sorger.R")
source("R/bindingdb.R")
source("R/iuphar.R")
source("R/targetome.R")
source("R/eval.R")
source("R/validate_pubmed_ids.R")
source("R/drugbank.R")

list(
  tar_target(
    uniprot_map,
    form.uniprot.map("HUMAN_9606_idmapping.dat.gz"),
    packages = c("data.table")
  ),

  # Begin processing

  tar_target(
    bioassay,
    parse.bioassay("pubchem/bioactivities.tsv.gz", uniprot_map),
    packages = c("data.table")
  ),
  tar_target(
    bioassay_w_values,
    remove.bioassay.nas(bioassay),
    packages = c("data.table")
  ),
  tar_target(
    sorger_pubchem_map,
    make.sorger.pubchem.map(
      "sorger/lsp_compound_dictionary.csv.gz",
      "sorger/lsp_biochem.csv.gz",
      "pubchem/CID-InChI-Key.gz",
      "pubchem/CID-Title.gz",
      "pubchem/CID-Parent.gz"
    ),
    packages = "data.table"
  ),
  tar_target(
    sorger_db,
    form.sorger.db(
      sorger_pubchem_map,
      "sorger/lsp_biochem.csv.gz",
      "sorger/lsp_target_mapping.csv.gz",
      "sorger/lsp_references.csv.gz",
      uniprot_map
    ),
    packages = c("data.table")
  ),
  tar_target(
    bioassay_sorger,
    add.sorger.to.targetome(bioassay_w_values, sorger_db),
    packages = c("data.table")
  ),
  tar_target(
    bdb,
    form.bindingdb.db(uniprot_map, "BindingDB_All.tsv"),
    packages = c("tidyverse", "data.table")
  ),
  tar_target(
    bioassay_sorger_bdb,
    add.bindingdb.to.targetome(bioassay_sorger, bdb),
    packages = c("data.table")
  ),
  tar_target(
    iuphar,
    form.iuphar.db(
      "iuphar/interactions.tsv",
      "iuphar/ligands.tsv",
      uniprot_map
    ),
    packages = "data.table"
  ),
  tar_target(
    bioassay_sorger_bdb_iuphar,
    add.iuphar.to.targetome(bioassay_sorger_bdb, iuphar),
    packages = "data.table"
  ),
  tar_target(
    pharm_mesh,
    combine.mesh.pubchem("pubchem/CID-MeSH", "pubchem/MeSH-Pharm"),
    packages = "tidyverse"
  ),
  tar_target(
    initial_targetome,
    form.initial.targetome(
      bioassay_sorger_bdb_iuphar, pharm_mesh,
      "pubchem/CID-Synonym-filtered.gz", "pubchem/CID-Parent.gz",
      "pubchem/CID-Title.gz",
      "pubchem/PubChem_compound_list_L2iLuvbvk1OkeZFgExjYS01YQzjcjwYSfDcdXmcmD19nPzM.csv.gz"
    ),
    packages = "data.table"
  ),
  tar_target(
    targetome_w_inchi,
    add.inchi(initial_targetome, "pubchem/CID-InChI-Key.gz"),
    packages = "data.table"
  ),
  tar_target(
    targetome_w_chembl,
    add.chembl.annots(targetome_w_inchi, "chembl_36/chembl_36_sqlite/chembl_36.db"),
    packages = c("data.table", "RSQLite")
  ),
  tar_target(
    targetome_w_syns,
    add.synonym.table(targetome_w_chembl),
    packages = c("data.table", "tidyverse")
  ),
  tar_target(
    targetome_annotated,
    add.nci.annots(targetome_w_syns, "Thesaurus.txt"),
    packages = c("data.table", "igraph")
  ),
  tar_target(
    targetome_pp,
    postprocess.targetome(
      targetome_annotated,
      "HUMAN_9606_idmapping.dat.gz",
      "Homo_sapiens.gene_info.gz"
    ),
    packages = c("data.table")
  ),

  # Read in predictions

  tar_target(
    tkg_preds,
    .add.symbols.targets(
      fread("TargetomeKG_Complex2_DTI_preds_03-19-2026.csv.gz"),
      "HUMAN_9606_idmapping.dat.gz",
      "Homo_sapiens.gene_info.gz"
    ),
    packages = "data.table"
  ),
  tar_target(
    saved_targetome_pp,
    write.full.drugs(targetome_pp, tkg_preds),
    packages = c("data.table"),
    format = "file"
  ),

  # Evaluation compared to (original) targetome

  tar_target(
    orig_targetome,
    preprocess.targetome("Targetome_FullEvidence_210618_All.txt", targetome_pp, tibble::as_tibble(uniprot_map)),
    packages = "tidyverse"
  ),
  tar_target(
    orig_targetome_v_pmid,
    subset.to.valid.pmids(orig_targetome),
    packages = "tidyverse"
  ),
  tar_target(
    orig_exp_eval_xl,
    compare.orig.w.expanded(tibble::as_tibble(targetome_pp$targets), orig_targetome_v_pmid),
    packages = c("tidyverse", "openxlsx"),
    format = "file"
  ),

  # Evaluation compared to originating databases

  tar_target(
    exp_db_eval_xl,
    compare.expanded.w.dbs(tibble::as_tibble(targetome_pp$targets),
      "pubchem/CID-Parent.gz",
      pubchem = tibble::as_tibble(bioassay_w_values),
      sorger = tibble::as_tibble(sorger_db),
      bindingdb = tibble::as_tibble(bdb),
      iuphar = tibble::as_tibble(iuphar)
    ),
    packages = "tidyverse",
    format = "file"
  )
  
  # Remove comments (and add comment above) if drugbank is available

  # tar_target(
  #   drugbank_inters,
  #   get.drugbank.interactions("drugbank_files/full\\ database.xml", "pubchem/CID-Parent.gz", tibble::as_tibble(uniprot_map)),
  #   packages = "tidyverse"
  # ),
  # tar_target(
  #   drugbank_st,
  #   .add.symbols.targets(
  #     drugbank_inters,
  #     "HUMAN_9606_idmapping.dat.gz",
  #     "Homo_sapiens.gene_info.gz"
  #   ),
  #   packages = "data.table"
  # ),
  # tar_target(
  #   drugbank_inters_annot,
  #   add.inchi.chembl.drugbank(drugbank_st, "pubchem/CID-InChI-Key.gz", 
  #                             "chembl_36/chembl_36_sqlite/chembl_36.db", pharm_mesh,
  #                             "pubchem/PubChem_compound_list_L2iLuvbvk1OkeZFgExjYS01YQzjcjwYSfDcdXmcmD19nPzM.csv.gz",
  #                             "pubchem/CID-Synonym-filtered.gz"),
  #   packages = c("data.table", "RSQLite")
  # ),
  # tar_target(
  #   saved_targetome_pp_w_db,
  #   write.full.drugs(targetome_pp, tkg_preds, drugbank_inters_annot),
  #   packages = c("data.table"),
  #   format = "file"
  # )
)
