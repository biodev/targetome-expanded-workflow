# Targetome-Extended Workflow

## Citation

Bottomly, D., Evans, N. and McWeeney, S.K. Expanding the set of high evidence drug-target interactions in the Cancer Targetome *Submitted*

## Download necessary files

``` bash

# PubChem (2026-01-30)

mkdir pubchem

cd pubchem

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Bioassay/Extras/bioactivities.tsv.gz

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-MeSH

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/MeSH-Pharm

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Parent.gz

wget https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Title.gz

## This is WHO ATC classification from limiting to drugs with ATC: 'PubChem: PubChem Compound TOC: WHO ATC Classification System':
PubChem_compound_list_L2iLuvbvk1OkeZFgExjYS01YQzjcjwYSfDcdXmcmD19nPzM.csv.gz

cd ..

# SMS

mkdir sorger

cd sorger

wget https://lsp.connect.hms.harvard.edu/smallmoleculesuite/_w_d85705b38b074ee2ae59bcf379a634a2/sms/assets/downloads/lsp_compound_dictionary.csv.gz

wget https://lsp.connect.hms.harvard.edu/smallmoleculesuite/_w_d85705b38b074ee2ae59bcf379a634a2/sms/assets/downloads/lsp_target_mapping.csv.gz

wget https://lsp.connect.hms.harvard.edu/smallmoleculesuite/_w_d85705b38b074ee2ae59bcf379a634a2/sms/assets/downloads/lsp_biochem.csv.gz

wget https://lsp.connect.hms.harvard.edu/smallmoleculesuite/_w_d85705b38b074ee2ae59bcf379a634a2/sms/assets/downloads/lsp_references.csv.gz

cd ..

# BindingDB

wget https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202601_tsv.zip

unzip BindingDB_All_202601_tsv.zip

# IUPHAR/BPS (2025.4)

mkdir iuphar

cd iuphar

wget https://www.guidetopharmacology.org/DATA/interactions.tsv

wget https://www.guidetopharmacology.org/DATA/ligands.tsv

cd ..

# Chembl

wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_36_sqlite.tar.gz

tar -xvzf chembl_36_sqlite.tar.gz

# NCI Thesaurus (2026-01-13)

wget https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Thesaurus.FLAT.zip

unzip Thesaurus.FLAT.zip

# UniProt (2026-01-28)

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

# Entrez (2026-01-30)

wget https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

# Need original targetome

wget https://raw.githubusercontent.com/ablucher/The-Cancer-Targetome/refs/heads/beta-V2/results_V2beta/Targetome_FullEvidence_210618_All.txt

# Need prediction data

Retrieve the DTI prediction results from https://github.com/biodev/targetome-expanded/ using git-lfs.

The file in question is: TargetomeKG_Complex2_DTI_preds_03-19-2026.csv.gz

**If desired**

# Drugbank

mkdir drugbank_files

cd drugbank_files

curl -Lfv -o filename.zip -u username:password https://go.drugbank.com/releases/VERSION/downloads/all-full-database

unzip filename.zip 

rm filename.zip

cd ..
```

## Run the main workflow

``` r
targets:tar_make()
```

## Requirements

### CRAN

-   targets
-   data.table
-   RSQLite
-   tidyverse
-   igraph
-   openxlsx

## Contact

Please open an Issue for questions or problems related to this workflow.

