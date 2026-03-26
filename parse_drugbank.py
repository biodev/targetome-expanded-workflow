import xml.etree.ElementTree as ET
import csv
import argparse
from pathlib import Path

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract DrugBank XML data and export to CSV."
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the DrugBank XML input file (e.g., full database.xml)"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the CSV output file (e.g., output_files/drugbank.csv)"
    )

    args = parser.parse_args()

    # Ensure output directory exists
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    # Parse the DrugBank XML input file into an ElementTree structure
    tree = ET.parse(args.input)
    root = tree.getroot()

    # Open the output CSV file for writing
    with open(args.output, 'w', newline='') as csvfile:

        # Define the CSV column headers
        fields = ['primary_drugbank_id', 'secondary_ids', 'pubchem_sid', 'pubchem_cid', 
                  'drug_name', 'status', 'type', 'action', 'target_name', 'uniprot_id', 
                  'uniprot_source', 'organism', 'pubmed_ids']
        
        # Create a CSV DictWriter object configured with the header fields
        drugwriter = csv.DictWriter(csvfile, fieldnames=fields)
        drugwriter.writeheader()

        # Iterate over every <drug> element in the XML
        for drugs in root.findall("{http://www.drugbank.ca}drug"):

            peptide_id=""
            source=""
            organism=""

            drug_type = drugs.attrib['type']
            name = drugs.find('{http://www.drugbank.ca}name').text

            # Collect primary and secondary DrugBank IDs
            drugbank_id = ""
            drugbank_ids = []
            for id in drugs.findall('{http://www.drugbank.ca}drugbank-id'):
                if ('primary' in id.attrib) and (id.attrib['primary'] == 'true'):
                    drugbank_id = id.text
                else:
                    drugbank_ids.append(id.text)
            
            # Collect drug approval/usage status from <groups>
            drug_status = []
            groups = drugs.find('{http://www.drugbank.ca}groups')
            for group in groups.findall('{http://www.drugbank.ca}group'):
                drug_status.append(group.text)
            
            # Collect PubChem external identifiers
            eids = drugs.find('{http://www.drugbank.ca}external-identifiers')
            pubchem_SID=[]
            pubchem_CID=[]
            for eid in eids.findall('{http://www.drugbank.ca}external-identifier'):
                if eid.find('{http://www.drugbank.ca}resource').text == "PubChem Substance":
                    pubchem_SID.append(eid.find('{http://www.drugbank.ca}identifier').text)
                if eid.find('{http://www.drugbank.ca}resource').text == "PubChem Compound":
                    pubchem_CID.append(eid.find('{http://www.drugbank.ca}identifier').text)

            # Process drug targets (one row written per target)
            tar = drugs.find("{http://www.drugbank.ca}targets")
            if (len(tar) > 0):
                for target in tar.findall('{http://www.drugbank.ca}target'):
                    
                    pids=[]
                    
                    # Collect actions describing how the drug interacts with the target
                    action_list = []
                    actions = target.find('{http://www.drugbank.ca}actions')
                    for action in actions.findall('{http://www.drugbank.ca}action'):
                        action_list.append(action.text)
                    
                    # Extract target polypeptide info (UniProt ID + source database)
                    polypeptide=target.find('{http://www.drugbank.ca}polypeptide')
                    if (polypeptide is not None):
                        peptide_id=polypeptide.attrib.get('id')
                        source=polypeptide.attrib.get('source')
                    
                    # Extract organism name for the target
                    organism=target.find('{http://www.drugbank.ca}organism').text

                    # Target name (if present)
                    t_name = target.find('{http://www.drugbank.ca}name')
                    if (t_name is not None):

                        # Collect PubMed IDs from target references (if present)
                        refs = target.find('{http://www.drugbank.ca}references')
                        if (refs is not None):
                            for articles in refs.findall('{http://www.drugbank.ca}articles'):
                                for article in articles.findall('{http://www.drugbank.ca}article'):
                                    pid=article.find('{http://www.drugbank.ca}pubmed-id').text
                                    if (pid is not None):
                                        pids.append(pid)

                    # Write one CSV row for this drug-target relationship
                    drugwriter.writerow({'primary_drugbank_id':drugbank_id, 
                            'secondary_ids':'|'.join(drugbank_ids),
                            'pubchem_sid':'|'.join(pubchem_SID),
                            'pubchem_cid': '|'.join(pubchem_CID),
                            'drug_name':name,
                            'status':'|'.join(drug_status), 
                            'type':drug_type, 
                            'action':'|'.join(action_list),
                            'target_name':t_name.text,
                            'uniprot_id':peptide_id,
                            'uniprot_source':source,
                            'organism':organism,
                            'pubmed_ids':'|'.join(pids)})
    
    output_path = Path(args.output)
    drug_info_filename = output_path.parent / (output_path.stem + "_drug_info" + output_path.suffix)

    with open(drug_info_filename, 'w', newline='') as csvfile:
        fields = ['primary_drugbank_id', 'secondary_ids', 'pubchem_sid', 'pubchem_cid', 'drug_name', 'status', 'type']
        drugwriter = csv.DictWriter(csvfile, fieldnames=fields)
        drugwriter.writeheader()
        for drugs in root.findall("{http://www.drugbank.ca}drug"):
                drug_type = drugs.attrib['type']
                tname=""
                peptide_id=""
                source=""
                organism=""
                name = drugs.find('{http://www.drugbank.ca}name').text
                drugbank_id = ""
                drugbank_ids = []
                for id in drugs.findall('{http://www.drugbank.ca}drugbank-id'):
                        if ('primary' in id.attrib) and (id.attrib['primary'] == 'true'):
                                drugbank_id = id.text
                        else:
                                drugbank_ids.append(id.text)
                drug_status = []
                groups = drugs.find('{http://www.drugbank.ca}groups')
                for group in groups.findall('{http://www.drugbank.ca}group'):
                        drug_status.append(group.text)
                eids = drugs.find('{http://www.drugbank.ca}external-identifiers')
                pubchem_SID=[]
                pubchem_CID=[]
                for eid in eids.findall('{http://www.drugbank.ca}external-identifier'):
                        if eid.find('{http://www.drugbank.ca}resource').text == "PubChem Substance":
                                pubchem_SID.append(eid.find('{http://www.drugbank.ca}identifier').text)
                        if eid.find('{http://www.drugbank.ca}resource').text == "PubChem Compound":
                                pubchem_CID.append(eid.find('{http://www.drugbank.ca}identifier').text)
                drugwriter.writerow({'primary_drugbank_id':drugbank_id, 
						  'secondary_ids':'|'.join(drugbank_ids),
						  'pubchem_sid':'|'.join(pubchem_SID),
						  'pubchem_cid': '|'.join(pubchem_CID),
						  'drug_name':name,
						  'status':'|'.join(drug_status), 
						  'type':drug_type})

        