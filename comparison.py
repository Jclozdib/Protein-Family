import requests
import pandas as pd

api_url = "https://www.ebi.ac.uk/interpro/api/protein/reviewed/entry/pfam/PF03060/"

def fetch_api_data(url):
    results = []
    while url:
        response = requests.get(url, headers={"Accept": "application/json"})
        
        if response.status_code != 200:
            print(f"Error: {response.status_code} - {response.text}")
            break
        
        data = response.json()
        results.extend(data.get("results", []))
        url = data.get("next")
        
    return results

# Domain search results parsing
def parse_results(results):
    parsed_data = []
    for entry in results:
        metadata = entry["metadata"]
        protein_accession = metadata["accession"]
        protein_name = metadata["name"]
        organism = metadata["source_organism"]["scientificName"]
        length = metadata["length"]
        
        # Extract entries and locations
        for domain_entry in entry["entries"]:
            domain_accession = domain_entry["accession"]
            for location in domain_entry["entry_protein_locations"]:
                fragments = location["fragments"]
                for fragment in fragments:
                    parsed_data.append({
                        "Protein Accession": protein_accession,
                        "Protein Name": protein_name,
                        "Organism": organism,
                        "Protein Length": length,
                        "Domain Accession": domain_accession,
                        "Start": fragment["start"],
                        "End": fragment["end"],
                        "Score": location.get("score")
                    })
    return parsed_data

api_results = fetch_api_data(api_url)
parsed_data = parse_results(api_results)

# df
df_gt = pd.DataFrame(parsed_data)

# Psiblast parse
def parse_psiblast(file_path):
    columns = [
        "Seq_Ignore", "Accession", "Percentage", "Length", "Unknown1", 
        "Unknown2", "Unknown3", "Unknown4", "Start", "End", "E-value", "Unknown5"
    ]
    df = pd.read_csv(file_path, sep="\t", names=columns)

    # Extract Accession (between | in the 2nd column)
    df[['Ignore1', 'Accession', 'Ignore2']] = df['Accession'].str.split('|', expand=True)
    df = df[['Accession', 'Percentage', 'Length', 'Start', 'End', 'E-value']]  # Relevant columns
    return df

# HMM parse
def parse_hmm(file_path):
    columns = ["Accession", "Ali_Coord_From", "Ali_Coord_To", "E-value"]
    parsed_rows = []

    # Open the file and process line by line
    with open(file_path, "r") as file:
        for line in file:
            # Skip header lines and empty lines
            if line.startswith("#") or not line.strip():
                continue

            fields = line.split()

            # Relevant fields
            accession = fields[0].split("|")[1]  # Extract value between "|"
            ali_coord_from = int(fields[19])    # Alignment start coordinate
            ali_coord_to = int(fields[20])      # Alignment end coordinate
            e_value = float(fields[6])          # E-value

            parsed_rows.append([accession, ali_coord_from, ali_coord_to, e_value])

    df = pd.DataFrame(parsed_rows, columns=columns)
    return df

psiblast_file = "psiblast_full_output.txt"
hmm_file = "hmmer_results_domain.txt"           

# Parse PSI-BLAST and HMM results
df_psiblast = parse_psiblast(psiblast_file)
df_hmm = parse_hmm(hmm_file)

# Display sample dataframes
print("PSI-BLAST Results:")
print(df_psiblast.head())

print("\nHMM Results:")
print(df_hmm.head())

print("\nGround Truth:")
print(df_gt.head())
