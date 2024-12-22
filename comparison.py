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

# Results parsing
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

print(df_gt.head())  # Display the first few rows
