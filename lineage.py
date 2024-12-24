import requests

input_file = "ground_truth_accessions.txt"

with open(input_file, "r") as f:
    accessions = set(line.strip() for line in f)

def get_request_url(accession_code):
    return f"https://rest.uniprot.org/uniprotkb/search?query=accession:{accession_code}&fields=lineage"


def get_lineage(url):
    response = requests.get(url, headers={"Accept": "application/json"})
    
    if response.status_code != 200:
        print(f"Error: {response.status_code} - {response.text}")
        return
        
    data = response.json()
    res = data.get("results")[0].get("lineages")
    return res

output = ''
print("Accession Code\t Lineage")

for a in accessions:
    lineage = get_lineage(get_request_url(a))
    parsed_lineage = " -> ".join(entry['scientificName'] for entry in lineage)
    output += f"{a}\t {parsed_lineage}\n"
print(output)

output_file = "lineage.txt"
with open(output_file, "w") as f:
    f.write(output)

print(f"Lineages saved to {output_file}")
