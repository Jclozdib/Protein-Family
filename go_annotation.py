import requests
import json


input_file = "ground_truth_accessions.txt"
output_file_path = "ground_truth_annotations.tsv"

with open(input_file, "r") as f:
    accessions = set(line.strip() for line in f)


def get_request_url(accession_code):
    return f"https://rest.uniprot.org/uniprotkb/{accession_code}"


def get_go_annotations(code):
    url = get_request_url(code)
    response = requests.get(url, headers={"Accept": "application/json"})

    if response.status_code != 200:
        print(f"Error: {response.status_code} - {response.text}")
        return None

    try:
        data = response.json()
        cross_references = data.get("uniProtKBCrossReferences", [])
        go_annotations = [
            entry for entry in cross_references if entry.get("database") == "GO"
        ]
        return go_annotations

    except json.JSONDecodeError:
        print(f"Error decoding JSON for {code}")
        return None


with open(output_file_path, "w") as outfile:
    outfile.write("accession\tgo_term\tgo_name\n")
    for accession in accessions:
        annotations = get_go_annotations(accession)

        if annotations is None:
            print(f"Skipping {accession} due to errors.")
            continue

        for i in annotations:
            go_id = i["id"]
            go_name = i["properties"][0]["value"].split(":")[1]
            outfile.write(f"{accession}\t{go_id}\t{go_name}\n")

print(f"GO Annotations saved in {output_file_path}")
