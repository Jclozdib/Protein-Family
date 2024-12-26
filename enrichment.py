import gzip
import json
import os
import shutil
import requests
from scipy.stats import fisher_exact
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

pd.set_option("display.max_columns", None)

# # Download SwissProt XML
compressed_file = "uniprot_sprot.xml.gz"
xml_file = "uniprot_sprot.xml"
url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz"

if not os.path.exists(xml_file):
    print(f"{xml_file} not found. Downloading and extracting...")
    print(f"This file is around 8GB. Get a coffee ☕️, it might take a while!")

    print(f"Downloading {compressed_file}...")
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(compressed_file, "wb") as f:
            shutil.copyfileobj(response.raw, f)
        print(f"Downloaded {compressed_file}.")
    else:
        print(
            f"Failed to download {compressed_file}. Status code: {response.status_code}"
        )
        exit(1)

    print(f"Extracting {compressed_file}...")
    with gzip.open(compressed_file, "rb") as f_in:
        with open(xml_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Extracted to {xml_file}.")

    os.remove(compressed_file)
    print(f"Removed {compressed_file}.")
else:
    print(f"{xml_file} already exists. No action needed.")


def extract_go_terms_from_xml_iterative(xml_file, go_terms_of_interest):
    go_annotations = {}
    context = ET.iterparse(xml_file, events=("start", "end"))
    context = iter(context)

    for event, elem in context:
        if event == "end" and elem.tag.endswith("entry"):
            accession = elem.find(".//{http://uniprot.org/uniprot}accession")
            if accession is not None:
                accession = accession.text
                # Extract GO terms
                for db_ref in elem.findall(
                    ".//{http://uniprot.org/uniprot}dbReference"
                ):
                    if db_ref.attrib.get("type") == "GO":
                        go_id = db_ref.attrib["id"]
                        go_name = db_ref.find(".//{http://uniprot.org/uniprot}property")
                        if go_id in go_terms_of_interest:
                            if go_id not in go_annotations:
                                go_annotations[go_id] = {
                                    "id": go_id,
                                    "name": (
                                        go_name.attrib["value"]
                                        if go_name is not None
                                        else "Unknown"
                                    ),
                                    "proteinCount": 0,
                                }
                            go_annotations[go_id]["proteinCount"] += 1

            # Clear the element to free memory
            elem.clear()

    return list(go_annotations.values())


swiss_prot_go_annotations_filename = "swiss_prot_go_annotations.json"

ground_truth = pd.read_csv("ground_truth_annotations.tsv", sep="\t")

# Load the ground truth annotations
go_terms_of_interest = ground_truth["go_term"].unique()

if not os.path.exists(swiss_prot_go_annotations_filename):
    print("Parsing XML incrementally...")
    go_annotations = extract_go_terms_from_xml_iterative(xml_file, go_terms_of_interest)
    # save go_annotations to a file
    with open(swiss_prot_go_annotations_filename, "w") as f:
        json.dump(go_annotations, f)
    print("Finished parsing XML.")
else:
    print("Loading pre-parsed GO annotations...")
    with open(swiss_prot_go_annotations_filename, "r") as f:
        go_annotations = json.load(f)
    print("GO annotations loaded.")

# Get unique proteins and total counts
total_proteins = sum(annotation["proteinCount"] for annotation in go_annotations)

# Prepare a lookup table from the refactored go_annotations
go_annotation_lookup = {
    annotation["id"]: annotation["proteinCount"] for annotation in go_annotations
}

# Analyze statistics for each go_term
results = []
for go_term in ground_truth["go_term"].unique():
    total_with_term = go_annotation_lookup.get(go_term, 0)

    # Access selected proteins for the current GO term
    selected_proteins = ground_truth.loc[
        ground_truth["go_term"] == go_term, "accession"
    ].tolist()

    set_with_term = len(selected_proteins)
    rest_with_term = total_with_term - set_with_term

    set_not_with_term = len(selected_proteins) - set_with_term
    rest_not_with_term = total_proteins - total_with_term - set_not_with_term

    go_name = ground_truth.loc[ground_truth["go_term"] == go_term, "go_name"].iloc[0]

    df = pd.DataFrame(
        {
            "name": [go_name],
            "id": [go_term],
            "total": [total_with_term],
            "set": [set_with_term],
            "rest": [rest_with_term],
            "not_set": [set_not_with_term],
            "not_rest": [rest_not_with_term],
        }
    )

    # Calculate ratios and fold enrichment
    df["set_ratio"] = df["set"] / (df["not_set"] + 1e-9)
    df["rest_ratio"] = df["rest"] / (df["not_rest"] + 1e-9)
    df["fold"] = df["set_ratio"] / (df["rest_ratio"] + 1e-9)
    # removing go terms that don't show up in swissprot
    # Cannot do fisher test with negative values, swissprot might not contain these annotations, hence rest will be negative
    df = df[(df["total"] > 0) & (df["rest"] > 0)]
    results.append(df)

final_results = pd.concat(results, ignore_index=True)

final_results["oddsratio"], final_results["p_value"] = zip(
    *final_results.apply(
        lambda row: fisher_exact(
            [[row["set"], row["not_set"]], [row["rest"], row["not_rest"]]],
            alternative="greater",
        ),
        axis=1,
    )
)

final_results["significant"] = final_results["p_value"] < 0.005
final_results = final_results.sort_values("p_value", ascending=True)
print(final_results.head(30))
final_results.to_csv("enrichment_results.tsv", sep="\t", index=False)
