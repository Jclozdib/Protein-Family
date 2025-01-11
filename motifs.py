import re
import csv
from Bio import SeqIO
import json
from collections import defaultdict


# step 1: Parse motifs from ProSite and ELM files
def parse_prosite_motifs(prosite_file):

    motifs = {}
    with open(prosite_file, "r") as file:
        current_id = None
        for line in file:
            if line.startswith("ID"):
                current_id = line.split()[1].strip()
            elif line.startswith("PA"):
                pattern = line.split()[1].strip()
                motifs[current_id] = pattern_to_regex(pattern)
    return motifs


def pattern_to_regex(pattern):
    regex = pattern.replace("{", "[^").replace("}", "]").replace("-", "")
    regex = regex.replace("(", "{").replace(")", "}")
    return regex


def parse_elm_motifs(elm_file):
    motifs = {}
    with open(elm_file, "r") as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            motifs[row["ELMIdentifier"]] = row["Regex"]
    return motifs


# step 2: Search motifs in protein sequences
def search_motifs(sequence, motifs):
    matches = []
    for motif_name, regex in motifs.items():
        for match in re.finditer(regex, sequence):
            matches.append((motif_name, match.start(), match.end()))
    return matches


# step 3: Main pipeline for motif analysis
def analyze_protein_family(fasta_file, prosite_file, elm_file):

    prosite_motifs = parse_prosite_motifs(prosite_file)
    elm_motifs = parse_elm_motifs(elm_file)
    all_motifs = {**prosite_motifs, **elm_motifs}

    sequences = {
        record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")
    }

    results = {}
    for protein_id, sequence in sequences.items():
        matches = search_motifs(sequence, all_motifs)
        results[protein_id] = matches

    return results


# step 4: Parse MobiDB disordered regions
def parse_mobidb_disordered_regions(mobidb_file):
    disordered_regions = {}
    with open(mobidb_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            protein_id = row[0]
            regions = json.loads(row[1])  # Parse the JSON string
            disordered_regions[protein_id] = regions
    return disordered_regions


# step 5: Check if a motif is within disordered regions
def is_within_disordered_regions(start, end, disordered_regions):
    for region_start, region_end in disordered_regions:
        if start >= region_start and end <= region_end:
            return True
    return False


# step 6: Cross-reference motifs with disordered regions
def filter_true_matches(motif_results, mobidb_disordered_regions):
    true_matches = {}
    for protein_id, matches in motif_results.items():
        if protein_id in mobidb_disordered_regions:
            disordered_regions = mobidb_disordered_regions[protein_id]
            true_matches[protein_id] = [
                match
                for match in matches
                if is_within_disordered_regions(
                    match[1] + 1, match[2], disordered_regions
                )
            ]
    return true_matches


def get_common_motifs(matches):

    # Step 1: Count motif occurrences across proteins
    motif_count = defaultdict(int)
    total_proteins = len(true_matches)

    for matches in true_matches.values():
        # Use a set to avoid counting duplicates within the same protein
        unique_motifs = set(match[0] for match in matches)
        for motif in unique_motifs:
            motif_count[motif] += 1

    # Step 2: Identify motifs present in all proteins
    shared_motifs = [
        motif for motif, count in motif_count.items() if count == total_proteins
    ]

    return shared_motifs


fasta_file = "PF03060_reviewed.fasta"
prosite_file = "prosite_motifs.txt"
elm_file = "elm_classes.tsv"
mobidb_file = "mobidb_lite_swissprot.csv"


motif_results = analyze_protein_family(fasta_file, prosite_file, elm_file)

with open("motif_matches.txt", "w") as file:
    for protein_id, matches in motif_results.items():
        file.write(f"Protein: {protein_id}\n")
        for motif_name, start, end in matches:
            file.write(f"  Motif: {motif_name}, Start: {start + 1}, End: {end}\n")

mobidb_disordered_regions = parse_mobidb_disordered_regions(mobidb_file)
true_matches = filter_true_matches(motif_results, mobidb_disordered_regions)

with open("true_motif_matches.txt", "w") as file:
    for protein_id, matches in true_matches.items():
        file.write(f"Protein: {protein_id}\n")
        for motif_name, start, end in matches:
            file.write(f"  Motif: {motif_name}, Start: {start + 1}, End: {end}\n")


common_motifs = get_common_motifs(motif_results)

with open("conserved_motifs_in_family.txt", "w") as file:
    for motif in common_motifs:
        print(motif)
        file.write(f"{motif}\n")
