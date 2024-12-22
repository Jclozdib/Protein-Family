import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIXML

blast_results = "blastp_output.xml"
query_file = "Q12723.fasta"

# Retrieve BLAST hits
with open(blast_results) as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    sequences = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sequences.append(hsp.sbjct)

# Write the sequences to a new FASTA file
with open("msa.fasta", "w") as output_fasta:
    for idx, seq in enumerate(sequences):
        output_fasta.write(f">sequence_{idx}\n{seq}\n")

# Perform MSA using ClustalW with output in FASTA format
try:
    msa_command = ["clustalw", "-INFILE=msa.fasta", "-OUTFILE=msa.fasta", "-OUTPUT=FASTA"]
    result = subprocess.run(msa_command, capture_output=True, text=True)
    
    # Check the command output for success/failure
    if result.returncode == 0:
        print("MSA completed. Results saved in 'msa.fasta'")
    else:
        print(f"Error in MSA execution. stderr: {result.stderr}")

except Exception as e:
    print(f"An error occurred while running ClustalW: {e}")

