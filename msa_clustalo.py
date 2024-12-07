from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline

blast_results = "blastp_output.xml"
record = SeqIO.read("Q12723.fasta", "fasta")

# Retrieve blastp hits
with open(blast_results) as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    sequences = []
    
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sequences.append(hsp.sbjct)

# Write the sequences to a new FASTA file
with open("blastp_hits.fasta", "w") as output_fasta:
    for idx, seq in enumerate(sequences):
        output_fasta.write(f">sequence_{idx}\n{seq}\n")

# Perform MSA using Clustal Omega
clustalw_cline = ClustalwCommandline("clustalw2", infile="blastp_hits.fasta")
stdout, stderr = clustalw_cline()

print("MSA completed. Results saved in 'blastp_hits.aln'")
