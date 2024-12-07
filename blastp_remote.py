from Bio.Blast.Applications import NcbiblastpCommandline

fasta_file = "Q12723.fasta"
output_file = "blastp_output.xml"

blastp_cline = NcbiblastpCommandline(query=fasta_file, db="nr", evalue=0.001, outfmt=5, out=output_file, remote=True)
stdout, stderr = blastp_cline()

print(f"BLASTP search complete. Results saved in {output_file}")
