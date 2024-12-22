#!/bin/bash

# Step 1: Build and run the Docker container
echo "Building the Docker image..."
docker build -t pf .

echo "Running the Docker container in detached mode..."
docker run -d --name pf pf tail -f /dev/null

# Step 2: Run the BLASTP search
echo "Running BLASTP search..."
# docker exec pf python blastp_remote.py
docker exec pf blastp -query base_sequence.fasta -db nr -evalue 0.0001 -out blastp_output.xml -outfmt 5 -remote

# Step 3: Run the multiple sequence alignment (Clustal Omega)
echo "Running multiple sequence alignment..."
docker exec pf python msa_clustalo.py

# Step 4: Perform PSI-BLAST to generate PSSM files
echo "Generating PSSM files using PSI-BLAST..."
docker exec pf ncbi-blast-2.16.0+/bin/psiblast -subject Q12723.fasta -in_msa msa.fasta -out_pssm model.pssm

# Step 5: Perform HMM construction
echo "Running HMMBUILD..."
docker exec pf hmmer-3.4/src/hmmbuild model.hmm msa.fasta

# Step 6: Run HMM-SEARCH against Swissprot
echo "Running HMM-SEARCH agaisnt sp"
docker exec pf hmmer-3.4/src/hmmsearch --domtblout /app/hmmer_results_domain.txt --cpu 4 -E 1e-5 -o /app/hmmer_full_output.txt /app/model.hmm /app/uniprot_sprot.fasta

# Step 7: Run PSI-BLAST against Swissprot
echo "Running PSI-BLAST against sp"
docker exec pf makeblastdb -in /app/uniprot_sprot.fasta -dbtype prot -out /app/swissprot
docker exec pf psiblast -db /app/swissprot -in_pssm /app/model.pssm -out /app/psiblast_results_domain.txt -evalue 1e-5 -num_threads 4 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qen
d sstart send evalue bitscore"

# Step : Copy output files to local system
echo "Copying output files to local system..."
docker cp pf:/app/model.pssm /home/jclozdib/Protein-family/
docker cp pf:/app/model.hmm /home/jclozdib/Protein-family/
docker cp pf:/app/hmmer_results_domain.txt /home/jclozdib/Protein-family/
docker cp pf:/app/psiblast_results_domain.txt /home/jclozdib/Protein-family/

echo "Workflow completed successfully!"
