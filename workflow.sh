#!/bin/bash

# Step 1: Build and run the Docker container
echo "Building the Docker image..."
docker build -t pf .

echo "Running the Docker container in detached mode..."
docker run -d --name pf pf tail -f /dev/null

# Step 2: Run the BLASTP search
echo "Running BLASTP search..."
docker exec pf python blastp_remote.py

# Step 3: Run the multiple sequence alignment (Clustal Omega)
echo "Running multiple sequence alignment..."
docker exec pf python msa_clustalo.py

# Step 4: Perform PSI-BLAST to generate PSSM files
echo "Generating PSSM files using PSI-BLAST..."
docker exec pf ncbi-blast-2.16.0+/bin/psiblast -subject Q12723.fasta -in_msa blastp_hits.fasta -out_ascii_pssm blastp_hits.pssm_ascii -out_pssm blastp_hits.pssm

# Step 5: Perform HMM construction
echo "Running HMMBUILD..."
docker exec pf hmmer-3.4/src/hmmbuild hmmer_clustal.hmm blastp_hits.fasta

# Step 6: Copy output files to local system
echo "Copying output files to local system..."
docker cp pf:/app/blastp_hits.pssm /home/jclozdib/Protein-family/
docker cp pf:/app/hmmer_clustal.hmm /home/jclozdib/Protein-family/

echo "Workflow completed successfully!"
