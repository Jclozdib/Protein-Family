# Protein-Family

This is an automated bioinformatics workflow designed for protein sequence analysis. The pipeline performs homology searches, multiple sequence alignment, motif discovery, taxonomic lineage retrieval, and functional annotation using various computational tools, including BLAST, HMMER, PSI-BLAST, and GO enrichment analysis.
This project enables the identification of conserved motifs, domain matches, and biological annotations, facilitating large-scale protein family analysis.

## Automatic workflow

1. Run the entire workflow:

```sh
./workflow.sh
```

## Manual step-by-step process

1. Build and run the docker:

```sh
docker build -t protein-family .
docker run -d --name pf pf tail -f /dev/null
```

2. Run the blastp search (already considers input sequence)

```sh
docker exec python blastp_remote.py
```

3. Run the multiple sequence alignment

```sh
docker exec python msa_clustalo.py
```

4. Insert command to perform PSSM (outputs both kinds of pssm files)

```sh
docker exec ncbi-blast-2.16.0+/bin/psiblast -subject Q12723.fasta -in_msa blastp_hits.fasta -out_ascii_pssm blastp_hits.pssm_ascii -out_pssm blastp_hits.pssm
```

5. Insert command to perform HMM

```sh
docker exec hmmer-3.4/src/hmmbuild hmmer_clustal.hmm blastp_hits.fasta
```

6. Run HMM-SEARCH against Swissprot

```sh
docker exec pf hmmer-3.4/src/hmmsearch --domtblout /app/hmmer_results_domain.txt --cpu 4 -E 1e-5 -o /app/hmmer_full_output.txt /app/model.hmm /app/uniprot_sprot.fasta
```

7. Run PSI-BLAST against Swissprot

```sh
docker exec pf makeblastdb -in /app/uniprot_sprot.fasta -dbtype prot -out /app/swissprot
docker exec pf psiblast -db /app/swissprot -in_pssm /app/model.pssm -out /app/psiblast_results_domain.txt -evalue 1e-5 -num_threads 4 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qen
d sstart send evalue bitscore"
```

8. Retrieve domain matches and compare with models

```sh
docker exec pf python comparison.py
```

9. Perform taxonomy lineage retreival:

```sh
docker exec python lineage.py
```

10. Plot and save taxonmic tree:

```sh
docker exec pf python tax_tree.py
```

11. Get GO annotations for each proteins

```sh
docker exec pf python go_annotation.py
```

12. Calculate enrichment for each GO term

```sh
docker exec pf python enrichment.py
```

13. Generate word cloud

```sh
docker exec pf python word_cloud.py
```

14. Get morifs

```sh
docker exec pf python motifs.py
```

15. Copy output files in local repository:

```sh
docker cp pf:/app/blastp_output.xml ./results/
docker cp pf:/app/msa.fasta ./results/
docker cp pf:/app/model.pssm ./results/
docker cp pf:/app/model.hmm ./results/
docker cp pf:/app/hmmer_results_domain.txt ./results/
docker cp pf:/app/psiblast_results_domain.txt ./results/
docker cp pf:/app/comparison_results.txt ./results/
docker cp pf:/app/ground_truth_accessions.txt ./results/
docker cp pf:/app/lineage.txt ./results/
docker cp pf:/app/taxonomic_tree.png ./results/
docker cp pf:/app/ground_truth_annotations.tsv ./results/
docker cp pf:/app/enrichment_results.tsv ./results/
docker cp pf:/app/word_cloud.png ./results/
docker cp pf:/app/motif_matches.txt ./results/
docker cp pf:/app/true_motif_matches.txt ./results/
docker cp pf:/app/conserved_motifs_in_family.txt ./results/
```

# Contributors

- José Carlos Lozano Dibildox
- Mohamed Louai Bouzaher
