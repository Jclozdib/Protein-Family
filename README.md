# Protein-Family

## To be cleaned later:
## Make sh file with listed commands

1) Build and run the docker:
```sh
docker build -t protein-family .
docker run -d --name pf pf tail -f /dev/null
```

2) Run the blastp search (already considers input sequence)
```sh
docker exec python blastp_remote.py
```

3) Run the multiple sequence alignment
```sh
docker exec python msa_clustalo.py
```

4) Insert command to perform PSSM (outputs both kinds of pssm files)
```sh
docker exec ncbi-blast-2.16.0+/bin/psiblast -subject Q12723.fasta -in_msa blastp_hits.fasta -out_ascii_pssm blastp_hits.pssm_ascii -out_pssm blastp_hits.pssm
```

5) Insert command to perform HMM
```sh
docker exec hmmer-3.4/src/hmmbuild hmmer_clustal.hmm blastp_hits.fasta
```

6) Copy output files in local repository:
```sh
docker cp pf:/app/blastp_hits.pssm /desired/path
docker cp pf:/app/hmmer_clustalo.hmm /desired/path
```