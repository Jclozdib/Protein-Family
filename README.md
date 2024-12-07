# Protein-Family

## To be cleaned later:

1) Build and run the docker:
```sh
docker build -t protein-family .
docker run -it protein-family bash
```

2) Run the blastp search (already considers input sequence)
```sh
python blastp_remote.py
```

3) Run the multiple sequence alignment
```sh
python msa_clustalo.py
```

4) Insert command to perform PSSM (outputs both kinds of pssm files)
```sh
ncbi-blast-2.16.0+/bin/psiblast -subject Q12723.fasta -in_msa blastp_hits.fasta -out_ascii_pssm blastp_hits.pssm_ascii -out_pssm blastp_hits.pssm
```

5) Insert command to perform HMM
```sh
hmmer-3.4/src/hmmbuild hmmer_clustal.hmm blastp_hits.fasta
```