FROM --platform=linux/amd64 python:3.12-slim

WORKDIR /app

RUN apt-get update && \
    apt-get install -y \
    gzip \
    wget \
    build-essential \
    libz-dev \
    graphviz \
    graphviz-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
    tar -xzf hmmer.tar.gz && \
    cd hmmer-3.4 && \
    ./configure && \
    make && \
    cd .. && \
    rm -rf hmmer.tar.gz

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    rm ncbi-blast-2.16.0+-x64-linux.tar.gz

RUN wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    gunzip uniprot_sprot.fasta.gz

RUN apt-get update && apt-get install -y clustalw

ENV PATH="/app/ncbi-blast-2.16.0+/bin:${PATH}"

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . /app

COPY workflow.sh /app/
RUN chmod +x /app/workflow.sh
