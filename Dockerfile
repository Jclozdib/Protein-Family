FROM python:3.8-slim

WORKDIR /app

RUN apt-get update && \
    apt-get install -y \
    wget \
    build-essential \
    libz-dev \
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

RUN apt-get update && apt-get install -y clustalw

ENV PATH="/app/ncbi-blast-2.16.0+/bin:${PATH}"

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . /app

COPY workflow.sh /app/
RUN chmod +x /app/workflow.sh
