# Python
FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies + BLAST/HMMER
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    curl \
    git \
    blast2 \
    hmmer \
    && rm -rf /var/lib/apt/lists/*

# Python libraries
RUN pip install --no-cache-dir \
    biopython \
    pandas \
    requests

RUN apt-get update && apt-get install -y git

RUN git clone https://github.com/your-username/Protein-family.git /app

# Working directory
WORKDIR /app

# Copy your local project files into the container's working directory
COPY . /app

# Set the default command to run your main Python script
CMD ["python", "main.py"]
