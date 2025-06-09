FROM continuumio/miniconda3:24.1.2-0 AS cheminf-python-ms

ENV PYTHON_VERSION=3.11 \
    INCLUDE_OCSR=true \
    JAVA_HOME=/usr/lib/jvm/java-11-openjdk-arm64/ \
    # Add default number of workers
    WORKERS=2 \
    # Add other Python configurations
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        openjdk-11-jre \
        curl \
        build-essential \
        gcc \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget -O /usr/bin/surge "https://github.com/StructureGenerator/surge/releases/download/v1.0/surge-linux-v1.0" && \
    chmod +x /usr/bin/surge

# Combine conda and pip operations to reduce layers
WORKDIR /code
COPY requirements.txt .
RUN conda install -c conda-forge python=${PYTHON_VERSION} sqlite --force-reinstall && \
    python3 -m pip install --no-cache-dir -U pip setuptools && \
    pip3 install --no-cache-dir -r requirements.txt && \
    # Install specific packages without dependencies
    pip3 install --no-cache-dir --no-deps \
        git+https://github.com/Kohulan/DECIMER-Image-Segmentation.git@bbox \
        decimer==2.7.1 \
        chembl_structure_pipeline


COPY ./app ./app

CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port 80 --workers ${WORKERS}"]