FROM continuumio/miniconda3:24.11.1-0 AS cheminf-python-ms

ENV PYTHON_VERSION=3.11 \
    INCLUDE_OCSR=true \
    WORKERS=2 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        openjdk-17-jre-headless \
        curl \
        build-essential \
        wget \
        # OpenCV runtime dependencies (needed for OCSR/DECIMER)
        libgl1 \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    # Create arch-independent JAVA_HOME symlink
    ln -sf /usr/lib/jvm/java-17-openjdk-$(dpkg --print-architecture) /usr/lib/jvm/java-17-openjdk && \
    wget "https://github.com/StructureGenerator/surge/releases/download/v2.0/surge-linux-x86_64.tar.gz" && \
    tar xzf surge-linux-x86_64.tar.gz && \
    mv surge /usr/bin/surge && \
    chmod +x /usr/bin/surge && \
    rm surge-linux-x86_64.tar.gz

ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk/

ARG RELEASE_VERSION=1.0
ENV RELEASE_VERSION=${RELEASE_VERSION}

# Combine conda and pip operations to reduce layers
WORKDIR /code
COPY requirements.txt .
RUN conda install -c conda-forge python=${PYTHON_VERSION} sqlite --force-reinstall && \
    python3 -m pip install --no-cache-dir -U pip setuptools && \
    pip3 install --no-cache-dir -r requirements.txt && \
    # Install specific packages without dependencies
    pip3 install --no-cache-dir --no-deps \
        git+https://github.com/Kohulan/DECIMER-Image-Segmentation.git \
        decimer==2.7.1 \
        chembl_structure_pipeline


COPY ./app ./app

RUN useradd -m -r appuser && chown -R appuser:appuser /code
USER appuser

CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port 80 --workers ${WORKERS}"]
