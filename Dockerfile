FROM continuumio/miniconda3:24.1.2-0 AS cheminf-python-ms

ENV PYTHON_VERSION=3.10
ENV RDKIT_VERSION=2023.09.4
ENV OPENBABEL_VERSION=v3.1.1

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    apt-get update -y && \
    apt-get install -y openjdk-11-jre && \
    apt-get install -y curl && \
    conda update -n base -c defaults conda

RUN wget -O surge "https://github.com/StructureGenerator/surge/releases/download/v1.0/surge-linux-v1.0"
RUN chmod +x surge
RUN mv surge /usr/bin

RUN conda install -c conda-forge python>=PYTHON_VERSION
#RUN conda install -c conda-forge rdkit==RDKIT_VERSION
RUN conda install -c conda-forge openbabel>=OPENBABEL_VERSION

RUN python3 -m pip install -U pip

ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME

WORKDIR /code
COPY ./requirements.txt /code/requirements.txt

RUN pip3 install --upgrade setuptools pip
RUN pip3 install --no-cache-dir -r /code/requirements.txt
RUN python3 -m pip uninstall -y imantics
RUN pip3 install imantics==0.1.12
RUN pip3 install rdkit
RUN pip3 install --no-deps decimer-segmentation==1.1.3
RUN pip3 install --no-deps decimer==2.3.0
RUN pip3 install --no-deps STOUT-pypi>=2.0.5
RUN python3 -m pip install uvicorn[standard]

RUN pip3 install --no-cache-dir chembl_structure_pipeline --no-deps

COPY ./app /code/app

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--workers", "2"]
