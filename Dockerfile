FROM debian:buster-slim AS nmrxiv-python-ms

# Install runtime dependencies
RUN apt-get update \
 && apt-get install -yq --no-install-recommends \
    libboost-atomic1.67.0 \
    libboost-chrono1.67.0 \
    libboost-date-time1.67.0 \
    libboost-iostreams1.67.0 \
    libboost-python1.67.0 \
    libboost-regex1.67.0 \
    libboost-serialization1.67.0 \
    libboost-system1.67.0 \
    libboost-thread1.67.0 \
    libcairo2-dev \
    python3-dev \
    python3-numpy \
    python3-cairo \
    python3-pip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    apt-get update -y && \
    apt-get install -y openjdk-11-jre
    
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME

WORKDIR /code

RUN python3 -m pip install -U pip 

COPY ./requirements.txt /code/requirements.txt
RUN pip3 install --upgrade setuptools pip
RUN pip3 install --no-cache-dir --upgrade -r /code/requirements.txt
RUN pip3 install --no-deps decimer-segmentation
RUN pip3 install --no-deps decimer>=2.2.0
RUN pip3 install --no-deps STOUT-pypi>=2.0.5

RUN python3 -m pip uninstall -y uvicorn

RUN python3 -m pip install uvicorn[standard]

RUN pip3 install --no-cache-dir chembl_structure_pipeline --no-deps

COPY ./app /code/app

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--reload"]