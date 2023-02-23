FROM debian:buster-slim AS rdkit-build-env

RUN apt-get update \
 && apt-get install -yq --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    wget \
    python3-dev \
    python3-setuptools \
    python-rdkit \
    librdkit1 \
    rdkit-data \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-python-dev \
    libboost-regex-dev \
    libboost-serialization-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libcairo2-dev \
    libeigen3-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

ARG RDKIT_VERSION=Release_2022_09_4
RUN wget --quiet https://github.com/rdkit/rdkit/archive/${RDKIT_VERSION}.tar.gz \
 && tar -xzf ${RDKIT_VERSION}.tar.gz \
 && mv rdkit-${RDKIT_VERSION} rdkit \
 && rm ${RDKIT_VERSION}.tar.gz

RUN mkdir /rdkit/build
WORKDIR /rdkit/build

RUN cmake \
    -D PYTHON_EXECUTABLE=/usr/bin/python3 \
    -D CMAKE_BUILD_TYPE=Release \
    -D CMAKE_INSTALL_PREFIX=/usr \
    -D RDK_BUILD_AVALON_SUPPORT=ON \
    -D RDK_BUILD_CAIRO_SUPPORT=ON \
    -D RDK_BUILD_CPP_TESTS=OFF \
    -D RDK_BUILD_INCHI_SUPPORT=ON \
    -D RDK_BUILD_FREESASA_SUPPORT=ON \
    -D RDK_INSTALL_INTREE=OFF \
    -D RDK_INSTALL_STATIC_LIBS=OFF \
    -D Boost_USE_STATIC_LIBS=OFF \
    ..

RUN make \
 && make install

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

# Copy rdkit installation from rdkit-build-env
COPY --from=rdkit-build-env /usr/lib/libRDKit* /usr/lib/
COPY --from=rdkit-build-env /usr/lib/cmake/rdkit/* /usr/lib/cmake/rdkit/
COPY --from=rdkit-build-env /usr/share/RDKit /usr/share/RDKit
COPY --from=rdkit-build-env /usr/include/rdkit /usr/include/rdkit
COPY --from=rdkit-build-env /usr/lib/python3/dist-packages/rdkit /usr/lib/python3.7/dist-packages/rdkit

WORKDIR /code

COPY ./requirements.txt /code/requirements.txt

RUN pip3 install --no-cache-dir --upgrade -r /code/requirements.txt

RUN pip3 install chembl_structure_pipeline --no-deps

COPY ./app /code/app

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80", "--reload"]