---
outline: deep
---

# Docker
The Cheminformatics Python Microservice project utilizes the containerized microservices approach to package chemistry toolkits and deep learning tools. It comes pre-packaged with toolkits [RDKit](https://github.com/rdkit/rdkit), [CDK](https://doi.org/10.1186/s13321-017-0220-4), [OpenBabel](https://github.com/openbabel/openbabel) and deep learning tools ([DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer), [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator)) for handling chemical data - OSR, format conversions, and descriptor calculation. This enables efficient handling of large data volumes and improved performance and development of cheminformatics applications that are scalable and interoperable.

It is containerized using Docker and is distributed publicly via the [Docker Hub](https://hub.docker.com/r/caffeinejena/cheminformatics-python-microservice), a cloud-based registry provided by Docker that allows developers to store, share, and distribute Docker images.

To use this image:

* Make sure you have Docker installed and configured on your target deployment environment.
* Pull the image by providing the appropiate tag.

```bash
docker pull caffeinejena/cheminformatics-python-microservice:[tag]

```
* Run the below command to run the image

```bash
docker run -d -p 8080:80 --build-arg RELEASE_VERSION= [release-version] --name [image-name] caffeinejena/cheminfo-microservice:[tag]

```