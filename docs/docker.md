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
docker run -d -p 8080:80 --name [image-name] caffeinejena/cheminfo-microservice:[tag]

```

# Docker Compose
[Docker Compose](https://docs.docker.com/get-started/08_using_compose/) is a tool that helps you define and share multi-container applications. With Compose, you can create a YAML file to define the services and with a single command, you can spin everything up or tear it all down. 
The Cheminformatics Python Microservice comes packaged with docker-compose file which you can use to deploy your application in your server. To deploy using Docker compose follow the steps as described below.

### Installation
1. Before you run the command make sure you've installed Docker including docker-compose support. If not then follow the link [here](https://docs.docker.com/compose/install/).
2. Clone the repository
```bash
git clone https://github.com/Steinbeck-Lab/cheminformatics-python-microservice.git
cd cheminformatics-python-microservice
```
3. For local deployment 
```bash
docker-compose up -d
```
4. For production deployment
```bash
docker-compose -f ops/docker-compose-prod.yml up -d
```
5. Scaling in case of performance issues - This service file supports the docker-compose builtin scaling load balanced by [Traefik](https://doc.traefik.io/traefik/). For example to add 3 additional application containers you can simply invoke:
```bash
docker-compose -f ops/docker-compose-prod.yml up -d --scale web=3 --no-recreate
```