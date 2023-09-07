---
outline: deep
---

# Docker
The Cheminformatics Microservice project utilizes the containerized microservices approach to package chemistry toolkits and deep learning tools. It comes pre-packaged with toolkits [RDKit](https://github.com/rdkit/rdkit), [CDK](https://doi.org/10.1186/s13321-017-0220-4), [OpenBabel](https://github.com/openbabel/openbabel) and deep learning tools ([DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer), [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator)) for handling chemical data - OSR, format conversions, and descriptor calculation. This enables efficient handling of large data volumes and improved performance and development of cheminformatics applications that are scalable and interoperable.

It is containerized using Docker and is distributed publicly via the [Docker Hub](https://hub.docker.com/r/nfdi4chem/cheminformatics-microservice), a cloud-based registry provided by Docker that allows developers to store, share, and distribute Docker images.

To use this image:

* Make sure you have Docker installed and configured on your target deployment environment.
* Pull the image by providing the appropriate tag.

```bash
docker pull nfdi4chem/cheminformatics-microservice:[tag]

```
* Run the below command to run the image

```bash
docker run -d -p 8080:80 --name [image-name] nfdi4chem/cheminformatics-microservice:[tag]

```

## Docker Compose

[Docker Compose](https://docs.docker.com/get-started/08_using_compose/) is a handy tool that allows you to define and manage multi-container applications. By creating a YAML file to define the services, you can use a simple command to start or stop everything at once.
Check out https://docs.docker.com/get-started/08_using_compose/ for more information.

The Cheminformatics Microservice is readily equipped with a docker-compose file that enables you to effortlessly deploy your application on your server. To deploy your application using Docker compose, simply follow the steps described below.The Cheminformatics Microservice comes packaged with a docker-compose file, which you can use to deploy your application on your server. To deploy using Docker compose, follow the steps described below.

### Installation
1. Before you run the command make sure you've installed Docker including docker-compose support. More details [here](https://docs.docker.com/compose/install/).
2. Clone the repository
```bash
git clone https://github.com/Steinbeck-Lab/cheminformatics-microservice.git
cd cheminformatics-microservice
```
3. Start the application
```bash
docker-compose -f ops/docker-compose-prod.yml up -d
```

## Scaling

For better performance you can scale your application using docker-compose builtin scaling support. Requests are load balanced by [Traefik](https://doc.traefik.io/traefik/). For example to scale upto 3 additional application containers you can simply invoke:
```bash
docker-compose -f ops/docker-compose-prod.yml up -d --scale web=3 --no-recreate
```