---
outline: deep
---

# Scaling 

Scalability in cheminformatics tools is crucial to ensure that researchers and scientists can process and analyze chemical data effectively, regardless of the size and complexity of their datasets and computational tasks. It involves considerations of both hardware and software architecture, as well as efficient algorithms and data management strategies.

## Use cases

### Molecular Diversity (Parallelisation)

Implementing parallelization in molecular diversity-related tasks can significantly enhance the efficiency and speed of cheminformatics and drug discovery workflows. It allows researchers to handle larger datasets, explore chemical space more comprehensively, and make more informed decisions in a shorter amount of time, ultimately facilitating the discovery of novel compounds and materials.

### DECIMER - OCSR / STOUT (Resourceloc intensive jobs)

Scaling neural network model predictions effectively involves a combination of model optimization, hardware utilization, and deployment strategies tailored to the specific use case and available computational resources. In the case of DECIMER - OCSR due to large model size, requesting and maintaining a large resource all the time is not effective both in terms of costs and environmental impact. Application needs to up scale on demand and downscale if the request load is low.

Here we presented only two use cases where scaling can impact significantly the time, resources and energy requirements, but we believe that Cheminformatics tasks can take full advantage of the hardware through optimisation and developmental strategies.

Currently CPM handles one molecule per request but in future releases we will also support batch processing in a single request enabling a larger through put. 

## Implementations

### Docker compose (manual)

Docker compose file distributed with CPM comes with a preconfigured [Traefik](https://github.com/traefik/traefik) serivce - a modern HTTP reverse proxy and load balancer that makes deploying microservices easy.

It uses round robin approach to distribute requests among the resources available.

Please follow the instructions [here](/docker.html#scaling) to configure docker compose deployments.

### Helm Chart (Auto scaling) - HPA

Helm is a package manager for Kubernetes that simplifies deploying, managing, and scaling applications. When it comes to scaling applications using Helm charts, there are several strategies and best practices to consider - Replica Scaling, Autoscaling, Horizontal Pod Autoscaler etc.

While each service provider (AWS, Azure, GCP) have their own configurations to configure scaling. Here we provide instructions on how to achieve scaling on a GKE deployment using helm chart.

### GKE

The Helm chart for CPM comes packaged with Horizontal Pod Autoscaler - more information [here](/cluster-deployment.html#scaling).





