---
outline: deep
---

# Architecture 

<p align="center">
  <img align="center" src="/architecture.png" alt="Logo" width="90%">
</p>


<p align="center">
  <img align="center" src="/abstract.png" alt="Logo" width="90%">
</p>


## Backend Architecture
- Written in **Python** using the **FastAPI** framework.
- Provides **RESTful APIs** for cheminformatics operations.
- **Hybrid toolkit integration**:
  - Native Python toolkits: **RDKit**, **Open Babel**.
  - Java toolkits (e.g., **CDK**, **SRU**, **OPSIN**) accessed via **JPype**.
- **Docker containerization** ensures reproducibility and consistent deployment.
- Designed for **scalability** and future extensibility without disrupting current services.
- Publicly available at https://api.naturalproducts.net/latest/docs.

## Frontend Architecture
- Built using **React** for strong community support and sustainability.
- **Component-based design** with modular service layers and **custom hooks**.
- Styled using **Tailwind CSS** for responsive UI.
- Communicates with backend via RESTful API using **Axios**.
- API interactions are abstracted through a **dedicated service layer**.
- Structured into **functionally specialized pages**.
- Publicly available at **https://app.naturalproducts.net/**.

## Toolkit Integration Challenge
- Toolkits use different languages:  
  - **RDKit/OpenBabel**: C++/Python  
  - **CDK**: Java
- Integrating multiple toolkits is complex due to **compatibility and environment setup**.
- Solution: use **containerization** and **microservices** for abstraction and flexibility.

## Microservices & Containerization
- Microservices: each function is a **separate, scalable, maintainable unit**.
- Containers provide isolated, reproducible environments.
- **Docker** is used for packaging and sharing toolkits.
- CM (Cheminformatics Microservice) is **publicly available via Docker Hub**.

## API Design
- Uses **REST API** based on **OpenAPI Specification 3.1.0**.
- Built with **FastAPI** for speed, efficiency, and automatic documentation.
- Promotes **interoperability**, **validation**, and **tool integration**.

## Functional Scope
- Provides:
  - **Format conversions**
  - **Optical Structure Recognition (OSR)**
  - **Chemical data standardization**
  - **Descriptor calculations**
- Integrated tools:
  - Cheminformatics: **RDKit**, **CDK**, **OpenBabel**
  - Deep learning: **DECIMER**, **STOUT**
- Enables **scalable**, **interoperable**, and **efficient** cheminformatics applications.

## Architectural Decisions
- Combines **FastAPI** + **Docker** for seamless deployment.
- Maintains flexibility via **microservice architecture**:
  - Individual services can be updated without affecting others.
- Contrary to usual practice, **all toolkits are packaged in one container**:
  - Simplifies deployment and avoids **complex orchestration**.

## Deployment & Monitoring
- **Dockerfile**, **docker-compose YAML**, and deployment scripts are available on GitHub.
- **HELM charts** provided for Kubernetes-based deployments.
- **Monitoring and visualization**:
  - **Prometheus** for metrics collection and alerting.
  - **Grafana** for real-time dashboards and usage visualization.






