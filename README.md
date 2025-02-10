# Cheminformatics Microservice

<p align="center">
  <a href="https://api.naturalproducts.net/" target="_blank">
    <img src="/public/img/logo.png" width="400" alt="CMS Logo">
  </a>
</p>

<p align="center">
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT%202.0-blue.svg" alt="License"></a>
  <a href="https://GitHub.com/Steinbeck-Lab/cheminformatics-microservice/graphs/commit-activity"><img src="https://img.shields.io/badge/Maintained%3F-yes-blue.svg" alt="Maintenance"></a>
  <a href="https://GitHub.com/Steinbeck-Lab/cheminformatics-microservice/issues/"><img src="https://img.shields.io/github/issues/Steinbeck-Lab/cheminformatics-microservice.svg" alt="GitHub issues"></a>
  <a href="https://GitHub.com/Steinbeck-Lab/cheminformatics-microservice/graphs/contributors/"><img src="https://img.shields.io/github/contributors/Steinbeck-Lab/cheminformatics-microservice.svg" alt="GitHub contributors"></a>
</p>

<p align="center">
  <a href="https://www.tensorflow.org"><img src="https://img.shields.io/badge/TensorFlow-2.10.1-FF6F00.svg?style=flat&logo=tensorflow" alt="tensorflow"></a>
  <a href="https://cdk.github.io"><img src="https://img.shields.io/badge/Powered%20by-CDK-blue.svg?style=flat&logo=chem" alt="Powered by CDK"></a>
  <a href="https://www.rdkit.org/"><img src="https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC" alt="RDKit badge"></a>
</p>

<p align="center">
  <img src="https://github.com/Steinbeck-Lab/cheminformatics-microservice/actions/workflows/dev-build.yml/badge.svg" alt="Dev Build">
  <img src="https://github.com/Steinbeck-Lab/cheminformatics-microservice/actions/workflows/prod-build.yml/badge.svg" alt="Prod Build">
  <img src="https://github.com/Steinbeck-Lab/cheminformatics-microservice/actions/workflows/release-please.yml/badge.svg" alt="Release Please">
</p>

<p align="center">
  <a href="https://fastapi.tiangolo.com/"><img src="https://img.shields.io/badge/Framework-FastAPI-blue?style" alt="framework"></a>
  <a href="https://api.naturalproducts.net/v1/docs#/"><img src="https://img.shields.io/badge/docs-fastapi-blue" alt="FastAPI Documentation"></a>
  <a href="https://cheminformatics-microservice.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/cheminformatics-microservice/badge/?version=latest" alt="Documentation Status"></a>
  <a href="https://codecov.io/gh/Steinbeck-Lab/cheminformatics-microservice"><img src="https://codecov.io/gh/Steinbeck-Lab/cheminformatics-microservice/graph/badge.svg?token=5BIQJPNCBA" alt="codecov"></a>
  <a href="https://doi.org/10.5281/zenodo.7745987"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7745987.svg" alt="DOI"></a>
</p>

## 🌐 Overview

A set of microservices to support cheminformatics through API calls. It primarily works with SMILES-based inputs and offers functionalities such as:

- Translating between different machine-readable representations
- Calculating Natural Product (NP) likeliness scores
- Visualizing chemical structures
- Generating descriptors

Additionally, it hosts instances of:
- [DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer): A deep learning model for optical chemical structure recognition
> Important Note: STOUT (Smiles-TO-iUpac-Translator) is currently unavailable due to legal issues.

## 📚 Documentation

Comprehensive documentation is available at: [https://docs.api.naturalproducts.net/](https://docs.api.naturalproducts.net/)

### API Swagger Docs

- Production: [https://api.naturalproducts.net/latest/docs](https://api.naturalproducts.net/latest/docs)
- Development: [https://dev.api.naturalproducts.net/latest/docs](https://dev.api.naturalproducts.net/latest/docs)

## 💻 Installation

You can run Cheminformatics Microservice in multiple ways:

1. As a standalone application using a Python virtual environment
2. Via Docker
3. Deploy to a Kubernetes cluster using [Helm charts](https://github.com/NFDI4Chem/repo-helm-charts/tree/main/charts)

For detailed instructions, please refer to:

- [Docker Installation Guide](https://docs.api.naturalproducts.net/docker.html)
- [Kubernetes Cluster Deployment Guide](https://docs.api.naturalproducts.net/cluster-deployment.html)
- [Standalone Installation Guide](https://docs.api.naturalproducts.net/standalone.html)

## 📈 Load Ramping Test Results

View the latest load ramping test results from September 29, 2023: [Test Results Discussion](https://github.com/Steinbeck-Lab/cheminformatics-microservice/discussions/413)

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](https://github.com/Steinbeck-Lab/cheminformatics-microservice/blob/main/LICENSE) file for details.

## 📰 Citation

### Paper
Chandrasekhar, V., Sharma, N., Schaub, J. et al. Cheminformatics Microservice: unifying access to open cheminformatics toolkits. J Cheminform 15, 98 (2023). [https://doi.org/10.1186/s13321-023-00762-4](https://doi.org/10.1186/s13321-023-00762-4)

### Software
Venkata, C., Sharma, N., & Rajan, K. (2023). Cheminformatics Microservice (Version v2.6.0) [Computer software]. [https://zenodo.org/records/13867839](https://zenodo.org/records/13867839)

## 🔧 Maintenance

Cheminformatics Microservice and [Natural Products Online](https://naturalproducts.net) are developed and maintained by the [Steinbeck group](https://cheminf.uni-jena.de) at the [Friedrich Schiller University](https://www.uni-jena.de/en/) Jena, Germany.

The code for this web application is released under the [MIT license](https://opensource.org/licenses/MIT). Copyright © CC-BY-SA 2023

<p align="center">
  <a href="https://cheminf.uni-jena.de/" target="_blank">
    <img src="https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png" width="800" alt="cheminf Logo">
  </a>
</p>

## 💡 Acknowledgments

Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Project number: **441958208** and [ChemBioSys](https://www.chembiosys.de/en/) (Project INF) - Project number: **239748522 - SFB 1127**.

<div style="display: flex; justify-content: space-between;">
  <a href="https://www.dfg.de/" target="_blank">
    <img src="./docs/public/dfg_logo_schriftzug_blau_foerderung_en.gif" width="30%" alt="DFG Logo">
  </a>
  <a href="https://nfdi4chem.de/" target="_blank">
    <img src="https://www.nfdi4chem.de/wp-content/themes/wptheme/assets/img/logo.svg" width="30%" alt="NFDI4Chem Logo">
  </a>
  <a href="https://www.chembiosys.de/en/welcome.html" target="_blank">
    <img src="https://github.com/Steinbeck-Lab/cheminformatics-microservice/assets/30716951/45c8e153-8322-4563-a51d-cbdbe4e08627" width="30%" alt="Chembiosys Logo">
  </a>
</div>
