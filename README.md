<p align="center"><a href="https://api.naturalproducts.net/" target="_blank"><img src="/public/img/logo.png" width="400" alt="CMS Logo"></a></p>

[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/Steinbeck-Lab/cheminformatics-microservice/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/Steinbeck-Lab/cheminformatics-microservice.svg)](https://GitHub.com/Steinbeck-Lab/cheminformatics-microservice/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/Steinbeck-Lab/cheminformatics-microservice.svg)](https://GitHub.com/Steinbeck-Lab/cheminformatics-microservice/graphs/contributors/)
[![tensorflow](https://img.shields.io/badge/TensorFlow-2.10.1-FF6F00.svg?style=flat&logo=tensorflow)](https://www.tensorflow.org)
[![Powered by CDK](https://img.shields.io/badge/Powered%20by-CDK-blue.svg?style=flat&logo=chem)](https://cdk.github.io)
[![RDKit badge](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
![Workflow](https://github.com/Steinbeck-Lab/cheminformatics-microservice/actions/workflows/dev-build.yml/badge.svg)
![Workflow](https://github.com/Steinbeck-Lab/cheminformatics-microservice/actions/workflows/prod-build.yml/badge.svg)
![Workflow](https://github.com/Steinbeck-Lab/cheminformatics-microservice/actions/workflows/release-please.yml/badge.svg)
[![framework](https://img.shields.io/badge/Framework-FastAPI-blue?style)](https://fastapi.tiangolo.com/)
[![FastAPI Documentation](https://img.shields.io/badge/docs-fastapi-blue)](https://api.naturalproducts.net/v1/docs#/)
[![Documentation Status](https://readthedocs.org/projects/cheminformatics-microservice/badge/?version=latest)](https://cheminformatics-microservice.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/Steinbeck-Lab/cheminformatics-microservice/graph/badge.svg?token=5BIQJPNCBA)](https://codecov.io/gh/Steinbeck-Lab/cheminformatics-microservice)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7745987.svg)](https://doi.org/10.5281/zenodo.7745987)
## Overview of Cheminformatics Microservice :globe_with_meridians:

This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator) and another instance of [DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer) (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).

## Documentation :book:

https://docs.api.naturalproducts.net/

### API Swagger Docs

Production: https://api.naturalproducts.net/latest/docs

Development: https://dev.api.naturalproducts.net/latest/docs

## Installation :computer:

You can run Cheminformatics Microservice as a standalone application using Python virtual environment or via Docker, or deploy to a Kubernetes cluster utilising [Helm charts](https://github.com/NFDI4Chem/repo-helm-charts/tree/main/charts). Please follow the links below for step-by-step instructions.

**Docker**
https://docs.api.naturalproducts.net/docker.html

**Kubernetes - Cluster deployment**
https://docs.api.naturalproducts.net/cluster-deployment.html

**Standalone**
https://docs.api.naturalproducts.net/standalone.html

## Load ramping test results - 29 September 2023 :chart_with_upwards_trend:
 https://github.com/Steinbeck-Lab/cheminformatics-microservice/discussions/413

## License :scroll:

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Steinbeck-Lab/cheminformatics-microservice/blob/main/LICENSE) file for details

## Citation :newspaper:

*Paper*
- Chandrasekhar, V., Sharma, N., Schaub, J. et al. Cheminformatics Microservice: unifying access to open cheminformatics toolkits. J Cheminform 15, 98 (2023). https://doi.org/10.1186/s13321-023-00762-4

*Software*
- Venkata, C., Sharma, N., & Rajan, K. (2023). Cheminformatics Microservice (Version v1.6.0) [Computer software]. https://doi.org/10.5281/zenodo.8336440

## Maintained by :wrench:

Cheminformatics Microservice and [Natural Products Online](https://naturalproducts.net) are developed and maintained by the [Steinbeck group](https://cheminf.uni-jena.de) at the [Friedrich Schiller University](https://www.uni-jena.de/en/) Jena, Germany.
The code for this web application is released under the [MIT license](https://opensource.org/licenses/MIT). Copyright © CC-BY-SA 2023
<p align="center"><a href="https://cheminf.uni-jena.de/" target="_blank"><img src="https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png" width="800" alt="cheminf Logo"></a></p>

## Acknowledgments :bulb:

Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Project number: **441958208** and  [ChemBioSys](https://www.chembiosys.de/en/) (Project INF) - Project number: **239748522 - SFB 1127**.

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
