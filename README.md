<p align="center"><a href="https://api.naturalproducts.net/" target="_blank"><img src="/public/img/logo.png" width="400" alt="CMS Logo"></a></p>

[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/Steinbeck-Lab/cheminformatics-python-microservice/graphs/commit-activity)
[![GitHub issues](https://img.shields.io/github/issues/Steinbeck-Lab/cheminformatics-python-microservice.svg)](https://GitHub.com/Steinbeck-Lab/cheminformatics-python-microservice/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/Steinbeck-Lab/cheminformatics-python-microservice.svg)](https://GitHub.com/Steinbeck-Lab/cheminformatics-python-microservice/graphs/contributors/)
[![tensorflow](https://img.shields.io/badge/TensorFlow-2.10.1-FF6F00.svg?style=flat&logo=tensorflow)](https://www.tensorflow.org)
[![Powered by CDK](https://img.shields.io/badge/Powered%20by-CDK-blue.svg?style=flat&logo=chem)](https://cdk.github.io)
[![RDKit badge](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
![Workflow](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/actions/workflows/dev-build.yml/badge.svg)
![Workflow](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/actions/workflows/prod-build.yml/badge.svg)
![Workflow](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/actions/workflows/release-please.yml/badge.svg)
[![framework](https://img.shields.io/badge/Framework-FastAPI-blue?style)](https://fastapi.tiangolo.com/)
[![FastAPI Documentation](https://img.shields.io/badge/docs-fastapi-blue)](https://api.naturalproducts.net/v1/docs#/)
[![Documentation Status](https://readthedocs.org/projects/cheminformatics-python-microservice/badge/?version=latest)](https://cheminformatics-python-microservice.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8112749.svg)](https://doi.org/10.5281/zenodo.8112749)
## Overview of Cheminformatics Python Microservice

This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator) and another instance of [DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer) (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).

## Documentation

https://steinbeck-lab.github.io/cheminformatics-python-microservice/

### API Swagger Docs:

Production: https://api.naturalproducts.net/latest/docs

Development: https://dev.api.naturalproducts.net/latest/docs

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/blob/dev-kohulan/LICENSE) file for details

## Citation

Venkata, C., Sharma, N., & Rajan, K. (2023). Cheminformatics Python Microservice (Version V1.0.0) [Computer software]. https://doi.org/10.5281/zenodo.7747862

## Maintained by

Cheminformatics Python Microservice and [Natural Products Online](https://naturalproducts.net) are developed and maintained by the [Steinbeck group](https://cheminf.uni-jena.de) at the [Friedrich Schiller University](https://www.uni-jena.de/en/) Jena, Germany. 
The code for this web application is released under the [MIT license](https://opensource.org/licenses/MIT). Copyright © CC-BY-SA 2023
<p align="center"><a href="https://cheminf.uni-jena.de/" target="_blank"><img src="https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png" width="800" alt="cheminf Logo"></a></p>

## Acknowledgments

Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Projektnummer **441958208**.

<p align="left"><a href="https://www.dfg.de/" target="_blank"><img src="./docs/public/dfg_logo_schriftzug_blau_foerderung_en.gif" width="50%" alt="DFG Logo"></a></p>