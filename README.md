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
[![FastAPI Documentation](https://img.shields.io/badge/docs-fastapi-blue)](https://api.naturalproducts.net/docs#/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7747862.svg)](https://doi.org/10.5281/zenodo.7747862)
## Overview of Cheminformatics Python Microservices

This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of [STOUT](https://github.com/Kohulan/Smiles-TO-iUpac-Translator) and another instance of [DECIMER](https://github.com/Kohulan/DECIMER-Image_Transformer) (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).

## Example usage 

- Convertors

  - SMILES to IUPAC name
  ```fastapi
  https://api.naturalproducts.net/convert/iupac?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
  ```
  - SMILES to SELFIES
  ```fastapi
  https://api.naturalproducts.net/convert/selfies?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
  ```
  - SMILES to mol (default: CDK)
  ```fastapi
  https://api.naturalproducts.net/convert/mol?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
  ```
  ```fastapi
  https://api.naturalproducts.net/convert/mol?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&generator=rdkit
  ```

- Chem

  - Calculate Descriptors
  ```fastapi
  https://api.naturalproducts.net/chem/descriptors?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
  ```
  - Depict molecule (default: CDK)
  ```fastapi
  https://api.naturalproducts.net/chem/depict?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
  ```
  - Depict molecule with settings
  ```fastapi
  https://api.naturalproducts.net/chem/depict?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C&generator=rdkit&width=256&height=256&rotate=75
  ```
  - Visualize molecule in 3D
  ```fastapi
  https://api.naturalproducts.net/chem/depict3D?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C
  ```

> **Note**
> For detailed documentation on how to use the API check [here](https://api.naturalproducts.net/docs#/)

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/blob/dev-kohulan/LICENSE) file for details

## Citation

Venkata, C., Sharma, N., & Rajan, K. (2023). Cheminformatics Python Microservice (Version v0.8.0 - prerelease) [Computer software]. https://doi.org/10.5281/zenodo.7747862

## Version information
<p align="center">
  <b>Cheminformatics Python Microservice: V0.8.0</b>
</p>

<p align="center">
  
| Tools      | Version   |
|------------|-----------|
| RDKit      | 2023.03.1 |
| CDK        | 2.8.0     |
| STOUT      | 2.0.0     |
| DECIMER    | 2.2.0     |
  
</p>



## Maintained by

Cheminformatics Python Microservice and [Natural Products Online](https://naturalproducts.net) are developed and maintained by the [Steinbeck group](https://cheminf.uni-jena.de) at the [Friedrich Schiller University](https://www.uni-jena.de/en/) Jena, Germany. 
The code for this web application is released under the [MIT license](https://opensource.org/licenses/MIT). Copyright © CC-BY-SA 2023
<p align="center"><a href="https://cheminf.uni-jena.de/" target="_blank"><img src="https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png" width="800" alt="cheminf Logo"></a></p>
