---
outline: deep
---

# Versions

<div style="text-align: justify;">

To ensure accurate data reproducibity, it is essential to have a strong tool versioning system in place. The best practices for research data management recommend documenting the software and its components' versions. This is especially important for tools like CM, which involve multiple dependencies. CM uses multi-level versioning to document the API and underlying software dependencies. This approach is user-friendly and avoids overwhelming users with confusing versions.

The CM codebase is updated twice a year (bi-annual release cycle) with documentation for the underlying toolkits, tools, and environment dependencies for each release. The REST API also has release cycles that run parallel to the software release cycle, allowing you to introduce changes and enhancements to your API without breaking existing client applications. However, the REST API releases are only based on changes in REST communication and will not be updated if there are no changes in REST endpoints. In principle, researchers can update the underlying cheminformatics toolkits as and when new releases are available without updating their code bases since the REST API would remain the same.

The Cheminformatics Microservice, framework is developed using Python and [FastAPI](https://fastapi.tiangolo.com/), leveraging the [RDKit](https://www.rdkit.org/), [OpenBabel](http://openbabel.org/wiki/Main_Page) and [Chemistry Development Kit (CDK)](https://cdk.github.io/) libraries in the backend. CDK is accessed through [jpype](https://jpype.readthedocs.io/en/latest/index.html), allowing the implementation of functions that can be seamlessly utilized within the Python environment. For detailed information about the specific versions employed, please refer to the list provided below.
</div>

## v1.0.0

**Cheminformatics toolkits Versions**

| Toolkits    | Version     |
|-------------|-------------|
| RDKit <sup>1</sup>    | 2023.09.4  |
| CDK <sup>2</sup>      | 2.9.0       |
| Open Babel <sup>3</sup> | 3.1.1       |


**External Tools Versions**

| Toolkits                  | Version   |
|---------------------------|-----------|
| STOUT<sup>4</sup>         | 2.0.5     |
| DECIMER<sup>5,6</sup>     | 2.3.0     |
| DECIMER-Segmentation<sup>7</sup> | 2.3.0     |
| Surge<sup>8</sup>         | 1.3.2     |
| SRU<sup>9</sup>           | 1.3.2     |
| ClassyFire<sup>10</sup>   | 1.0.0     |


## References

1. Landrum G, Others (2016) RDKit: Open-Source Cheminformatics Software.(2016). URL http://www.rdkit.org/, https://github.com/rdkit/rdkit
2. Willighagen EL, Mayfield JW, Alvarsson J, et al (2017) The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. Journal of Cheminformatics. https://doi.org/10.1186/s13321-017-0220-4
3. O'Boyle, N.M., Banck, M., James, C.A. et al. Open Babel: An open chemical toolbox. J Cheminform 3, 33 (2011). https://doi.org/10.1186/1758-2946-3-33
4. Rajan, K., Zielesny, A. & Steinbeck, C. STOUT: SMILES to IUPAC names using neural machine translation. J Cheminform 13, 34 (2021). https://doi.org/10.1186/s13321-021-00512-4
5. Rajan, K., Zielesny, A. & Steinbeck, C. DECIMER: towards deep learning for chemical image recognition. J Cheminform 12, 65 (2020). https://doi.org/10.1186/s13321-020-00469-w
6. Rajan K, Brinkhaus HO, Agea MI, Zielesny A, Steinbeck C (2023) DECIMER.ai - An open platform for automated optical chemical structure identification, segmentation and recognition in scientific publications. ChemRxiv (2023). https://doi.org/10.26434/chemrxiv-2023-xhcx9
7. Rajan, K., Brinkhaus, H.O., Sorokina, M. et al. DECIMER-Segmentation: Automated extraction of chemical structure depictions from scientific literature. J Cheminform 13, 20 (2021). https://doi.org/10.1186/s13321-021-00496-1
8. McKay, B.D., Yirik, M.A. & Steinbeck, C. Surge: a fast open-source chemical graph generator. J Cheminform 14, 24 (2022). https://doi.org/10.1186/s13321-022-00604-9
9. Schaub, J., Zielesny, A., Steinbeck, C. et al. Too sweet: cheminformatics for deglycosylation in natural products. J Cheminform 12, 67 (2020). https://doi.org/10.1186/s13321-020-00467-y
10. Djoumbou Feunang, Y., Eisner, R., Knox, C. et al. ClassyFire: automated chemical classification with a comprehensive, computable taxonomy. J Cheminform 8, 61 (2016). https://doi.org/10.1186/s13321-016-0174-y
