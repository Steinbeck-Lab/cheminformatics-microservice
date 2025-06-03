---
outline: deep
---

<p align="center">
  <img align="center" src="/logo_big.png" alt="Logo" style="filter: drop-shadow(0px 0px 10px rgba(0, 0, 0, 0.5));" width="400">
</p>

<br/>

<div style="text-align: justify;">
The Cheminformatics Microservice offers a collection of versatile functions accessible via REST endpoints that can handle chemical data and perform various cheminformatics tasks. These tasks include (but are not limited to) generating chemical structure depictions, 3D conformers, descriptors, IUPAC names, and converting machine-readable formats. Researchers and developers can effectively access open-source cheminformatics toolkits such as CDK, RDKit, and OpenBabel through this microservice and extend them easily.

This microservice packaged as a docker image (container) enables effortless deployment and scalability, making it suitable for academic, research and industry applications. Because of their modular nature, these microservices can be customized and combined to meet various Cheminformatics research and chemical data analysis needs.
</div>


## Front-end 

<div style="text-align: justify;">

The Cheminformatics Microservice features a modern React frontend that offers an intuitive user interface for a wide range of chemical informatics tasks. Key capabilities include chemical structure analysis and validation, format conversion (e.g., SMILES, InChI, IUPAC), and both 2D and interactive 3D molecular visualization. It also supports Optical Chemical Structure Recognition (OCSR), allowing conversion of chemical structure images into SMILES. Advanced analysis tools include structure validation, molecular descriptor calculation, stereoisomer generation, tautomer standardization, natural product likeness scoring, functional group detection, HOSE code generation, and Tanimoto similarity calculation. Users can perform batch operations and generate 2D/3D coordinates, making the platform robust for cheminformatics research and applications. The platform is publicly accessible at https://app.naturalproducts.net.

### Key Features

#### Chemical Analysis Tools
- Structure validation and standardization  
- Molecular descriptor calculation  
- Stereoisomer generation  
- Natural product likeness scoring  
- Functional group detection (**Ertl method**)  
- Tautomer standardization  
- HOSE code generation  
- Tanimoto similarity calculation  

#### Format Conversion
- **SMILES** to various chemical formats (InChI, InChIKey, SELFIES, etc.)  
- **IUPAC name** to SMILES conversion  
- 2D/3D coordinate generation (Molblock)  

#### Molecule Depiction
- 2D molecule visualization with customizable settings  
- Interactive 3D molecule visualization  
- Batch depiction for multiple molecules  

#### OCSR (Optical Chemical Structure Recognition)
- Convert chemical structure images to **SMILES** notation


### Interactive Examples
Visit the [public API documentation](https://api.naturalproducts.net/latest/docs) to explore interactive examples and test the endpoints directly in your browser.


<div style="text-align: justify;">

::: info

## Public Instance of Cheminformatics Microservice
Both the  public instance of the **Cheminformatics Microservice API server** & **frontend interface** is hosted by [Friedrich Schiller University Jena](https://www.uni-jena.de) in Germany and is accessible at:

👉 **Cheminformatics Microservice API server** - [https://api.naturalproducts.net/latest/docs](https://api.naturalproducts.net/latest/docs)

👉 **Frontend interface** [https://app.naturalproducts.net](https://app.naturalproducts.net)


:::

</div>

</div>

## Citation guidelines
<div style="text-align: justify;">

It is strongly recommended that readers cite both the software and the corresponding paper when utilizing this work. By providing appropriate citations for the Cheminformatics Microservice, readers gain a convenient means to precisely track the original source of the utilized source code and data.
</div>


- Citing software:

```
Venkata, C., Sharma, N., Schaub, J., Steinbeck, C., & Rajan, K. (2023). cheminformatics-microservice (Version v1.5.0) [Computer software].
https://doi.org/10.5281/zenodo.7745987
```

## Acknowledgments and Maintainence

<div style="text-align: justify;">

Cheminformatics Microservice and [Natural Products Online](https://naturalproducts.net) are developed and maintained by the [Steinbeck group](https://cheminf.uni-jena.de) at the [Friedrich Schiller University](https://www.uni-jena.de/en/) Jena, Germany.

Funded by the [Deutsche Forschungsgemeinschaft (DFG, German Research Foundation)](https://www.dfg.de/) under the [National Research Data Infrastructure – NFDI4Chem](https://nfdi4chem.de/) – Project number: **441958208** and  [ChemBioSys](https://www.chembiosys.de/en/) (Project INF) - Project number: **239748522 - SFB 1127**.
</div>

<p align="center">
<a href="https://www.dfg.de/" target="_blank"><img src="./public/dfg_logo_schriftzug_blau_foerderung_en.gif" width="50%" alt="DFG Logo"></a>
<a href="https://cheminf.uni-jena.de/" target="_blank"><img src="https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true" width="50%" alt="cheminf Logo"></a></p>
