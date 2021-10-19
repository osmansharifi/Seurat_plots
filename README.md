
# Single-nucleus RNA-seq 5’ reveals cell-to-cell heterogeneity in Rett Syndrome mouse brain temporally

This pipeline is developed to analyze 10X single nucleus RNA-seq 5' data from mouse brain. It includes mapping of reads via cellranger, preprocessing, DEG analysis, gene ontology analysis, trajectory analysis and velocity analysis.


## Study Design

#### Test Conditions
  
| Disease State      | In utero          | Pre-symptomatic    | Early symptomatic | Late symptomatic  |
|:------------------:|:-----------------:|:------------------:|:-----------------:|:-----------------:|
| **Mouse Age**      | E18               | P30                | P60               |      P120/P150    |
| **Brain Region**   | whole brain       |Cortex/Hypothalamus |Cortex/Hypothalamus|Cortex/Hypothalamus|
| **Sex**            | M/F               | M/F                |    M/F            |         M/F       |
| **Mecp2-e1 Genotype**| -/y  +/y  +/+  -/+|-/y  +/y  +/+  -/+  |-/y  +/y  +/+  -/+ |-/y  +/y  +/+  -/+ |     

This table shows the different conditions of samples that were collected to gain insight into Rett syndrome progression overtime. 
## Pipeline
![Pipeline](https://github.com/osmansharifi/snRNA-seq-pipeline/blob/master/figures/snRNA-seq%20Pipeline.png)

## Run Locally

Clone the project

```bash
  git clone https://github.com/osmansharifi/snRNA-seq-pipeline
```
  
## Contributors

- [@osmansharifi](https://github.com/osmansharifi)[![twitter](https://img.shields.io/badge/twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/osmansharifi3)
- [@keithfraga](https://github.com/xperthunter)[![linkedin](https://img.shields.io/badge/linkedin-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](https://linkedin.com/in/keith-fraga-56b025102)
- [@vikihaghani](https://github.com/vhaghani26)[![twitter](https://img.shields.io/badge/twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/vikihaghani26)
- [@iankorf](https://github.com/iankorf)

  
