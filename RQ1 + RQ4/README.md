# RQ1 and RQ4 — Network Analysis of FDA-Approved Cancer Drugs and Protein Interactions

This folder contains all data, code, and results related to **Research Question 1 (RQ1)** of the project: *Which human-body proteins serve as central hubs in the FDA-approved cancer drug-target interaction network, and what biological significance do they hold in terms of sensitivity or resistance?*, as well as **Research Question 4 (RQ4)** of the project: *Which proteins are co-targeted by multiple FDA-approved cancer drugs, and what does their network connectivity suggest about potential synergistic effects or unintended cross-reactivity?*


## Data

- **CandidateDrugs4Cancer Dataset**:
  - Original file is **not included in the repo** due to size (>1GB) but can be downloaded/accessed from [CandidateDrugs4Cancer](https://drive.google.com/file/d/1gXpGc5UhAYB9zVYnSl6E2MEaWfl5w98O/view)

- **Oncology.xlsx**: Key oncology drug dataset used for network mapping.

- **BIOGRID protein interaction dataset**:  
  - Original file (`BIOGRID-ALL-5.0.251.tab3.txt`) is **not included in the repo** due to size (>1GB).  
  - You can download it from [BioGRID](https://thebiogrid.org/downloads.php).

---

## R Code
  - [Rcode file](Rcode/projecrtRefined%20R%20console.txt)
---

## Results

- **Figures**:
  - [Image for RQ1](results/results-rq1.png): Top Hubs
  - [Image for RQ4](results/results-rq4.png): Co-targeted proteins

---

## Usage

 - Ensure you have the required R packages installed (e.g., `igraph`, `readr`, `dpylr`).

---

## Notes

- Large datasets (>100MB) are **excluded from the repository**.  
- This README provides a high-level overview to help replicate and understand the analysis workflow for RQ2.

---
