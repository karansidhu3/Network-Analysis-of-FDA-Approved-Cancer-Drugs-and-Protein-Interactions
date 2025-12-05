# RQ2 — Network Analysis of FDA-Approved Cancer Drugs and Protein Interactions

This folder contains all data, code, and results related to **Research Question 2 (RQ2)** of the project: *Are there clusters or communities of drugs that target similar regions of the protein interaction network, suggesting common mechanisms of action or potential cross-reactivity between treatments?*

## Data

- **FDA datasets** (`datasatfda/`):
  - `Applications.txt`, `Products.txt`, etc.  
  - Contains information on FDA-approved drugs and submissions.

- **Oncology.xlsx**: Key oncology drug dataset used for network mapping.

- **BIOGRID protein interaction dataset**:  
  - Original file (`BIOGRID-ALL-5.0.251.tab3.txt`) is **not included in the repo** due to size (>1GB).  
  - You can download it from [BioGRID](https://thebiogrid.org/downloads.php).

---

## R Code

- **graph_creation.R**: Scripts to construct and visualize drug-protein networks.
- **network_rcode.R**: Functions for network metrics, clustering, and community detection.
- **project_rq.R**: Main analysis script for RQ2.

---

## Results

- **CSV outputs**:
  - `drug_communities_projection.csv`
  - `drug_target_network.csv`
  - `drug_target_oncology.csv`
  - `oncology_drugs.csv`

- **Figures**:
  - `ppi_community_snapshot_white.png`: Community structure of the protein interaction network.
  - `plot_drug_projection_top80.png`: Visualization of top drug communities.

- **RData files**:
  - `final_networks.RData`: Saved network objects for reproducibility.

---

## Usage

1. Ensure you have the required R packages installed (e.g., `igraph`, `tidyverse`).
2. Download the **BIOGRID dataset** and place it in the `RQ2/data/BIOGRID-ALL-LATEST.tab3/` folder if needed.
3. Run the scripts in the following order:
   1. `graph_creation.R`
   2. `network_rcode.R`
   3. `project_rq.R`
4. Results will be automatically generated in the `results/` folder.

---

## Notes

- Large datasets (>100MB) are **excluded from the repository**.  
- This README provides a high-level overview to help replicate and understand the analysis workflow for RQ2.

---
