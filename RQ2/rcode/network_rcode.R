library(dplyr)
library(readxl)
library(readr)
library(stringr)
library(DBI)
library(RSQLite)
library(igraph)

# Path to ChEMBL SQLite

con <- dbConnect(SQLite(), "C:/Users/karan/Downloads/chembl_36_sqlite/chembl_36/chembl_36_sqlite/chembl_36.db")

setwd("C:/Users/Karan/OneDrive/documents/Uni/cosc421/project")

# -----------------------------

# STEP 1 — Load FDA & Oncology

# -----------------------------

products <- read_tsv("data/datasatfda/Products.txt", show_col_types = FALSE) %>%
  mutate(
    DrugName = str_trim(toupper(DrugName)),
    ActiveIngredient = str_trim(toupper(ActiveIngredient))
  )

onc <- read_excel("data/oncology.xlsx", sheet = "Sheet1")
desc_col <- colnames(onc)[1]

onc_drugs <- onc %>%
  mutate(
    raw = .data[[desc_col]],
    DrugName = str_extract(
      raw,
      "(?i)(?<=approves |grants )[A-Za-z0-9\\-\\s;]+"
    ),
    DrugName = toupper(str_trim(DrugName))
  ) %>%
  filter(!is.na(DrugName)) %>%
  distinct(DrugName)



oncology_products <- products %>%
  filter(
    DrugName %in% onc_drugs$DrugName |
      ActiveIngredient %in% onc_drugs$DrugName
  )

if (!dir.exists("results")) dir.create("results")
write.csv(oncology_products, "results/oncology_drugs.csv", row.names = FALSE)

oncology_drugs <- oncology_products %>%
  mutate(
    CanonName = ifelse(ActiveIngredient != "", ActiveIngredient, DrugName),
    CanonName = toupper(str_trim(CanonName))
  ) %>%
  distinct(CanonName)

# -----------------------------------------------

# STEP 2 — ChEMBL: Drug ↦ Target ↦ Protein (UniProt)

# -----------------------------------------------

# 1. Drug → Target (tid)

drug_mech <- dbGetQuery(con, "
SELECT molregno, tid, site_id
FROM drug_mechanism
")

# 2. Molecule names

mol_dict <- dbGetQuery(con, "
SELECT molregno, pref_name
FROM molecule_dictionary
")

# 3. Target metadata (optional, not required for UniProt)

target_dict <- dbGetQuery(con, "
SELECT tid, pref_name AS target_name
FROM target_dictionary
")

# 4. Site → component mapping

site_comp <- dbGetQuery(con, "
SELECT site_id, component_id
FROM site_components
")

# 5. Component → UniProt

comp_seq <- dbGetQuery(con, "
SELECT component_id, accession, organism
FROM component_sequences
")

comp_seq_human <- comp_seq %>% filter(organism == 'Homo sapiens')

# -----------------------------

# JOIN ALL DATA CORRECTLY

# -----------------------------

drug_target <- drug_mech %>%
  left_join(mol_dict, by = "molregno") %>%
  left_join(target_dict, by = "tid") %>%
  left_join(site_comp, by = "site_id") %>%
  left_join(comp_seq_human, by = "component_id") %>%
  select(
    DrugName = pref_name,
    Target = target_name,
    UniProt = accession
  ) %>%
  filter(!is.na(UniProt)) %>%
  distinct()

# Save final file

write.csv(drug_target, "results/drug_target_network.csv", row.names = FALSE)

# Normalize drug names for matching
drug_target_clean <- drug_target %>%
  mutate(
    DrugName = toupper(str_trim(DrugName))
  )

oncology_drugs_clean <- oncology_drugs %>%
  mutate(
    CanonName = toupper(str_trim(CanonName))
  )

# ============================================================
# NEW STEP — Map FDA brand names ➝ ChEMBL generic names
# ============================================================

# Load synonyms from ChEMBL (brand names, aliases, INNs, etc.)
syn <- dbGetQuery(con, "
SELECT molregno, syn_type, synonyms
FROM molecule_synonyms
")

syn_clean <- syn %>%
  mutate(
    synonyms = toupper(str_trim(synonyms))
  )

# Match FDA drug names to any synonym in ChEMBL
matched_synonyms <- syn_clean %>%
  filter(synonyms %in% oncology_drugs_clean$CanonName) %>%
  left_join(mol_dict, by = "molregno") %>%
  select(
    CanonName = synonyms,      # FDA brand name
    ChEMBLName = pref_name     # generic name used in drug_target
  ) %>%
  mutate(
    CanonName = toupper(str_trim(CanonName)),
    ChEMBLName = toupper(str_trim(ChEMBLName))
  ) %>%
  distinct()

# Merge synonyms into the FDA oncology drug list
oncology_drugs_clean <- oncology_drugs_clean %>%
  left_join(matched_synonyms, by = "CanonName") %>%
  mutate(
    FinalName = ifelse(!is.na(ChEMBLName), ChEMBLName, CanonName)
  )

drug_target_oncology <- drug_target_clean

write.csv(drug_target_oncology, "results/drug_target_oncology.csv", row.names = FALSE)

# STEP 4 — Load and clean BioGRID PPI
biogrid <- read_tsv("data/BIOGRID-ALL-LATEST.tab3/BIOGRID-ALL-5.0.251.tab3.txt")

# Use backticks because column names contain spaces
biogrid_human <- biogrid %>%
  filter(
    `Organism Name Interactor A` == "Homo sapiens",
    `Organism Name Interactor B` == "Homo sapiens",
    `Experimental System Type` == "physical"
  ) %>%
  select(
    A = `Official Symbol Interactor A`,
    B = `Official Symbol Interactor B`
  ) %>%
  distinct()

#STEP 5
bipartite_edges <- drug_target_oncology %>%
  select(DrugName, UniProt)

bipartite_graph <- graph_from_data_frame(
  bipartite_edges,
  directed = FALSE
)

ppi_edges <- biogrid_human

ppi_graph <- graph_from_data_frame(
  ppi_edges,
  directed = FALSE
)

combined_graph <- igraph::union(bipartite_graph, ppi_graph)

deg <- degree(combined_graph)
bet <- betweenness(combined_graph)
eig <- eigen_centrality(combined_graph)$vector
clust <- transitivity(combined_graph, type = "local")

comm <- cluster_louvain(combined_graph)

save(
  bipartite_graph,
  ppi_graph,
  combined_graph,
  drug_target_oncology,
  file = "results/final_networks.RData"
)

