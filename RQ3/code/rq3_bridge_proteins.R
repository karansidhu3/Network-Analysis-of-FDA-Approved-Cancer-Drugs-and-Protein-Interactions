# RQ3 --> Identify the protein targets that act as bridges between 
# otherwise distinct drug communities. 


library(igraph)
# Paths relative to project root
data_dir <- "RQ3/data"
out_dir  <- "RQ3/results"
# load the CandidateDrug4Cancer file
candidate <- read.csv(
  file.path(data_dir, "cancer_0.2.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# First - Build the drug-target bipartite graph 
tmp <- data.frame(
  drug = toupper(candidate$ChEMBL_ID),
  target = toupper(candidate$Target_ID),
  lab = candidate$label,
  stringsAsFactors = FALSE
)
# remove duplicate rows 
tmp <- unique(tmp)

# keep edges where label = 1 (which are FDA approved cancer drugs)
# and drug and target are not missing 
tmp <- tmp[tmp$lab == 1 & !is.na(tmp$drug) & !is.na(tmp$target), ] 

# keep only the drug/target columns for the bipartite graph 
drug_target_edges <- tmp[, c("drug", "target")]
# this result gives us clean and unique drug-target pairs 



# Next - Build the bipartite graph and project to drug-drug network 
drug_vertices <- data.frame(name = unique(drug_target_edges$drug),
  type = TRUE,
  stringsAsFactors = FALSE
)
protein_vertices <- data.frame( name = unique(drug_target_edges$target),
  type = FALSE,
  stringsAsFactors = FALSE
)
# combine vertex sets 
vertices <- rbind(drug_vertices, protein_vertices)

# build the bipartite drug-protein graph
g_bipartite <- graph_from_data_frame(d = drug_target_edges,
  directed = FALSE,
  vertices = vertices
)
proj <- bipartite_projection(g_bipartite)
g_drug <- proj$proj2  # drugs network (linked if sharing a target)
vcount(g_drug) 
ecount(g_drug) 



# Next - Find drug communities using Louvain in the drug-drug network
set.seed(219)
drug_co <- cluster_louvain(g_drug)

drug_communities <- data.frame(
  drug = V(g_drug)$name, # use vertex names from g_drug for durg IDs
  drug_co = membership(drug_co), # drug_co = drug community ID 
  stringsAsFactors = FALSE
)


# Next - for each protein target we count how many drug communities hit it 
# join the drug-target pairs with drug communities by drug
dt_with_commu <- merge(
  x = drug_target_edges,
  y = drug_communities,
  by = "drug"
)

# get the unique protein targets
targets <- unique(dt_with_commu$target)
n_targets <- length(targets)

# n_drug_communities = how many distinct drug communities contain drugs that hit this protein
# n_drugs = how many cancer drugs hit this protein overall 
n_drug_communities <- integer(n_targets) # number of drug communities per target
n_drugs <- integer(n_targets) # number of drugs per target
drug_communities_str <- character(n_targets) # a list of community IDs per target

# loop over each target to get the counts 
for(i in seq_along(targets)) {
  cur_target <- targets[i]
  
  # rows per the specific target we want 
  subset_rows <- dt_with_commu[dt_with_commu$target == cur_target, ]
  
  # all the distinct community IDs of drugs that hit the target
  comms <- sort(unique(subset_rows$drug_co))
  
  # all distinct drugs that hit the target
  drugs_for_target <- unique(subset_rows$drug)
  
  # then fill in the vectors
  n_drug_communities[i] <- length(comms)
  n_drugs[i] <- length(drugs_for_target)
  drug_communities_str[i] <- paste(comms, collapse = ",")
}

# create the data frame that describes each protein/target & connection to drug communities
protein_drug_comm <- data.frame(
  target = targets,                         # target ID CHEMBL
  n_drug_communities = n_drug_communities,  # number of distinct drug communities
  n_drugs = n_drugs,                        # number of drugs targeting protein 
  drug_communities = drug_communities_str,  # community IDs 
  stringsAsFactors= FALSE
)
head(protein_drug_comm)


# Next - extract bridge proteins 
# proteins targeted by drugs from at least 2 different communities
# i.e. bridge distinct drug communities 
bridge_proteins_of_2_or_more <- protein_drug_comm[
  protein_drug_comm$n_drug_communities >= 2,
]
# sort by (number of communities desc, then number of drugs desc)
bridge_proteins_of_2_or_more <- bridge_proteins_of_2_or_more[
  order(
    -bridge_proteins_of_2_or_more$n_drug_communities,
    -bridge_proteins_of_2_or_more$n_drugs
   ),
]



# Finally - add a column to use the gene symbols/protein names for consistency 
# read the new/clean target name file that contains the Protein/gene symbol
target_info <- read.csv(
  file.path(data_dir, "target_name_clean.csv"),
  stringsAsFactors = FALSE
)

# keep only the columns we want (target id and gene symbol)
target_map <- target_info[, c("target_id", "Protein")]
target_map$target_id <- toupper(target_map$target_id)


# merge bridge proteins using Protein names
bridge_upd <- merge(
  x = bridge_proteins_of_2_or_more,
  y = target_map,
  by.x = "target",
  by.y = "target_id",
  all.x = TRUE  #keep all the bridge targets 
)
bridge_upd

write.csv(
  bridge_upd, 
  file.path(out_dir, "bridge_targets_q3.csv"),
  row.names = FALSE
)
