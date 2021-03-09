## Compare delta M/Z of drug metabolites across groups defined in molecular networking job

library(tidyverse)
library(magrittr)

# files
## drug annotations
an <- "data/Team2_approvedlist_hmr.csv"
## tagged drugs
tag_temp <- "data/GNPS Tag Template BETA For class project - GNPS Tag - Batch Tag Template.tsv"
## result_specnets_DB
specnets <- "data/306eebe6c5bc481d912727a08d7e80a5.tsv"
## networkedges_selfloop
edges <- "data/d582caf1956b419cb27f80cdbbf7db66..selfloop"
## spectral counts per group
counts <- "data/METABOLOMICS-SNETS-V2-160284a2-view_all_clusters_withID_beta-main.tsv"

# steps
## read in files
an <- read_csv(an, guess_max = 7000) %>%
  set_colnames(make.names(colnames(.))) # fix column names to remove special characters
tag_temp <- read_delim(tag_temp, delim = "\t") %>%
  set_colnames(make.names(colnames(.))) %>%
  filter(grepl("drug", Current_Tags, ignore.case = TRUE)) # only keep drugs
specnets <- read_delim(specnets, delim = "\t") %>%
  set_colnames(make.names(colnames(.))) %>%
  mutate(InChIKey = gsub("N/A", NA, InChIKey),
         InChIKey.Planar = gsub("N/A", NA, InChIKey.Planar)) # make "N/A" -> NA
edges <- read_delim(edges, delim = "\t") %>%
  set_colnames(make.names(colnames(.)))
counts <- read_delim(counts, delim = "\t") %>%
  set_colnames(make.names(colnames(.)))

## match cluster IDs to drug annotations
table(unique(specnets$InChIKey) %in% unique(an$InChIKey)) # only 2 matches
table(unique(specnets$SpectrumID) %in% unique(an$SpectrumID)) # only 2 matches
table(unique(specnets$InChIKey.Planar) %in% unique(an$InChIKey.Planar)) # only 3 matches

# get unique compounds (n=146)
an_unique <- an %>%
  group_by(InChIKey) %>%
  summarise(InChIKey.Planar = unique(InChIKey.Planar),
            SpectrumID = unique(SpectrumID),
            Compound_Name = unique(Compound_Name),
            name = unique(name))
# add annotations to specnets to identify clusters that are drugs
df <- inner_join(specnets, an_unique, by = "InChIKey.Planar") %>%
  select(InChIKey.Planar, X.Scan.) %>%
  bind_rows(tag_temp %>% # add tagged drugs (n=136)
              select(INPUT.InChIKey.Planar, INPUT..Scan.) %>%
              set_colnames(c("InChIKey.Planar", "X.Scan.")))

## match edge table to drug annotations by cluster index
### determine which cluster is the drug - remove edges where both clusters are drugs and delta M/Z == 0
edges_w_an <- edges %>%
  mutate(Cluster_drug = case_when(all(c(CLUSTERID1, CLUSTERID2) %in% df$X.Scan.) ~ 0, # if both clusters are drugs
                                  (CLUSTERID1 %in% df$X.Scan.) & !(CLUSTERID2 %in% df$X.Scan.) ~ 1, # if only cluster1 is a drug
                                  (CLUSTERID2 %in% df$X.Scan.) & !(CLUSTERID1 %in% df$X.Scan.) ~ 2, # if only cluster2 is a drug
                                  !all(c(CLUSTERID1, CLUSTERID2) %in% df$X.Scan.) ~ 0), # if neither clusters are drugs
         Cluster_metabolite = case_when(Cluster_drug == 1 ~ CLUSTERID2, # pick the non-drug column
                                        Cluster_drug == 2 ~ CLUSTERID1)) %>%
  filter(Cluster_drug != 0,
         DeltaMZ != 0)

## add groups to edge table
edges_w_an <- left_join(edges_w_an, counts, by = c("Cluster_metabolite" = "cluster.index"))

write_csv(edges_w_an, "drug_metabolites_delta_MZ_by_group.csv")
