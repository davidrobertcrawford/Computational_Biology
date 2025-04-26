# "fgsea_tables_plots.R"
# 2022_03_28


## There are 5 steps in this script:
#  (1) Load gene sets (AKA pathways) and gene ranks
#  (2) Run GSEA
#  (3) Remove redundant results using "collapse" function
#  (4) Build ranked results tables
#  (5) Plot ranked results tables

# fgsea on Bioconductor: https://bioconductor.org/packages/release/bioc/html/fgsea.html
# fgsea bioRxiv paper: https://doi.org/10.1101/060012


# Forematter --------------------------------------------------------------

rm(list=ls())
gc()

library(tidyverse)
library(fgsea)
library(data.table)


# Load inputs -------------------------------------------------------------

## Load gene sets

# These files were downloaded directly from MSigDB
gobp_sigs = gmtPathways("MSigDB/c5.go.bp.v7.4.symbols.gmt")
loc_sigs = gmtPathways("MSigDB/c1.all.v7.4.symbols.gmt")
full_sigs = c(gobp_sigs, loc_sigs)

# Add custom gene set
custom_sig = fread("example_gene_set.txt", header=FALSE)
full_sigs$example_gene_set = custom_sig$V1

length(full_sigs) # 7760
rm(gobp_sigs, loc_sigs, custom_sig)


## Load gene rank data
# > These example ranks are correlation coefficients
load("example_gene_ranks.RData")


# Run GSEA ----------------------------------------------------------------

## Set parameters

# Precision of p-value 
# >> A larger value like 10^5 saves time
eps_set = 0           

# Permutations 
# >> Increase if there's a warning later on about needing more permutations
perm_set = 1e5     

# Random seed
seed_start = 2022

# Minimum gene set size 
# >> I tend to ignore really small sets
min_size = 20

# Maximum gene set size
#  >> I tend to find really large sets uninformative
#  >> Larger sets also take more time to process (permutations...)
max_size = 200        

# Cutoff for adjusted p-values 
p_cutoff = 0.05

# Build named list of gene ranks
gene_names = example_gene_ranks$Gene_name
ranks = example_gene_ranks$Gene_score
names(ranks)=gene_names

# Run GSEA function
# >> scoreType is "std" because our ranks include positive and negative values
# >> Standard output will include positively and negatively enriched gene sets 
set.seed(seed_start)
gsea_temp = fgsea(full_sigs, ranks, scoreType="std",
                  eps=eps_set, nPermSimple=perm_set,
                  minSize=min_size, maxSize=max_size)

file_name = "Output/gsea_full.RData"
save(gsea_temp, file=file_name)


# Collapse results --------------------------------------------------------

# Remove results that don't pass FDR
gsea_cut = gsea_temp[padj<=p_cutoff,]

## Set parameters for collapse function
# >> This removes redundant gene sets (see bioRxiv paper pp.31-32)
# >> For two sets A and B, A is non-redundant given B only if either:
#     (1) A is enriched even if we don't consider genes it shares with B
#     (2) A is enriched even if we consider only genes shared with B

# Set threshold for redundancy tests just described
collapse_pval = 0.05

# Set permutations
# > The ratio given is what is recommended in the fgsea documentation
collapse_nperm = 10/collapse_pval

# Run collapse function
collapse_temp = collapsePathways(gsea_cut, pathways=full_sigs,
                                 stats=ranks, pval.threshold=collapse_pval,
                                 nperm=collapse_nperm)

file_name = "Output/gsea_collapsed.RData"
save(collapse_temp, file=file_name)

# Keep FDR-passed results only if they are included in the not-redundant list
collapse_gsea_temp = gsea_cut[pathway %in% collapse_temp$mainPathways,]


# Build ranked tables -----------------------------------------------------

# We are ranking by the normalized enrichment score (NES)

# Separate collapsed results into (+) and (-) for ranking
temp_pos = collapse_gsea_temp[NES>0,]
temp_neg = collapse_gsea_temp[NES<0,]

temp_pos[,NES_rank:=frankv(temp_pos$NES, order=-1, ties.method="min")]
temp_neg[,NES_rank:=frankv(temp_neg$NES, order=1, ties.method="min")]

setkey(temp_pos, NES_rank)
setkey(temp_neg, NES_rank)

file_name = "Output/gsea_collapsed_ranked.RData"
save(temp_pos, temp_neg, file=file_name)


# Plot results ------------------------------------------------------------

# Number of top hits to plot
top_number = 5

top_pos = temp_pos[NES_rank<=top_number,]
top_neg = temp_neg[NES_rank<=top_number,]

top_pos_picks = full_sigs[top_pos$pathway]
top_neg_picks = full_sigs[top_neg$pathway]

# I manually save these plots in RStudio using the "Export... PDF" option.
# There might be a way to automatically save these, but a simple "ggsave()" call doesn't work.

dev.off()
plotGseaTable(top_pos_picks, ranks, fgseaRes=top_pos, gseaParam=1,
              colwidths=c(8,5,1,1.6,1.6)) # Manually save as 18"x3" PDF

dev.off()
plotGseaTable(top_neg_picks, ranks, fgseaRes=top_neg, gseaParam=1,
              colwidths=c(8,5,1,1.6,1.6)) # Manually save as 18"x3" PDF

