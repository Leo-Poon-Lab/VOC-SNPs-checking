library(Biostrings)
library(tidyverse)
library(rjson)

# read reference 
ref_data <- fromJSON(file = "./constellations/constellations/data/SARS-CoV-2.json")
ref_data_gene <- lapply(ref_data$genes, function(x){
	tibble(name = tolower(x$name), from = as.numeric(x$coordinates$from), to = as.numeric(x$coordinates$to))
})
ref_data_gene <- bind_rows(ref_data_gene)
write_csv(ref_data_gene, "../data/ref_data_gene.csv")
ref_data_proteins <- lapply(ref_data$proteins, function(x){
	tibble(name = tolower(x$name), gene = tolower(x$gene), from = as.numeric(x$coordinates$from), to = as.numeric(x$coordinates$to))
})
ref_data_proteins <- bind_rows(ref_data_proteins)
write_csv(ref_data_proteins, "../data/ref_data_proteins.csv")

source("./helper/translate_mod.r")

# read constellations
files_cons <- list.files("./constellations/constellations/definitions/", "json", full.names = T)
sites_cons <- lapply(files_cons, function(x){
	tmp <- fromJSON(file = x)
	name <- tmp[[1]]
	sites <- tmp$sites
	return(list(name, sites))
})
df_defining_snps <- lapply(seq_along(sites_cons), function(i){
	lineage_i <- sites_cons[[i]]
	name_i <- lineage_i[[1]]
	lineage_i <- lineage_i[[2]]
	return(tibble(lineage = name_i, snps = paste(lineage_i)))
})
df_defining_snps <- bind_rows(df_defining_snps)
write_csv(df_defining_snps, "../results/df_defining_snps.csv")

# wave 3 and 4
# info from https://doi.org/10.1016/j.lanwpc.2021.100130, table 1
sites_cons_w3 <- list(list("B.1.1.63", paste0("nuc:", c("G28881A", "G28882A", "G28883C", "C2973T", "C20312T", "C21597T", "C28308G", "C29144T"))))
sites_cons_w4 <- list(list("B.1.36.27", paste0("nuc:", c("G25563T", "G922A", "T1947C", "G3431T", "T5653C", "G5950A", "C6255T", "C7504T", "C18877T", "C22444T", "T24175C", "C26060T", "C26735T", "C28854T"))))

sites_cons <- c(sites_cons_w3, sites_cons_w4, sites_cons) 

# get lineage stat
check_voc_snps <- function(seqs, mc.cores = 2){
	stopifnot(class(seqs)[1] == "DNAStringSet")
	stopifnot(width(seqs) == 29903)
	df_check <- lapply(seq_along(sites_cons), function(i){
		lineage_i <- sites_cons[[i]]
		name_i <- lineage_i[[1]]
		print(name_i)
		lineage_i <- lineage_i[[2]]
		df_check_i <- mclapply(seq_along(lineage_i), function(j){
			print(j)
			y <- lineage_i[j]
			if(grepl("+", y, fixed = T)){
				return(rep(NA, length(seqs))) # ignore insertion
			}
			gene_t <- tolower(strsplit(y, ":", fixed = T)[[1]][1])
			gene_t <- gsub("orf", "", gene_t)
			pos_t <- strsplit(y, ":", fixed = T)[[1]][2]
			pos_t_sim <- as.numeric(gsub("\\D", "", pos_t))
			sub_t <- strsplit(pos_t, "")[[1]]
			sub_t <- sub_t[length(sub_t)]
			length_t <- as.numeric(strsplit(y, ":", fixed = T)[[1]][3])
			
			if(is.na(length_t)){
				length_t <- 1
			}
			
			check <- (grepl("del", gene_t, fixed = T)) | (grepl("nuc", gene_t, fixed = T))
			if(check){ # nucleotide position
				seq_t <- as.character(subseq(seqs, pos_t_sim, pos_t_sim+length_t-1))
				if(grepl("del", gene_t, fixed = T)){  # deletion 
					rst <- grepl("-", seq_t)
				} else {
					rst <- seq_t == sub_t
				}
			} else {  # AA position
				seq_gene_1a <- subseq(seqs, 266, 13468)
				seq_gene_1b <- subseq(seqs, 13468, 21555)
				seq_gene_1ab <- xscat(seq_gene_1a, seq_gene_1b)
				if(grepl("1ab", gene_t)){
					seq_gene_t <- subseq(seq_gene_1ab, pos_t_sim*3-2, pos_t_sim*3)
				} else if(grepl("1a", gene_t)){
					seq_gene_t <- subseq(seq_gene_1a, pos_t_sim*3-2, pos_t_sim*3)
				} else if(grepl("1b", gene_t)){
					seq_gene_t <- subseq(seq_gene_1b, pos_t_sim*3-2, pos_t_sim*3)
				} else if(grepl("nsp", gene_t)){
					ref_data_proteins_t <- ref_data_proteins %>% filter(name == gene_t)
					stopifnot(nrow(ref_data_proteins_t)>0)
					seq_gene_t <- subseq(seq_gene_1ab, as.numeric(ref_data_proteins_t$from)*3-2, as.numeric(ref_data_proteins_t$to)*3)
					seq_gene_t <- subseq(seq_gene_t, pos_t_sim*3-2, pos_t_sim*3)
				} else {
					ref_data_proteins_t <- ref_data_proteins %>% filter(gene == gene_t)
					stopifnot(nrow(ref_data_proteins_t)>0)
					ref_data_gene_t <- ref_data_gene %>% filter(name == ref_data_proteins_t$name)
					seq_gene_t <- subseq(seqs, as.numeric(ref_data_gene_t$from), as.numeric(ref_data_gene_t$to))
					seq_gene_t <- subseq(seq_gene_t, pos_t_sim*3-2, pos_t_sim*3)
				}
				
				seq_gene_aa_t <- translateGappedAln(seq_gene_t)
				if(grepl("-", pos_t, fixed = T)){  # deletion 
					rst <- grepl("-", seq_gene_aa_t)
				} else {
					rst <- as.character(seq_gene_aa_t) == sub_t
				}
			}
			return(rst)
		}, mc.cores = mc.cores)
		df_check_i <- as_tibble(as.data.frame(df_check_i))
		names(df_check_i) <- lineage_i
		num_mut <- apply(df_check_i, 1, sum, na.rm = T)
		name_i <- paste0(name_i, "(total n=", length(lineage_i), ")")
		return(tibble(id = seq_along(seqs), name = name_i, value = num_mut))
	})

	df_out <- bind_rows(df_check)
	df_out <- df_out %>% pivot_wider(names_from = name, values_from = value)
	df_out$id <- names(seqs)
	return(df_out)
}

