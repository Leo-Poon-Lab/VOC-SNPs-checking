library(Biostrings)
library(tidyverse)
library(rjson)

source("./helper/update_data.r")
update_data()

df_orf <- read_csv("../data/ORF_SCoV2.csv")

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

df_defining_snps_nt <- lapply(seq_len(nrow(df_defining_snps)), function(i){
	y <- df_defining_snps$snps[i]
	gene_t <- tolower(strsplit(y, ":", fixed = T)[[1]][1])
	gene_t <- gsub("orf", "", gene_t)
	pos_t <- strsplit(y, ":", fixed = T)[[1]][2]
	pos_t_sim <- as.numeric(gsub("\\D", "", pos_t))
	sub_t <- strsplit(pos_t, "\\d")[[1]]
	ref_t <- sub_t[1]
	alt_t <- sub_t[length(sub_t)]

	check <- (grepl("del", gene_t, fixed = T)) | (grepl("nuc", gene_t, fixed = T))
	if(check){ # nucleotide position
		pos_nt_start <- pos_t_sim
		pos_nt_stop <- pos_t_sim
	} else {  # AA position
		# seq_gene_1a <- subseq(seqs, 266, 13468)
		# seq_gene_1b <- subseq(seqs, 13468, 21555)
		# seq_gene_1ab <- xscat(seq_gene_1a, seq_gene_1b)
		if(grepl("1ab", gene_t)){
			check_1a <- pos_t_sim*3-2 <= (13468-266+1)
			if(check_1a){
				pos_nt_start <- pos_t_sim*3-2 + 266 -1
				pos_nt_stop <- pos_t_sim*3 + 266 -1
			} else {
				pos_nt_start <- pos_t_sim*3-2 + 266 
				pos_nt_stop <- pos_t_sim*3 + 266 
			}
		} else if(grepl("1a", gene_t)){
			pos_nt_start <- pos_t_sim*3-2 + 266 -1
			pos_nt_stop <- pos_t_sim*3 + 266 -1
		} else if(grepl("1b", gene_t)){
			pos_nt_start <- pos_t_sim*3-2 + 266
			pos_nt_stop <- pos_t_sim*3 + 266
		} else {
			if(gene_t=="nsp12"){
				if(pos_t_sim<=9){gene_t <- "nsp12_1"}else{gene_t <- "nsp12_2"}
			} 
			if(gene_t=="spike"){
				gene_t <- "s"
			} 
			if(sum(tolower(df_orf$sequence) == gene_t)>0){
				ref_data_proteins_t <- df_orf %>% filter(tolower(sequence) == gene_t)
			} else {
				gene_t <- paste0("orf",gene_t)
				ref_data_proteins_t <- df_orf %>% filter(tolower(sequence) == gene_t)
			}
			print(gene_t)
			stopifnot(nrow(ref_data_proteins_t)>0)
			pos_nt_start <- pos_t_sim*3-2 + ref_data_proteins_t$start -1
			pos_nt_stop <- pos_t_sim*3 + ref_data_proteins_t$stop -1
		}
	}	
	tibble(variant=df_defining_snps$lineage[i], orf=gene_t, pos_orf=pos_t_sim, pos_nt_start=pos_nt_start, pos_nt_stop=pos_nt_stop, ref=ref_t, alt=alt_t)
})
df_defining_snps_nt <- bind_rows(df_defining_snps_nt)
write_csv(df_defining_snps_nt, "../results/df_defining_snps_nt.csv")
