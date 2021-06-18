library(Biostrings)
library(tidyverse)

source("./helper/update_data.r")
source("./helper/translate_mod.r")
update_data()

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

# read seqs
ref_seq <- readDNAStringSet("../data/reference.fasta")
# system("mafft --6merpair --thread -10 --keeplength --addfragments ../data/gisaid_hcov-19_2021_06_02_03.fasta ../data/reference.fasta > ../data/gisaid_hcov-19_2021_06_02_03.fasta_29903.fasta")

# system("pangolin ../data/gisaid_hcov-19_2021_06_02_03.fasta_29903.fasta -o ../results")
seqs <- readDNAStringSet(sample)

date_tmp="2021-06-04"
orf <- read_csv("../../2020-09-01_COVID_NGS_pipeline/scripts/ORF_SCoV2.csv")
seq_s <- subseq(seqs, orf$start[orf$sequence=="S"], orf$stop[orf$sequence=="S"])
seq_s <- chartr("-", "N", seq_s)
seq_s_aa <- translateGappedAln(seq_s)

## refine deletion and gaps with VOC info
pos_deletion <- read_csv("../../2021-01-25_COVID_NGS_analysis/data/voc_deletion.csv")

seqs_voc_nt_del <- seq_s
seqs_voc_aa_del <- seq_s_aa
sapply(seq_len(nrow(pos_deletion)), function(i){
	print(i)
	pos_i <- pos_deletion$Deletion[i]
	check1 <- subseq(seqs_voc_aa_del, pos_i, pos_i) == "-"
	check2 <- subseq(seqs_voc_aa_del, pos_i-2, pos_i+2) == "-----"
	check <- check1 & (!check2)
	if(sum(check)==0){return("")}
	seq_i <- seqs_voc_aa_del[check]

	at_i <- IRanges(pos_i, pos_i)
	seq_i <- replaceAt(seq_i, at_i, "-") 
	seqs_voc_aa_del[check] <<- seq_i

	seq_nt_i <- seqs_voc_nt_del[check]
	at_i <- IRanges((pos_i*3-2), (pos_i*3))
	seq_nt_i <- replaceAt(seq_nt_i, at_i, "---") 
	seqs_voc_nt_del[check] <<- seq_nt_i
})

writeXStringSet(seqs_voc_nt_del, paste0("../results/Thai_spike_nt_", date_tmp, ".fasta"))
writeXStringSet(seqs_voc_aa_del, paste0("../results/Thai_spike_aa_", date_tmp, ".fasta"))

# check lineage specific mutations
source("./helper/check_voc_snps.r")
df_snps_check <- check_voc_snps(seqs, mc.cores = 30)

df_snps_check$date <- sapply(df_snps_check$id, function(x){
	tmp <- strsplit(x, "|", fixed = T)[[1]]
	tmp <- tmp[length(tmp)]
	check <- strsplit(tmp, "-")[[1]]
	if(length(check) == 1){
		tmp <- paste0(tmp, "-01-01")
	} else if(length(check) == 2){
		tmp <- paste0(tmp, "-01")
	}
	return(tmp)
})
df_snps_check <- df_snps_check %>% select(id, date, everything())
df_snps_check$date <- lubridate::ymd(df_snps_check$date)
df_snps_check <- df_snps_check %>% arrange(desc(date))
file_out <- strsplit(sample, "/", fixed = T)[[1]]
file_out <- file_out[length(file_out)]
writexl::write_xlsx(df_snps_check, paste0("../results/df_check", file_out, ".xlsx"))

