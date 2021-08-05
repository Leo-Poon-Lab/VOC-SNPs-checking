update_data <- function(x){
	if(file.exists("./constellations/data/SARS-CoV-2.json")){
		setwd("./constellations")
		system("git pull")
		setwd("../")
	} else {
		system("git clone https://github.com/cov-lineages/constellations.git")
	}
}

