update_data <- function(x){
	if(dir.exists("./constellations")){
		setwd("./constellations")
		system("git pull")
		setwd("../")
	} else {
		system("git clone https://github.com/cov-lineages/constellations.git")
	}
}

