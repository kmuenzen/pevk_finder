####### Find Orthologous Exons ######

# Read CSV files into R
all_blast_files = list.files(pattern=".csv")

# Check for empty files, store in list "empty"
info = file.info(all_blast_files)
empty = rownames(info[info$size==0,])

# Create an empty list to hold only files containing information
blast_files <- vector(length=0)

# Loop through files and create a new list with only non-empty CSV files
for (i in 1:length(all_blast_files)) {
	if (is.element(all_blast_files[i], empty) == FALSE) {
		blast_files <- c(blast_files, all_blast_files[i])
	}
}

# Create list of species names
species_names <- c("Phascolarctos_cinereus", "Homo_sapiens", "Ictidomys_tridecemlineatus", "Loxodonta_africana", "Manis_javanica", "Mus_musculus", "Mus_pahari", "Myotis_brandtii", "Ochotona_princeps", "Orcinus_orca", "Oryctolagus_cuniculus", "Pteropus_vampyrus", "Rattus_norvegicus", "Rhinolophus_sinicus", "Sorex_araneus", "Sus_scrofa", "Trichechus_manatus", "Tupaia_chinensis", "Tursiops_truncatus", "Vicugna_pacos", "Acinonyx_jubatus", "Bos_taurus","Ceratotherium_simum", "Condylura_cristata", "Dasypus_novemcinctus", "Eptesicus_fuscus", "Equus_caballus", "Erinaceus_europaeus", "Felis_catus", "Galeopterus_variegatus", "Hipposideros_armiger") 

# Loop through species names and find all corresponing csv files
for (i in 1:length(species_names)) {
	
	######### Initialize Lists and Data Frames #########

	# Initialize data frame that will hold exon data
	match_data <- data.frame(matrix(nrow=1, ncol=3))
	colnames(match_data) <- c("exon_start","exon_end","hits")

	# Initialize vector that will hold exons already entered into data frame
	matched_exons <- vector(length=0)
	
	# Initialize counter for rows in match data
	match_data_row = 0
	
	######### Get species data ##########
	species <- species_names[i]
	
	# Define sequence identifier
	species_identifier <- paste(species,"_and", sep="")
	
	# Define same sequence identifier
	same_identifier <- paste(species,"_and_",species, sep="")
	
	# Initialize species csv files
	species_files <- vector(length=0)
	
	# Loop through csv files and read in data
	for (j in 1:length(blast_files)) {
		if (grepl(species_identifier, blast_files[j]) == TRUE) {
			if (grepl(same_identifier, blast_files[j]) == FALSE) {
				species_files <- c(species_files, blast_files[j])
			}
		}
	}
	
	# Read species csv file data into R
	for (k in 1:length(species_files)) {
		assign(species_files[k], read.csv(species_files[k]))
	}

	# Loop through each csv file and compile data on exon starts/stops, excluding the self-blast
	for (m in 1:length(species_files)) {
		# Read csv file and add column names
		blast.data <- read.csv(species_files[m], header=F)
		colnames(blast.data) <- c("query_exon","subject_exon","%_identity","alignment_length","mismatch","gap_openings","query_start","query_end","subject_start","subject_end","evalue","bit_score")
		
		# Loop through current blast data
		for (n in 1:nrow(blast.data)) {
			
			##### Get info on exon starts/ends
			exon_name <- as.character(blast.data[n,1])
			split <- strsplit(exon_name, "_")
			location <- split[[1]][3]
			location <- strsplit(location, ":")
			exon_start <- as.numeric(location[[1]][1])
			exon_end <- as.numeric(location[[1]][2])
			
			# If the current exon does not exist in the list of already matched exons
			if (is.element(exon_start, matched_exons) == FALSE) {
				
				# Increment match_data row by one
				match_data_row = match_data_row + 1
				
				# Add starts/ends to match_data
				match_data[match_data_row,1] <- exon_start
				match_data[match_data_row,2] <- exon_end
				
				# Add one to number of hits
				match_data$hits[match_data_row] = 1
				
				# Add the current exon start to matched_exons
				matched_exons <- c(matched_exons, exon_start)
			}
			
			# If the current exon DOES exost in the list of already matched exons
			if (is.element(exon_start, matched_exons) == TRUE) {
				
				# Find the corresponding row in match_data
				same_row <- 0
				for (p in 1:nrow(match_data)) {
					if (exon_start == match_data[p,1]) {
						same_row <- p
					}
				}
				
				# Incremement the hit count by 1
				match_data[same_row,3] <- match_data[same_row,3] + 1
			}
		}	
	}
	
	# Write data frame to a csv file
	write.csv(match_data, file = paste(species,"_exon_hits",sep=''))
}


species_names <- c("Phascolarctos_cinereus", "Homo_sapiens", "Ictidomys_tridecemlineatus", "Loxodonta_africana", "Manis_javanica", "Mus_musculus", "Mus_pahari", "Myotis_brandtii", "Ochotona_princeps", "Orcinus_orca", "Oryctolagus_cuniculus", "Pteropus_vampyrus", "Rattus_norvegicus", "Rhinolophus_sinicus", "Sorex_araneus", "Sus_scrofa", "Trichechus_manatus", "Tupaia_chinensis", "Tursiops_truncatus", "Vicugna_pacos", "Acinonyx_jubatus", "Bos_taurus","Ceratotherium_simum", "Condylura_cristata", "Dasypus_novemcinctus", "Eptesicus_fuscus", "Equus_caballus", "Erinaceus_europaeus", "Felis_catus", "Galeopterus_variegatus", "Hipposideros_armiger") 


