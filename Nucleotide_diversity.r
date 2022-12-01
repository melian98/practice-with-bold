#Introduction
#Frogs are a group of tailless amphibians that make up the order of “Anura”. These amphibians can be found all over the world, in both tropical regions and closer to the hemispheres of the earth. There exist over 7100 species of frogs all around the world. Upon inspecting the BOLD database, it was observed that while over half of frogs existed along the equator of the earth, there were still many frog species in the northern and southern parts of the world. This assignment aims to assess the genetic variation in the barcode gene “COI-5P” for a chosen species of frog that exists in both the warmer tropical regions of the earth as well as the colder “hemispheres” and to compare the amount of genetic diversity for the frogs in the tropical regions as opposed to those in the hemispheres. 

#Studies have shown that genetic mutations increase with temperature (Chu et al., 2018). It is therefore hypothesized that genetic diversity will be higher for the frogs that exist in the tropical region of the planet. This question is of interest as it can be further expanded for not just a single species of frogs but for the entire order and potentially for other orders to detect patterns in genetic diversity with respect to temperature.  

#initialize all needed variables 
tropics_length <- NULL
hemisphere_length <- NULL
joint_length <- NULL
Bin <- NULL
difference <- NULL
best_choice <- NULL
BOLD_AAT9369_BINS <- NULL
AAT9369_nucleotides <- NULL
tropics_string_set <- NULL
hemisphere_string_set <- NULL
number_of_known_bins <- 0

#create a string of the required libraries
libraries <- c('tidyverse', 'DescTools', 'BiocManager', 'pegas', 'rmarkdown')

#install all libraries that aren't present
install.packages(setdiff(libraries, rownames(installed.packages())))
BiocManager::install("muscle")

lapply(libraries, library, character.only = TRUE)
library(muscle)

#accessing data for anura from BOLD public database 
amphibian <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Anura&format=tsv")

#write tsv using data obtained from BOLD database and save to working directory (uncomment next line if not yet done on your machine)
write_tsv(amphibian, "Anura_BOLD_data.tsv")

#Read data from the file that was previously written and assign to variable
anura <- read_tsv("Anura_BOLD_data.tsv")

#Confirm that all the data is from the order "Anura" and no unwanted data was downloaded
unique(anura$order_name)

#create subsets of your "anura" dataframe to more easily parse data:

#subset containing information regarding the location in which the data was collected
anura_location_data <- anura[ , c(1, 2, 47:59)]
#subset containing genetic data for each of the samples (BIN, nucleotide sequence etc.)
anura_genetic_information <- anura[ , c(1, 2, 69:80)]

#Check for any data entry errors by confirming that latitude and longitude do not fall below or above the min/max possible values:

#latitude data must be between -90 and 90
#histogram is used to get a visual distribution of the samples in terms of latitudes
hist(anura_location_data$lat)
min(anura_location_data$lat, na.rm = TRUE)
max(anura_location_data$lat, na.rm = TRUE)

#longitude data must be between -180 and 180
min(anura_location_data$lon, na.rm = TRUE)
max(anura_location_data$lon, na.rm = TRUE)

#Using the histogram, it is clear that there are more samples of frogs near the equator as opposed to farther away. The national geographic considers the tropics to be the regions of the earth between -23.5 degrees and 23.5 degrees. This was deemed to be too narrow a range but was used as a base point. Using figure 1, it is clear that the "anura" order exist more along the equatorial line than in the hemispheres. The tropics were chosen as the latitudes existing between -30 and 30 degrees. 

#Use summation to see exactly how many samples exist in both the tropics and the hemispheres. 
hemispheres <- sum(anura_location_data$lat > 30 | anura_location_data$lat < -30, na.rm = TRUE)
tropics <- sum(anura_location_data$lat <= 30 | anura_location_data$lat >= -30, na.rm = TRUE)

#define a dataframe which will contain the latitudes and the nucleotide sequence of all tropic anura samples with marker "COI-5P"
anura_tropics_gene_info <- subset(anura, (lat < 30 & lat > -30 & marker_codes == 'COI-5P|COI-5P' ) , select = c(processid, sampleid, recordID, bin_uri, species_name, lat, nucleotides, marker_codes))

#Confirm that all existing samples fall within the required "tropic" latitude range and only the marker_code COI-5P was obtained
#This will check which samples contain a latitude value outside of the "tropics" region. If it returns "integer(0)", then all values are in range
which(anura_tropics_gene_info$lat > 30 | anura_tropics_gene_info$lat < -30)
#This will print all unique marker codes to make sure only "COI-5P" was obtained
unique(anura_tropics_gene_info$marker_codes)

#define a dataframe which will contain the latitudes and the nucleotide sequence of all hemisphere anura samples
anura_hemisphere_gene_info <- subset(anura, ((lat >= 30 | lat <= -30) & marker_codes == 'COI-5P|COI-5P' ) , select = c(processid, sampleid, recordID, bin_uri,species_name, lat, nucleotides, marker_codes))

#Confirm that all existing samples fall within the required "hemisphere" latitude range and only the marker_code COI-5P was obtained
#This will check which samples contain a latitude value outside of the "hemisphere" region. If it returns "integer(0)", then all values are in range
which(anura_hemisphere_gene_info$lat < 30 & anura_hemisphere_gene_info$lat > 30)
#This will print all unique marker codes to make sure only "COI-5P" was obtained
unique(anura_hemisphere_gene_info$marker_codes)

#create a list that contains all BIN's that exist in both tropics and hemisphere regions
Bin <- intersect(anura_hemisphere_gene_info$bin_uri,anura_tropics_gene_info$bin_uri)

#create a dataframe that will contain the number of records for each of the BIN's present in "list". In addition, the dataframe will calculate the optimal "BIN" to use for the nucleotide diversity analysis by assessing which BIN has the most number of samples in both tropics and hemisphere regions. This will be done by adding the total number of samples in both regions together and subtracting by the difference in number of samples between the regions to calculate "Best choice"
for (i in 1:length(Bin))
{
  #This vector will determine how many samples of this particular BIN exist in the "tropics" region
   tropics_length <- c(append(tropics_length, length(grep(Bin[i], anura_tropics_gene_info$bin_uri))))
  #This vector will determine how many samples of this particular BIN exist in the "hemisphere" region
   hemisphere_length <- c(append(hemisphere_length, length(grep(Bin[i], anura_hemisphere_gene_info$bin_uri))))
  #This vector will determine how many samples of this particular BIN exist in both regions combined
   joint_length <- c(append(joint_length, length(grep(Bin[i], anura_hemisphere_gene_info$bin_uri))+ length(grep(Bin[i], anura_tropics_gene_info$bin_uri)))) 
   #number_of_known_bins will sum up all the samples with a BIN that is not n.a (i.e a known BIN). This will be used in the barplot which will form     figure 2
   if (!is.na(Bin[i])) {number_of_known_bins <- number_of_known_bins + joint_length[i]}
  #This vector will determine how great the difference in number of samples for this particular BIN is between the regions
   difference <- c(abs(append(difference, length(grep(Bin[i], anura_hemisphere_gene_info$bin_uri))- length(grep(Bin[i], anura_tropics_gene_info$bin_uri))))) 
  #The best choice for the genetic diversity analysis will be the BIN with a high "joint_length" (i.e many overall samples) and low "difference" (i    .e both regions have a relatively close number of samples) 
   best_choice <- c(append(best_choice, (joint_length[i] - difference[i])))
}

#This data frame is composed of all the results from the previous for loop and will determine which "BIN" will be used for further analysis
joint <- data.frame(Bin, tropics_length, hemisphere_length, difference,  joint_length, best_choice)

#create a barplot to show that the number of unkown BINs (where the value of the BIN is N.A) is much greater than the sum of all known BINs.
barplot_data <- c(joint$joint_length[17], number_of_known_bins)

barplot(barplot_data, names.arg = c("Unknown_Bins", "Known_BINS"), col = c("red", "green"), ylab = "Number of samples", main = "Number of known vs unknown BINS")

#using the "joint" dataframe, it is clear that the best choice is "BOLD: AAT9369" as it has the largest value in the "best_choice" column

#remove all unneeded variables and dataframes from global environment to clear memory 
#rm (i, list, tropics, hemisphere, Bin, tropics_length, hemispheres, hemisphere_length, joint_length, difference, best_choice)

#check the list location of "BOLD: AAT9369" in the tropics gene info dataframe
list <- which(anura_tropics_gene_info$bin_uri == "BOLD:AAT9369", arr.ind = TRUE)

#prepare a string which will contain the nucleotide sequence of the COI-5P gene from "BOLD: AAT9369" in the tropics region by looping through each instance of these in the tropics gene info dataframe
for (i in 1: length(list))
{
  BOLD_AAT9369_BINS <- c(append(BOLD_AAT9369_BINS, anura_tropics_gene_info$bin_uri[list[i]]))
  AAT9369_nucleotides <- c(append(AAT9369_nucleotides, anura_tropics_gene_info$nucleotides[list[i]]))
}

#prepare a dataframe which contains the sequences of the obtained COI-5P genes to be aligned in the tropics region
tropics_AAT9369 <- data.frame(BOLD_AAT9369_BINS, AAT9369_nucleotides)

#check the list location of "BOLD: AAT9369" in the hemisphere gene info dataframe
list <- which(anura_hemisphere_gene_info$bin_uri == "BOLD:AAT9369", arr.ind = TRUE)

#reset the values of the variables used in the previous for loop so they can be reused to prepare the hemisphere dataframe
BOLD_AAT9369_BINS <- NULL
AAT9369_nucleotides <- NULL

#prepare a string which will contain the nucleotide sequence of the COI-5P gene from "BOLD: AAT9369" in the hemisphere region by looping through each instance of these in the hemisphere gene info dataframe
for (i in 1: length(list))
{
  BOLD_AAT9369_BINS <- c(append(BOLD_AAT9369_BINS, anura_hemisphere_gene_info$bin_uri[list[i]]))
  AAT9369_nucleotides <- c(append(AAT9369_nucleotides, anura_hemisphere_gene_info$nucleotides[list[i]]))
}

#prepare a dataframe which contains the sequences of the obtained COI-5P genes to be aligned in the hemisphere region
hemisphere_AAT9369 <- data.frame(BOLD_AAT9369_BINS, AAT9369_nucleotides)

#create an object of class "DNAStringSet" which will contain the nucleotide of the COI_5P gene from the tropics region to be aligned
for (i in 1: length(tropics_AAT9369$AAT9369_nucleotides))
{
  tropics_string_set <- c(append(tropics_string_set, DNAStringSet(x =  tropics_AAT9369$AAT9369_nucleotides[i])))
}

#create an object of class "DNAStringSet" which will contain the nucleotide of the COI_5P gene from the hemisphere region to be aligned
for (i in 1: length(hemisphere_AAT9369$AAT9369_nucleotides))
{
  hemisphere_string_set <- c(append(hemisphere_string_set, DNAStringSet(x =  hemisphere_AAT9369$AAT9369_nucleotides[i])))
}

#use the muscle library to create a "DNAMultipleAlignment" object which contains the aligned sequences 
aligned_tropics <- muscle::muscle(tropics_string_set)
aligned_hemisphere <- muscle::muscle(hemisphere_string_set) 

#convert the aligned sequences to classes of type "DNAbin" which can be used to find the nucleotide divergence

#first convert to object of class "matrix"
matrix_tropics <- as.matrix(aligned_tropics)
matrix_hemisphere <- as.matrix(aligned_hemisphere)

#then convert to object of class "alignment"
alignment_tropics <- as.alignment(matrix_tropics)
alignment_hemisphere <- as.alignment(matrix_hemisphere)

#finally, convert to object of class "DNAbin"
DNAbin_tropics <- as.DNAbin(alignment_tropics)
DNAbin_hemisphere <- as.DNAbin(alignment_hemisphere)

#Use the "nuc.div" class from the pegas library to calculate the nucleotide diversity for both the tropics and the hemisphere
tropics_divergence <- nuc.div(DNAbin_tropics, variance = FALSE)
hemisphere_divergence <- nuc.div(DNAbin_hemisphere, variance = FALSE)

#create a second barplot. This one will compare the nucleotide diversity values from both the tropics and the hemisphere.
barplot_data <- c(tropics_divergence, hemisphere_divergence)

barplot(barplot_data, names.arg = c("tropics", "hemisphere "), col = c("blue", "orange"), ylab = "nucleotide diversity", main = "nucleotide diversity in tropics and hemisphere regions")

#Results and discussion

#By analyzing the nucleotide diversity in both the tropics and the hemisphere regions, it is noted that in both cases, diversity is low with values of 0.00202 and 0.00239 respectively as is shown in figure 3. This is also contrary to the hypothesis which claimed that diversity would be greater in the tropics due to the greater average temperature in that region. There are several potential explanations for the unexpected results. A primary explanation would be the overall low number of samples tested. While the "best_choice" vector was used to pick the species of frog with the greatest number of overall samples as well as the smallest difference between the regions, there still existed only 50 samples in the hemisphere region and 21 in the tropics (Phillips et al., 2019). concluded that while the number of samples required for genetic diversity prediction is taxon specific, an estimated 2000 individual samples would be needed to determine the genetic diversity with an interval of 95% confidence.

#The original dataframe with anura contained well over 10 000 samples but filtering for a BIN that was present in both geographic regions (hemisphere and tropics) as well as excluding data in which the BIN was not registered (N.A accounted for the highest number of samples by far in the "joint" dataframe even when all known samples were summed as shown in figure 2), the number of samples was greatly reduced. An expansion to this study may be done with a taxon that contains a larger dataframe that is more spread out geographically overall than the Anura. In addition, an expansion to this study could alter the code to include not just the BIN that was deemed to be the best choice but all BINs that are good choices in order to see the difference between individual species in a group with regards to nucleotide diversity. 

#9. References

#Phillips, J. D., Gillis, D. J., & Hanner, R. H. (2019). Incomplete estimates of genetic diversity within species: Implications for DNA barcoding. In Ecology and Evolution (Vol. 9, Issue 5, pp. 2996–3010). John Wiley and Sons Ltd. https://doi.org/10.1002/ece3.4757

#Rutledge, K., Costa, H., Boudreau, D., Hunt, J., Sprout, E., Ramroop, T., Hall, H., Teng, S., &amp; McDaniel, M. (n.d.). Tropics. National Geographic Society. Retrieved October 3, 2022, from https://education.nationalgeographic.org/resource/tropics  

#Chu, X. L., Zhang, B. W., Zhang, Q. G., Zhu, B. R., Lin, K., & Zhang, D. Y. (2018). Temperature responses of mutation rate and mutational spectrum in an Escherichia coli strain and the correlation with metabolic rate. BMC Evolutionary Biology, 18(1). https://doi.org/10.1186/s12862-018-1252-8

