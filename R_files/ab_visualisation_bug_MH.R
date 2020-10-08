# Grouped Antibiotics gene--------------------------
library(tidyverse)
library(reshape2)
load("visualisation.RData")

#we use card, The CARD is a rigorously curated collection of characterized, peer-reviewed resistance determinants 
#and associated antibiotics, organized by the Antibiotic Resistance Ontology (ARO) and AMR gene detection models.
card_faecis$strain <- names_faecis

head(card_faecis)
card_coli$strain <- names_coli

faecis_tidy <- melt(card_faecis) %>% 
  pivot_longer(c(-"strain"),names_to = "resistance_gene", values_to= "presence")
coli_tidy <- melt(card_coli) %>% 
  pivot_longer(c(-"strain"),names_to = "resistance_gene", values_to= "presence")

#add antibiotic groups
Aminoglycosides <- c("aad.6.","aadA5","acrD","AAC.6...Ii","APH.3....Ib","APH.3...IIIa","APH.6..Id","baeR","baeS","cpxA","kdpE","tolC")
Glycopeptide <-  c("vanA","vanHA","vanRA","vanSA","vanXA","vanYA","vanZA")
Tetracycline <- c("acrB","acrS","emrK","emrY","Escherichia_coli_acrA","evgA","evgS","H.NS","marA","msrC","tet.B.","tetM","tolC")
Aminocoumarin <-c("baeR","baeS","cpxA","mdtA","mdtB","mdtC","tolC")
Polypeptide_ab <- c("bacA","eptA","pmrF","ugd","yojI")
Fluoroquinolone <- c("acrB","acrE","acrF","acrS","CRP","efmA","emrA","emrB","emrR","Escherichia_coli_acrA","evgA","evgS","gadW","gadX","H.NS","marA","mdtE","mdtF","mdtH","tolC")
beta_lactam <- c("acrB","acrE","acrF","acrS","CRP","CTX.M.1","Escherichia_coli_acrA","Escherichia_coli_ampC","Escherichia_coli_ampC1_beta.lactamase","Escherichia_coli_ampH","evgA","evgS","gadW","gadX","H.NS","marA","mdtE","mdtF","TEM.1","tolC")
broad_spectrum <- c("SAT.4","acrB","acrS","Escherichia_coli_acrA","Escherichia_coli_mdfA","marA","mdtM","mdtN","mdtO","mdtP","msrC","tolC")
Protein_synthese_inhibitor <- c("CRP","efmA","ErmB","Escherichia_coli_emrE","evgA","evgS","gadW","gadX","H.NS","mdtE","mdtF","mdtM","mphB","msrC","tolC")
Diaminopyrimidine <- c("dfrG","dfrA17")

annotation<-function(tidy_data){

  #create an empty vector of lists where we will store each abgene's annotation
  groups <- vector("list",length=nrow(tidy_data))


  for (i in 1:nrow(tidy_data)){
    
    #empty vector to store each abgene's annotation
    anno <- c()
    
    
    #creating the annotation for each abgene
    if(tidy_data$resistance_gene[i] %in% Aminoglycosides ){
      anno <- c(anno,"Aminoglycosides") 
    } 
    if (tidy_data$resistance_gene[i] %in% Glycopeptide){
      anno <- c(anno,"Glycopeptide")
    } 
    if (tidy_data$resistance_gene[i] %in% Tetracycline){
      anno <- c(anno,"Tetracycline")
    } 
    if (tidy_data$resistance_gene[i] %in% Aminocoumarin){
      anno <- c(anno,"Aminocoumarin")
    } 
    if (tidy_data$resistance_gene[i] == "mdtG"){
      anno <-c(anno,"Fosfomycin")
    } 
    if (tidy_data$resistance_gene[i] == "msbA"){
      anno <- c(anno,"Nitroimidazole")
    } 
    if (tidy_data$resistance_gene[i] == "sul2"){
      anno <- c(anno,"Sulfonamide")
    } 
    if (tidy_data$resistance_gene[i] %in% Diaminopyrimidine ){
      anno <- c(anno,"Diaminopyrimidine")  
    } 
    if (tidy_data$resistance_gene[i] %in% Polypeptide_ab){
      anno <- c(anno,"Polypeptide_ab")
    } 
    if (tidy_data$resistance_gene[i] %in% Fluoroquinolone){
      anno <- c(anno,"Fluoroquinolone")
    }
    if (tidy_data$resistance_gene[i] %in% beta_lactam){
      anno <- c(anno,"Beta_lactame")
    }
    if (tidy_data$resistance_gene[i] %in% broad_spectrum){
      anno <- c(anno,"Broad_spectrum_ab")
    } 
    if (tidy_data$resistance_gene[i] %in% Protein_synthese_inhibitor){
      anno <- c(anno,"Protein_synthese_inhibitor_ab")
    } 
    
    #storing the annotation 
    groups[[i]] <- anno
  
   }

  #return the annotated dataset
  return(tidy_data %>% mutate(Group =groups))
}
## problem, when gene is active against multiple groups of antibiotic? Bcs if-else loop jumps to next i-entry when test expression is = TRUE
##Possible solution: (example aminoglycosides)
#n <- "Aminoglycosides"
#faecis_tidy$group[i] <- "paste(faecis_tidy$group[i],n") but problem of if-else loop still appears


faecis_tidy <- annotation(faecis_tidy)
warcoli_tidy <- annotation(coli_tidy)


#note: if you wamt to access say the groups of the 5th data entry as just a character vector use: faecis_tidy$Group[[5]] 
#the double brackets are important because I have stored the annotations as elements of a list data structure
#therefore if you just use single brakets R considers that element as a list. 

#a useful function to check what data type an object is:
typeof(faecis_tidy$Group[5])
typeof(faecis_tidy$Group[[5]])

