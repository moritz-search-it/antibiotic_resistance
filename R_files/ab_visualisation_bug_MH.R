# Grouped Antibiotics gene--------------------------
library(tidyverse)
library(reshape2)
load("visualisation.RData")

#we use card, The CARD is a rigorously curated collection of characterized, peer-reviewed resistance determinants 
#and associated antibiotics, organized by the Antibiotic Resistance Ontology (ARO) and AMR gene detection models.
card_faecis$strain <- names_faecis
card_coli$strain <- names_coli

faecis_tidy <- melt(card_faecis) %>% 
  pivot_longer(c(-"strain"),names_to = "resistance_gene", values_to= "presence")
coli_tidy <- melt(card_coli) %>% 
  pivot_longer(c(-"strain"),names_to = "resistance_gene", values_to= "presence")

faecis_tidy$group <- NA
coli_tidy$group <- NA


##add antibiotic groups
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


annotation <- function(faecis_tidy,...){
  for (i in 1:length(faecis_tidy$resistance_gene)){
    if(faecis_tidy$resistance_gene[i] %in% Aminoglycosides ){
      faecis_tidy$group[i] <- "Aminoglycosides" 
    } else if (faecis_tidy$resistance_gene[i] %in% Glycopeptide){
      faecis_tidy$group[i] <- "Glycopeptide"
    } else if (faecis_tidy$resistance_gene[i] %in% Tetracycline){
      faecis_tidy$group[i] <- "Tetracycline"
    } else if (faecis_tidy$resistance_gene[i] %in% Aminocoumarin){
      faecis_tidy$group[i] <- "Aminocoumarin"
    } else if (faecis_tidy$resistance_gene[i] == "mdtG"){
      faecis_tidy$group[i] <-"Fosfomycin"
    } else if (faecis_tidy$resistance_gene[i] == "msbA"){
      faecis_tidy$group[i] <- "Nitroimidazole"
    } else if (faecis_tidy$resistance_gene[i] == "sul2"){
      faecis_tidy$group[i] <- "Sulfonamide"
    } else if (faecis_tidy$resistance_gene[i] %in% Diaminopyrimidine ){
      faecis_tidy$group[i] <- "Diaminopyrimidine"  
    } else if (faecis_tidy$resistance_gene[i] %in% Polypeptide_ab){
      faecis_tidy$group[i] <- "Polypeptide_ab"
    } else if (faecis_tidy$resistance_gene[i] %in% Fluoroquinolone){
      faecis_tidy$group[i] <- "Fluoroquinolone"
    } else if (faecis_tidy$resistance_gene[i] %in% beta_lactam){
      faecis_tidy$group[i] <- "Beta_lactame"
    } else if (faecis_tidy$resistance_gene[i] %in% broad_spectrum){
      faecis_tidy$group[i] <- "Broad_spectrum_ab"
    } else if (faecis_tidy$resistance_gene[i] %in% Protein_synthese_inhibitor){
      faecis_tidy$group[i] <- "Protein_synthese_inhibitor_ab"
    } 
  }
  return(faecis_tidy)
}
## problem, when gene is active against multiple groups of antibiotic? Bcs if-else loop jumps to next i-entry when test expression is = TRUE
##Possible solution: (example aminoglycosides)
#n <- "Aminoglycosides"
#faecis_tidy$group[i] <- "paste(faecis_tidy$group[i],n") but problem of if-else loop still appears


faecis_tidy <- annotation(faecis_tidy)
coli_tidy <- annotation(coli_tidy)
coli_tidy
faecis_tidy
