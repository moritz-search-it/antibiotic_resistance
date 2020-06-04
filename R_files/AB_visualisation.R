library(grDevices)
library(ggplot2)
library(tidyverse)


#Database comparison-------------------------------------------------------------------------------------

##ncbi file
ncbi <- read.csv2("ncbi.csv",header = TRUE)
###ncbi enterococcus faecis
ncbi_faecis<- apply(ncbi[c(1:4),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
names_faecis <- c("59161","59162","59167","59168")
rownames(ncbi_faecis) <- names_faecis  #ncbi_faecis_logical[0,] gives the ABs, ncbi_faecis_logical[,0]  gives the names
ncbi_faecis <- as.data.frame(ncbi_faecis) 
### ncbi E.coli
ncbi_coli<- apply(ncbi[c(5:6),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
names_coli <-c("HV114-1","HV292") 
ncbi_coli <- as.data.frame(ncbi_coli,row.names = names_coli)


##card file
card <- read.csv2("card.csv",header = TRUE)
card_faecis<- apply(card[c(1:4),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
rownames(card_faecis) <- names_faecis
card_faecis <- as.data.frame(card_faecis) 
card_coli<- apply(card[c(5:6),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
card_coli <- as.data.frame(card_coli,row.names = names_coli)

##argannot file
argannot<- read.csv2("argannot.csv",header = TRUE)
argannot_faecis<- apply(argannot[c(1:4),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
rownames(argannot_faecis) <- names_faecis
argannot_faecis <- as.data.frame(argannot_faecis) 
argannot_coli<- apply(argannot[c(5:6),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
argannot_coli <- as.data.frame(argannot_coli,row.names = names_coli)

##resfinder file
resfinder <- read.csv2("resfinder.csv",header = TRUE)
resfinder_faecis<- apply(resfinder[c(1:4),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
rownames(resfinder_faecis) <- names_faecis
resfinder_faecis <- as.data.frame(resfinder_faecis) 
resfinder_coli<- apply(resfinder[c(5:6),-c(1:2)],c(1,2), function(x){if(x == "100"){x <- TRUE}else{x <- FALSE}})
rownames(resfinder_coli) <- names_coli
resfinder_coli <- as.data.frame(resfinder_coli,row.names = names_coli)

apply(resfinder_coli,1,sum) #number of positive antibiotic genes respond in  E. coli strains from resfinder database
sum(resfinder_coli[1,]== "TRUE") #same but only first  (HV114-1) E. coli strain 

# Positive AB response-----------------------------------------------------------------------------------------------
##enterococcus faecis
a <- apply(resfinder_faecis,1,sum)
b <- apply(card_faecis,1,sum)
c <- apply(argannot_faecis,1,sum)
d <- apply(ncbi_faecis,1,sum)
faecis_t <- data.frame("resfinder"=a,"card"=b,"argannot"=c,"ncbi"=d)

library(reshape2)
df <- melt(faecis_t)
df$rowid <- names_faecis
df

##Eschericia coli
a <- apply(resfinder_coli,1,sum)
b <- apply(card_coli,1,sum)
c <- apply(argannot_coli,1,sum)
d <- apply(ncbi_coli,1,sum)
coli_t <- data.frame("resfinder"=a,"card"=b,"argannot"=c,"ncbi"=d)
df2 <- melt(coli_t)
df2$rowid <- names_coli
df2

# Plot--------------------------------------------------

## Database comparisaison

my_title <- expression("Screening of different databases for ab-resistant \n genes shows unequal results")

### faecis
database <- ggplot(df,aes(rowid,value,fill=variable))+
  geom_col(position = "dodge",col="black",width = 0.5)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
    labs(x="Enteroccous faecis strains",y= " # of antibiotic resistance genes",fill="Databases", title = my_title)+
     theme(legend.position="bottom",plot.title = element_text(hjust = 0.2,vjust = -2.5))
database
ggsave("result_databases_faecis.png",width = 5, height = 5)

### e.coli

database2 <- ggplot(df2,aes(rowid,value,fill=variable))+
  geom_col(position = "dodge",col="black",width = 0.5)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(x="Escherichia coli strains",y= " # of antibiotic resistance genes",fill="Databases", title = my_title)+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.2,vjust = -2.5))
database2
ggsave("result_databases_coli.png",width = 5, height = 5)

# Grouped Antibiotics gene--------------------------

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
