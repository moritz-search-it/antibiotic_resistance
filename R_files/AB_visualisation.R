library(grDevices)
library(ggplot2)
library(tidyverse)
library(reshape2)

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
ggsave(path= "/Users/moritzherrmann/Projects/Antibiotic_resistance/","result_databases_faecis.png",width = 5, height = 5)

### e.coli

database2 <- ggplot(df2,aes(rowid,value,fill=variable))+
  geom_col(position = "dodge",col="black",width = 0.5)+
  scale_fill_brewer(palette="Blues")+
  theme_minimal()+
  labs(x="Escherichia coli strains",y= " # of antibiotic resistance genes",fill="Databases", title = my_title)+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.2,vjust = -2.5))
database2
ggsave(path= "/Users/moritzherrmann/Projects/Antibiotic_resistance/","result_databases_coli.png",width = 5, height = 5)

# Grouped Antibiotics gene--------------------------

#we use card, The CARD is a rigorously curated collection of characterized, peer-reviewed resistance determinants 
#and associated antibiotics, organized by the Antibiotic Resistance Ontology (ARO) and AMR gene detection models.
card_faecis$strain <- names_faecis
card_coli$strain <- names_coli

faecis_tidy <- melt(card_faecis) %>% 
pivot_longer(c(-"strain"),names_to = "resistance_gene", values_to= "presence")
coli_tidy <- melt(card_coli) %>% 
  pivot_longer(c(-"strain"),names_to = "resistance_gene", values_to= "presence")


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

faecis_tidy <- annotation(faecis_tidy)
coli_tidy <- annotation(coli_tidy)
print(coli_tidy,n=20)
print(faecis_tidy,n=20)

faecis_tidy %>% filter(presence==TRUE) %>% filter(strain=="59161") %>% 

colnames(tibble1[-1])
faecis_tidy %>% filter(presence==TRUE) %>%  filter(strain=="59161") %>% summarise(resistance=Group) %>%
  flatten %>% rapply(as.vector) %>% str_count("Fluoroquinolone") %>% sum() #change str_count()



tibble1 <- tibble( strain = c("59161","59162","59167","59168"),
        Aminocoumarin= c(0,0,0,0),
        Aminoglycosides = c(1,3,1,1),
        beta_lactam = c(0,0,0,0),
        broad_spectrum = c(0,2,1,1),
        Diaminopyrimidine = c(0,1,0,0),
        Fluoroquinolone=c(1,1,1,1),
        Fosfomycin=c(0,0,0,0),
        Glycopeptide=c(7,7,7,7),
        Nitroimidazole=c(0,0,0,0),
        Polypeptide= c(0,0,0,0),
        Protein_synthese_inhibitor = c(1,2,2,2),
        Sulfonamide= c(0,0,0,0),
        Tetracycline=c(1,1,2,2))

#tibble1 %>% melt %>% 
 # ggplot(aes(variable,value, group =strain)) +
  #geom_col(aes(fill=strain),position="dodge") + coord_flip() + labs( x = "resistance", y= "nr of genes")


tibble1 %>% melt %>% 
  ggplot(aes(strain,variable,fill=value)) +
  geom_tile() + xlab(label= "strains") +
  ylab(label= "Antibiotic resistance") +
  scale_fill_gradient(name= "nr. resistance genes",
                      low = "#FFFFFF", high = "#012345") +
  ggtitle(label = "Screening for antibiotic resistance genes",subtitle = "Vancomycin resistant enteroccoci (VRE)") +
  theme_bw()

ggsave("endscreening_ab_faecis.png", width=6,height = 5)

#could be good to do a function to go through all antibiotic groups
coli_tidy %>% filter(presence==TRUE) %>%  filter(strain=="HV114-1") %>% summarise(resistance=Group) %>%
  flatten %>% rapply(as.vector) %>% str_count("Fluoroquinolone") %>% sum() #change str_count()
coli_tidy %>% filter(presence==TRUE) %>%  filter(strain=="HV292") %>% summarise(resistance=Group) %>%
  flatten %>% rapply(as.vector) %>% str_count("Fluoroquinolone") %>% sum()

tibble2 <- tibble( strain = c("HV114-1","HV292-1"),
                   Aminocoumarin= c(7,7),
                   Aminoglycosides = c(8,8),
                   beta_lactam = c(19,20),
                   broad_spectrum = c(10,10),
                   Diaminopyrimidine = c(1,1),
                   Fluoroquinolone=c(19,19),
                   Fosfomycin=c(1,1),
                   Glycopeptide=c(0,0),
                   Nitroimidazole=c(1,1),
                   Polypeptide= c(3,4),
                   Protein_synthese_inhibitor = c(12,12),
                   Sulfonamide= c(1,1),
                   Tetracycline=c(11,11)
                   )

tibble2 %>% melt %>% 
  ggplot(aes(strain,variable,fill=value)) +
  geom_tile() + xlab(label= "strains") +
  ylab(label= "Antibiotic resistance") +
  scale_fill_gradient(name= "nr. resistance genes",
                      low = "#FFFFFF", high = "#012345") +
  ggtitle(label = "Screening for antibiotic resistance genes",subtitle = "Extended-spectrum-betalactamase") +
  theme_bw()
ggsave("endscreening_ab_coli.png", width=6,height = 5)



tibble2 %>% melt %>% 
  ggplot(aes(strain,variable,size=value))+
  geom_count() + scale_size_area() +scale_size(range = c(0,10))
