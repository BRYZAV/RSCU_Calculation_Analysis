
#Install Packages
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

if (!require("reshape2")) install.packages("reshape2"); library(reshape2)

if (!require("dplyr")) install.packages("dplyr"); library(dplyr)

if (!require("readxl")) install.packages("readxl"); library(readxl)

if (!require("readr")) install.packages("readr"); library(readr)

if (!require("caper")) install.packages("caper"); library(caper)

if (!require("nlme")) install.packages("nlme"); library(nlme)

if (!require("purrr")) install.packages("purrr"); library(purrr)

#Inputting Tree and csv files

#joining the tips to the assembly_IDs
yeast_key<-read_excel("1154yeasts_21outrgoups_info_20220408.xlsx")
yeast_key<-yeast_key[,c("assembly_fullID_updated","Tip_ID_on_tree","clade")]

all_tree<-read.nexus("1175taxa_1403OGs_iqtree_ML_rooted.nex")
fixed_env<-read.delim("updated_env_plant_vict_animal_env_abbe_10_5.txt")
yeast_data <- read.csv("Yeast_Matrix.csv")


#Adding the TipIDs to yeast dataframe
tree.tip <- yeast_key$Tip_ID_on_tree
yeast_data['Tip_ID_on_tree'] <- c(yeast_key$Tip_ID_on_tree)
yeast_data <- as.data.frame(yeast_data)



#Empty values for the variables
summ <- c() #empty list for linear model summary
pvalue <- c() #empty list for p values
coeff <- c() #empty list for slope
lam <- c() #empty list for lambda values
codvs <- c() #empty list for codon labels (AAA vs AAC)
rSqr <- c() #empty list for r squared
adSqr <- c() #empty list for adjusted r squared

#lets see if the seed changes things
eff_seed <- sample(1:2^15, 1)
print(sprintf("Seed for session: %s", eff_seed))
set.seed(eff_seed)


#create a safe version of pgls to bypass any errors that occur
pgls.possible=possibly(.f=pgls, otherwise=NULL)

#Index of where the column of the codon starts
index<-3

#Goes through all the indices that contains the codon values
while(index<=66){
  
  column<-3
  #Goes through all the indices that contains the codon values to test against each other  
  while (column<=66){
  #Test codons that are not the same e.g(AAA vs AAA)
    if(column != index)
      #Run PGLS
      summ <- pgls.possible( yeast_data[,c(index)] ~ yeast_data[,c(column)], comparative.data(all_tree,yeast_data,"Tip_ID_on_tree"), lambda = "ML", bounds = list(lambda=c(0.001,1), kappa=c(1e-6,3), delta=c(1e-6,3)) )
      #Append the codons that are being tested against each other
      codvs <- append(codvs, print(paste0(names(yeast_data)[index]," vs ", names(yeast_data)[column])))
      
      #If an error occurs append "NULL" in the results, else continue appending the output results
      if(is.null(summ)==T){
        coeff <- append(coeff, "NULL") 
        
        pvalue <- append(pvalue, "NULL")
        
        lam <- append(lam, "NULL")

        rSqr <- append(rSqr, "NULL")

        adSqr <- append(adSqr, "NULL")
      }else{
          coeff <- append(coeff, signif(summary(summ)$coefficients[2,1], 5)) 
    
          pvalue <- append(pvalue, signif(summary(summ)$coefficients[2,4], 5))
      
          lam <- append(lam, summ$param[2])

          rSqr <- append(rSqr, signif(summary(summ)$r.squared, 5))

          adSqr <- append(adSqr, signif(summary(summ)$adj.r.squared, 5))

      }
      #Print the column that is being tested
      print("column is")
      print(column)
    }
    
    
    column<-column+1
  }
  
  index<-index+1
  print("index is")
  print(index)
}


PGLS_df <- data.frame(codons=codvs, slope=coeff, p_value=pvalue, lambda=lam, r_squared=rSqr, adusted_r_squared=adSqr)

#Output results to a csv file
write.csv(PGLS_df,"RSCU_PGLS_profile.csv", row.names=FALSE)
