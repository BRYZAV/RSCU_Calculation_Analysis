if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

if (!require("reshape2")) install.packages("reshape2"); library(reshape2)

if (!require("dplyr")) install.packages("dplyr"); library(dplyr)

if (!require("readxl")) install.packages("readxl"); library(readxl)

if (!require("readr")) install.packages("readr"); library(readr)

if (!require("caper")) install.packages("caper"); library(caper)

if (!require("nlme")) install.packages("nlme"); library(nlme)

if (!require("purrr")) install.packages("purrr"); library(purrr)

all_tree<-read.nexus("1175taxa_1403OGs_iqtree_ML_rooted.nex")
fixed_env<-read.delim("updated_env_plant_vict_animal_env_abbe_10_5.txt")
yeast_key<-read_excel("1154yeasts_21outrgoups_info_20220408.xlsx")
yeast_data <- read.csv("Yeast_Matrix.csv")
yeast_features <- read.delim("y1000_features.final.12_9_22.tsv")

yeast_data['Tip_ID_on_tree'] <- c(yeast_key$Tip_ID_on_tree)
yeast_data <- as.data.frame(yeast_data)

yeast_features <- yeast_features %>% relocate(clade)
yeast_features <- yeast_features %>% relocate(Tip_ID)
yeast_features <- yeast_features[,-c(12,13)]


#Empty values for the variables
summ <- c() #empty list for linear model summary
pvalue <- c() #empty list for p values
coeff <- c() #empty list for slope
lam <- c() #empty list for lambda values
codvs <- c() #empty list for codon labels (AAA vs AAC)
rSqr <- c() #empty list for r squared
adSqr <- c() #empty list for adjusted r squared

#Run possibly to bypass any errors
pgls.possible=possibly(.f=pgls, otherwise=NULL)

#STarting on the Features file
#When doing it forreal change to 4
index <- 4

while(index <= 23) {
  
  #Running agaonst RSCU Files
  #When doing it foreal change to 3
  column <- 3
  
  
  while(column <= 66) {
    #it is Yeast_data(Y-axis) vs Yeast_features(X-axis)
    
    if (index == 6 || index == 10 || index == 11){
      
      #Extract column with NA values
      na_yeastFeatures <- c(yeast_features[index])
      #Turn it into a dataframe
      na_yeastFeatures <- data.frame(na_yeastFeatures)
      #Get species/assembly ID to identify the species with/without NA values
      na_yeastFeatures["assembly_fullID_updated"] <- c(yeast_features$assembly_fullID_updated)
      #Reoder the columns
      na_yeastFeatures = na_yeastFeatures[,c(2,1)]
      #Remove NA values
      na_yeastFeatures <- na.omit(na_yeastFeatures)
      
      #Remove NA values
      yeast_Features_mergeDF <- yeast_features[!is.na(yeast_features[index]),]
      #Reorder by clade (alphabetically)
      yeast_Features_mergeDF <- yeast_Features_mergeDF[order(yeast_Features_mergeDF$clade),]
      
      #Merge the dataframe of the yeast features with omitted na values with the yeast RSCU values dataframe
      #In thic case I want to try it out with yeast_Features_mergeDF
      #yeast_Data_mergeDF <- merge(yeast_Features_mergeDF,yeast_data,by="assembly_fullID_updated")
      yeast_Data_mergeDF <- merge(na_yeastFeatures,yeast_data,by="assembly_fullID_updated")
      
      #Reorder by clade (alphabetically)
      yeast_Data_mergeDF <- yeast_Data_mergeDF[order(yeast_Data_mergeDF$clade),]
      
      #Increase values in column by 1, NITROGEN BREADTH ONLY
      #yeast_Features_mergeDF[index] <- yeast_Features_mergeDF[index] + 1
      
      
      #Do PGLS (I do column+1 as it adds a column to the dataframe when merging so I am sure to get the right index)
      summ <- pgls.possible( yeast_Data_mergeDF[,c(column+1)] ~ yeast_Features_mergeDF[,c(index)], comparative.data(all_tree,yeast_Data_mergeDF,"Tip_ID_on_tree"), lambda = "ML" )
      codvs <- append(codvs, print(paste0(names(yeast_Data_mergeDF)[column+1]," Vs ", names(yeast_Features_mergeDF)[index])))

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
      
      column <- column + 1
      
    }
    
    summ <- pgls.possible( yeast_data[,c(column)] ~ yeast_features[,c(index)], comparative.data(all_tree,yeast_data,"Tip_ID_on_tree"), lambda = "ML" )
    codvs <- append(codvs, print(paste0(names(yeast_data)[column]," Vs ", names(yeast_features)[index])))
    
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
    
    column <- column + 1
    
  }
  
  index <- index + 1 
}

Features_df <- data.frame(features_vs_codon=codvs, slope=coeff, p_value=pvalue, lambda=lam, r_squared=rSqr, adusted_r_squared=adSqr)
write.csv(Features_df,"RSCU_FEATURES_profile.csv", row.names=FALSE)

