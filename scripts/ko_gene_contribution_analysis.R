#dir <- 'C:/Users/ctataru/Documents/Lab/data'
#file <- 'ko_metagenome_contributions2.tab'
#filename <- paste(dir, '/', file, sep='')
####################
####################
####################
##############DONT READ###################
collapseTables <- function(data){
  library(data.table)
  DT <- data.table(data)
  final <- unique(DT[,list( sum(CountContributedByOTU)), by=.(Gene, Sample)])
  numeric_treatment <- ifelse(mapping$Treatment[match(final$Sample, mapping$SampleID)] == "Control", 0, 1)
  aut_control_multiplier <- ifelse(mapping$Treatment[match(final$Sample, mapping$SampleID)] == "Control", -1, 1)
  final_plus_treatment <- cbind(final, mapping$Treatment[match(final$Sample, mapping$SampleID)], numeric_treatment, rep(1, nrow(final)),aut_control_multiplier )
  gene_abundances <- final_plus_treatment[,list(sum(V4)), by = Gene]
  num_autism_samples = sum(mapping$Treatment == "Aut")
  num_control_samples = sum(mapping$Treatment == "Control")
  #Get those rows with genes that are present in more than 5 samples
  common_genes <- final_plus_treatment[final_plus_treatment$Gene %in% gene_abundances$Gene[gene_abundances$V1 > 5] , ]
  return(common_genes)
}



analyzeTreatmentAbundances <- function(common_genes){
  common_genes_treatments <- common_genes[, list(num_autism = sum(numeric_treatment), total = sum(V4),
                                                 prop_samples_aut =sum(numeric_treatment) / sum(V4), 
                                                 prop_samples_control = 1 - (sum(numeric_treatment) / sum(V4)),
                                                 abundance_by_cond= sum(V1*aut_control_multiplier) / sum(V1)), 
                                          by=Gene]
  
  #Abundance by condition: positive = autism, negative = control. Close to 1 means very dichotomous, close to 0.5 means evenly split
  #Of all the samples this gene is present in, how many are autism samples
  sum(common_genes_treatments$prop_samples_aut > 0.75)
  common_genes_treatments[ common_genes_treatments$prop_samples_aut > 0.75,]
  common_genes_treatments[ common_genes_treatments$prop_samples_control > 0.75,]
  highly_abundant_autism <- common_genes_treatments[common_genes_treatments$abundance_by_cond > 0.75,]
  highly_abundant_autism[highly_abundant_autism$total > 10]
}

#common_genes <- collapseTables(ko_metagenome_contributions2)
#analyzeTreatmentAbundances(common_genes)


#######YOU WANT THIS PART####################################################################################################################

getCountsPerSamplePerKo <- function(data){
  DT <- data.table(data)
  final2 <- unique(DT[,list( sum(CountContributedByOTU)), by=.(Sample, Gene)])
  final2_plus_treatment <- cbind(final2,  treatment = mapping$Treatment[match(final2$Sample, mapping$SampleID)])
  
  counts_per_sample_per_ko <- lapply(unique(final2$Gene), function(ko) return(data.frame(Samples =final2_plus_treatment$Sample[final2_plus_treatment$Gene == ko],
                                                                                         Counts =final2_plus_treatment$V1[final2_plus_treatment$Gene == ko],
                                                                                         Treatment =final2_plus_treatment$treatment[final2_plus_treatment$Gene == ko])))
  
  names(counts_per_sample_per_ko) <- unique(final2_plus_treatment$Gene)
  
 
  return(counts_per_sample_per_ko)
}


filtering <- function(data){
  #Filter out those KOs that appear in less than 6 samples
  counts_per_sample_per_ko <- getCountsPerSamplePerKo(data)
  common_kos <- sapply(counts_per_sample_per_ko, function(dataframe) return(nrow(dataframe) > 6))
  common_count_dist_per_ko <- counts_per_sample_per_ko[common_kos]
  
  #Filter out those KOs that have SD less than the mean (or 3 times the mean)
  standard_devs <- sapply(common_count_dist_per_ko, getSD)
  highly_varying_ko_dataframes <- common_count_dist_per_ko[standard_devs > mean(standard_devs)]
  very_high_varying_ko_dataframes <- common_count_dist_per_ko[standard_devs > 3 * mean(standard_devs)]
  
  return(highly_varying_ko_dataframes)
}


getSD <- function(dataframe){
  colnames(dataframe) = c("Samples", "Counts", "Treatment")
  aut_dist <- dataframe$Counts[dataframe$Treatment == "Aut"]
  control_dist <- dataframe$Counts[dataframe$Treatment == "Control"]
  return (sd(dataframe$Counts))
}

getAutVControlDist <- function(dataframe){
  colnames(dataframe) <- c("Samples", "Counts", "Treatment")
  aut_dist <- dataframe$Counts[dataframe$Treatment=="Aut"]
  control_dist <- dataframe$Counts[dataframe$Treatment == "Control"]
  return(list(aut_dist, control_dist))
}

checkNormalcyKS <- function(dataframe){
  aut_control_dist <- getAutVControlDist(dataframe)
  aut_dist <- unlist(aut_control_dist[1])
  control_dist <- unlist(aut_control_dist[2])
  #Now you just have to check if each is normal
  aut_pvalue <- ks.test(aut_dist,"pnorm",mean=mean(aut_dist),sd=sd(control_dist))$p.value
  control_pvalue <- ks.test(control_dist,"pnorm",mean=mean(control_dist),sd=sd(control_dist))$p.value
  return(c(aut_pvalue, control_pvalue))
}

checkNormalcyShapiro <- function(dataframe){
  aut_control_dist <- getAutVControlDist(dataframe)
  #Now you just have to check if each is normal
  #TODO: check if both groups have greater than 3 samples
  return(c(shapiro.test(unlist(aut_control_dist[1]))$p.value, shapiro.test(unlist(aut_control_dist[2]))$p.value))
}

wilcoxonTest <- function(dataframe){
  aut_control_dist <- getAutVControlDist(dataframe)
  return(wilcox.test(unlist(aut_control_dist[1]), unlist(aut_control_dist[2]))$p.value)
}
tTest <- function(dataframe){
  aut_control_dist <- getAutVControlDist(dataframe)
  return(t.test(unlist(aut_control_dist[1]),unlist(aut_control_dist[2]))$p.value)
}
ksTest <- function(dataframe){
  aut_control_dist <- getAutVControlDist(dataframe)
  return(ks.test(unlist(aut_control_dist[1]),unlist(aut_control_dist[2]))$p.value)
}

#Run Ks Test
runKsTest <- function(ko_dataframes){
  pvalues_ksTest <- sapply(ko_dataframes, ksTest)
  pvalues_ksTest_corrected <- p.adjust(pvalues_ksTest, method="fdr")
  imp_ind_ksTest <-  match(pvalues_ksTest_corrected[pvalues_ksTest_corrected <threshold], pvalues_ksTest_corrected)
  names(ko_dataframes[imp_ind_ksTest])
}


#Run t Test
runtTest <- function(ko_dataframes){
  pvalues_tTest <- sapply(ko_dataframes, tTest)
  pvalues_tTest_corrected <- p.adjust(pvalues_tTest,method="fdr" )
  plot(density(pvalues_tTest_corrected))
  plot(density(pvalues_tTest))
  sum(pvalues_tTest_corrected <threshold)
  imp_ind_tTest <- match(pvalues_tTest_corrected[pvalues_tTest_corrected <threshold], pvalues_tTest_corrected)
  names(ko_dataframes[imp_ind_tTest])
  
}


#Run Wilcoxon test
runWilcoxonTest <- function(ko_dataframes){
  pvalues_wilcoxon <- sapply(ko_dataframes, wilcoxonTest)
  plot(density(pvalues_wilcoxon))
  sum(pvalues_wilcoxon < threshold)
  imp_ind_wilcoxon <- match(pvalues_wilcoxon[pvalues_wilcoxon <threshold], pvalues_wilcoxon)
  names(ko_dataframes[imp_ind_wilcoxon])
}

############Main#############
###Main object is ko_dataframes, which is a list where each element represents one KO. Each element is a data.frame with column 1
### the sample id, column2 the counts of that KO in that sample, and column3 the treatment(aut/control) of that sample
data <- ko_metagenome_contributions2
ko_dataframes <- filter(data)


#Check normality
aut_control_normal_shapiro <- lapply(ko_dataframes, checkNormalcyShapiro)
plot(density(unlist(aut_control_normal_shapiro)))

aut_control_normal_ks <- lapply(ko_dataframes, checkNormalcyKS)
pvalues_ks <- unlist(aut_control_normal_ks)
plot(density(pvalues_ks))

threshold <- 0.1
sum(pvalues_ks > threshold) / length(pvalues_ks) #__% have P values over threshold

#Run Tests:
imp_kos_tTest <- runtTest(ko_dataframes)
imp_kos_wilcoxon <- runWilcoxonTest(ko_dataframes)



