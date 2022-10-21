
#First, we will look at the distribution of the taxon-specific LB and tip-to-root values for the supermatrix to see how the distribution is
LBscore_Taxon_density <- density(LB_scores_perTaxon$Matrix_original_supermatrix.fas.treefile)
plot(LBscore_Taxon_density)
TR_Taxon_density <- density(TR_scores_perTaxon$Matrix_original_supermatrix.fas.treefile)
plot(TR_Taxon_density)

#To more easily determine a threshold we will display them as histograms as well
LBscore_Taxon_hist <- hist(LB_scores_perTaxon$Matrix_original_supermatrix.fas.treefile)
plot(LBscore_Taxon_hist)
TR_Taxon_hist <- hist(TR_scores_perTaxon$Matrix_original_supermatrix.fas.treefile)
plot(TR_Taxon_hist)

#We will look now on the correlation of the different loci-measurements
# conduct a correlation tests between LB_score_Heterogeneity to all other columns
cor.test (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$LB_score_upper_quartile)
cor.test (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$Tip_to_Root_upper_quartile)
cor.test (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$Tip_to_Root_Heterogeneity)
cor.test (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$Average_PD)
# plot the same
plot (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$LB_score_upper_quartile)
plot (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$Tip_to_Root_upper_quartile)
plot (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$Tip_to_Root_Heterogeneity)
plot (LB_scores_summary_perPartition$LB_score_Heterogeneity, LB_scores_summary_perPartition$Average_PD)

#As LB and TR are measurements are highly correlated, we will concentrate now on LB_score_Heterogeneity and Average_PD for the loci
#We will generate density plots and histograms of both to determine thresholds again
LBscore_Loci_density <- density(LB_scores_summary_perPartition$LB_score_Heterogeneity)
plot(LBscore_Loci_density)
AvePD_Loci_density <- density(LB_scores_summary_perPartition$Average_PD)
plot(AvePD_Loci_density)
LBscore_Loci_hist <- hist(LB_scores_summary_perPartition$LB_score_Heterogeneity, breaks = 16)
plot(LBscore_Loci_hist)
AvePD_Loci_hist <- hist(LB_scores_summary_perPartition$Average_PD, breaks = 16)
plot(AvePD_Loci_hist)

#generate a heatmap of the values of taxon per locus for the LB scores
#delete the first column with the supermatrix values from the table
LB_scores_perTaxon_woSupermatrix <- LB_scores_perTaxon[,-(1)]
LB_scores_perTaxon_woSupermatrix <- LB_scores_perTaxon_woSupermatrix[,-(101)]

#install gplots for heatmap.2
require(gplots)
library(gplots)
#generate a heat map with hierarchical clustering and key legend
dev.off()
LB_scores_perTaxon_woSupermatrix_matrix <- heatmap.2(as.matrix(LB_scores_perTaxon_woSupermatrix), Rowv=TRUE,Colv=TRUE, col=rainbow(256,s=1,v=1,start=0.1,end=0,alpha=1),scale="none",margins=c(20,20),hclustfun=hclust,dendrogram="both",symm=FALSE,key=TRUE,trace="none",tracecol="black",density.info="density")
