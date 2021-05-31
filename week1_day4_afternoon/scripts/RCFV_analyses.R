
#Delete first column in taxon specific values
summarized_taxon_basefrequencies <- Matrix_original_supermatrix_mod.fas_summarized_taxon_basefrequencies[,-1]

#Plot density and histogram of taxonomic-specific RCFV values
taxon_RCFV_density <- density(summarized_taxon_basefrequencies$RCFV.Value..taxon.specific.)
plot(taxon_RCFV_density)
taxon_RCFV_hist <- hist(summarized_taxon_basefrequencies$RCFV.Value..taxon.specific.)
plot(taxon_RCFV_hist)

#Plot density and histogram of locus-specific RCFV values
locus_RCFV_density <- density(summarized_frequencies$RCFV.Value)
plot(locus_RCFV_density)
locus_RCFV_hist <- hist(summarized_frequencies$RCFV.Value)
plot(locus_RCFV_hist)

#generate a heatmap of the values of taxon per locus
#delete the first column with the supermatrix values from the table
RCFV_frequencies_all_partitions <- RCFV_frequencies_all_partitions[,-(1)]
#install gplots for heatmap.2
install.packages("gplots")
require(gplots)
library(gplots)
#generate a heat map with hierarchical clustering and key legend
RCFV_frequencies_all_partitions_Matrix <- heatmap.2(as.matrix(RCFV_frequencies_all_partitions), Rowv=TRUE,Colv=TRUE, col=rainbow(256,s=1,v=1,start=0.1,end=0,alpha=1),scale="none",margins=c(20,20),hclustfun=hclust,dendrogram="both",symm=FALSE,key=TRUE,trace="none",tracecol="black",density.info="density")
