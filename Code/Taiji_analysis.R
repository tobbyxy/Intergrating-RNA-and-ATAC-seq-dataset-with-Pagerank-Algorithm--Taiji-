
# Author : Oluwatobiloba Aminu
#Taiji pagerank data analysis


import("ggpubr")
#####
#read in scores and p values from taiji output
generank_e7.5 <- read.table("./Data/e7.5/GeneRanks.tsv")
generank_e7.5.pval <- read.table("./Data/e7.5/GeneRanks_PValues.tsv")
generank_e8.5 <- read.table("./Data/e8.5/GeneRanks.tsv")
generank_e8.5.pval <- read.table("./Data/e8.5/GeneRanks_PValues.tsv")
generank_e9.5 <- read.table("./Data/e9.5/GeneRanks.tsv")
generank_e9.5_pval<- read.table("./Data/e9.5/GeneRanks_PValues.tsv")


#Combine all the timepoints, filter based on 0.01
pval.com <- cbind(generank_e7.5.pval, generank_e8.5.pval, generank_e9.5_pval)
pval.filter <- pval.com %>% filter(e7.5 <= 0.01, e8.5 <= 0.01, e9.5 <= 0.01)

#get the tfs pagerank scores using the filtered p values
new.mat <- cbind(generank_e7.5, generank_e8.5, generank_e9.5)

#matrix contains duplicates
new.mat <- new.mat[, c("e7.5", "e8.5", "e9.5")]
mat.com <-  new.mat[rownames(new.mat) %in% rownames(pval.filter),]

#modify the tf names for access
tf.names <- strsplit(rownames(mat.com), "_")
tf.names2 <- data.frame(t(as.data.frame(tf.names)))

colnames(tf.names2)<- c("tf_name", "id")

#tf matrix contians all version of ids
tf.mat <- cbind(mat.com, tf.names2)
tf.mat$tf_name <- str_to_title(tf.mat$tf_name)


###write.table(tf.mat$tf_name, "./all_sig_tfs.txt", quote = F, col.names = F, row.names = F)



######## 
#Expression scores (TPM)

tpm <- cbind(e7.5abundance, e8.5abundance, e9.5abundance)
colnames(tpm) <- toupper(colnames(tpm))
mean_tpm <- data.frame(trans=rownames(tpm),
                       E7.5_mean=NA,
                       E8.5_mean=NA,
                       E9.5_mean=NA)

#Some replicates are outliers so they are removed from analysis
#get mean tpm of all time points

mean_tpm$E7.5_mean <- rowMeans(tpm[, c("E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
mean_tpm$E8.5_mean <- rowMeans(tpm[, c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
mean_tpm$E9.5_mean <- rowMeans(tpm[, c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5")])
row.names(mean_tpm) <- tpm$ENS_GENE

mean_tpm$trans <- rownames(mean_tpm)
colnames(mean_tpm)[1] <- "ens_gene"


#convert ensenmble ids to gene names
#different ensemble gene ids(ens_gene) can have the same gene name (ext_gene)
#different transcript can also come from the same gene

t2g <- read.table("./Data/t2g.txt", header = T)
t2g.tmp <- t2g %>% arrange(ens_gene)

t2g.tmp <- t2g.tmp %>% distinct(ext_gene, ens_gene)
merged.mean.tpm <- merge(mean_tpm, t2g.tmp, by="ens_gene")



#get top driver transcription factors at different time points
#use only trancription factors that show variation across time points, as top drivers
#get sd, means and covariance

#b method
tf.mat$sd <- rowSds(as.matrix(tf.mat[, -c(4,5)]))
tf.mat$means <- rowMeans2(as.matrix(tf.mat[, -c(4,5)]))
tf.mat$cov <- tf.mat$sd/tf.mat$means

#check normality for all data generated. (e7.5, e8.5, e9.5)
#1. is the data normal

shapiro.test(tf.mat$e9.5)

# p < 2.2e-16 data is not normal
# get top drivers  (n=100) cov highest

#merge pageranks and expression values
#tf.mat.all.merge contains alls tfs variable that has mean values(important)
#some tfs names are replicated because the result from pagerank had motifs names with multiple version, hence multiple tfs
tf.mat.means <- merged.mean.tpm[merged.mean.tpm$ext_gene %in% tf.mat$tf_name,]
tf.mat.all.merge <- merge(tf.mat.means, tf.mat, by.x = "ext_gene", by.y="tf_name")


#check why nrows of tf.mat.all.merge and tf.mat.means don't match
#get a unique tf or the tf with the higher pagerank score

tf.mat.uniq <- tf.mat %>% distinct(tf_name, .keep_all = T)
tf.mat.all.merge_2 <- merge(tf.mat.means, tf.mat.uniq, by.x = "ext_gene", by.y="tf_name")

# tf.mat.all.merge_2 should be used since it contains no duplicate tfs 
#307 driver transcription factors passed the filter


#######
#what is the correlation between the expression scores and pagerank scores

ggscatter(tf.mat.all.merge, x = "E7.5_mean", y = "e7.5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mean TPM", ylab = "Page Rank Score")
ggscatter(tf.mat.all.merge, x = "E8.5_mean", y = "e8.5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mean TPM", ylab = "Page Rank Score")
ggscatter(tf.mat.all.merge, x = "E9.5_mean", y = "e7.5", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Mean TPM", ylab = "Page Rank Score")

ggsave("E9.5_corrplot.pdf")

#########
#read in marker genes for comparison
#How does the marker genes change with expression levels and page rank scores
e7.5_epc_sgt <- read.table("./Data/0_e7.5_TSCspecificationMarkers.txt")


epc.marker.merge <- tf.mat.all.merge[tf.mat.all.merge$ext_gene %in% e7.5_epc_sgt$V1, c(1, 3,4,5,6,7,8)]

#take the log of both expression and page ranks
#center and scale the log of only the expression values 
log.epc.marker.merge <- log2(epc.marker.merge[,-1]+1)
zscore.epc.merge <- t(scale(t(log.epc.marker.merge[, 4:6]), scale=TRUE, center=TRUE))

#update with (+ve and _ve expression values)
log.epc.marker.merge[, 4:6] <- zscore.epc.merge
log.epc.marker.merge$tfs <- epc.marker.merge$ext_gene

#modify column names
colnames(log.epc.marker.merge)[4:6] <- c("E7.5_pagerank", "E8.5_pagerank", "E9.5_pagerank")

#convert to tidy data
tidy.epc.marker <- log.epc.marker.merge %>% 
  pivot_longer(-tfs, 
               names_to = c("timepoint", ".value"), 
               names_sep="_" )
colnames(tidy.epc.marker) <- c("tfs", "timepoint","Log Scaled TPM", "Normalised Page Rank")

# plot: dot plot
ggplot(data = tidy.epc.marker, aes(x = timepoint, y = tfs, 
                                   color = `Normalised Page Rank`, size = `Log Scaled TPM`)) + 
  geom_point() +
  scale_color_gradient(low = "white", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("Markers for E7.5 TSC Specification")
ggsave("./Figures/pagerank_e7.5_tsc_dotplot.pdf")

#results: we see that some marker genes have increased expression and pagerank scores at E7.5 (Highest)


#now we can write a function that takes in any number of markers of interest


marker_genes_file <- "./Data/0_e9.5_branching-labyrinthMarkers.txt"
output_file <- "./Figures/pagerank_e9.5_branching_labyrinth.pdf"
plot_marker_genes <- function(marker_genes_file, output_file) {
  data <- read.table(marker_genes_file)
  
  epc.marker.merge <- tf.mat.all.merge[tf.mat.all.merge$ext_gene %in% data$V1, c(1, 3,4,5,6,7,8)]
  
  log.epc.marker.merge <- log2(epc.marker.merge[,-1]+1)
  zscore.epc.merge <- t(scale(t(log.epc.marker.merge[, 4:6]), scale=TRUE, center=TRUE))
  
  log.epc.marker.merge[, 4:6] <- zscore.epc.merge
  log.epc.marker.merge$tfs <- epc.marker.merge$ext_gene
  
  colnames(log.epc.marker.merge)[4:6] <- c("E7.5_pagerank", "E8.5_pagerank", "E9.5_pagerank")
  
  tidy.epc.marker <- log.epc.marker.merge %>% 
    pivot_longer(-tfs, 
                 names_to = c("timepoint", ".value"), 
                 names_sep="_" )
  colnames(tidy.epc.marker) <- c("tfs", "timepoint","Log Scaled TPM", "Normalised Page Rank")
  
  p <- ggplot(data = tidy.epc.marker, aes(x = timepoint, y = tfs, 
                                          color = `Normalised Page Rank`, size = `Log Scaled TPM`)) + 
    geom_point() +
    scale_color_gradient(low = "white", high = "red") +
    theme_bw() + 
    ylab("") + 
    xlab("") + 
    ggtitle("Markers for E9.5 TSC Specification")
  
  ggsave(output_file, plot = p)
}
plot_marker_genes(marker_genes_file, output_file)

############
#can we identify novel driver transcription factors


top.drivers.all.cov <- tf.mat.all.merge_2 %>% arrange(desc(cov)) 

top.drivers.all.cov[top.drivers.all.cov$ext_gene == 'Ddit3',]
#remove one of the Ddit3 gene because it is a duplicate and choose one with row means value
top.drivers.all.cov2 <- top.drivers.all.cov[!(top.drivers.all.cov$ens_gene == 'ENSMUSG00000116429'
),]

top.drivers.all.cov.mat <- top.drivers.all.cov2 %>% select("E7.5_mean", "E8.5_mean", "E9.5_mean", 
                                                           "e7.5", "e8.5", "e9.5")

#take log of all values
#take zcore of only pagerank
log.top.drivers.all.cov.mat <- log2(top.drivers.all.cov.mat+1)
zscore.all.cov.mat <- t(scale(t(log.top.drivers.all.cov.mat[, 4:6]), scale=TRUE, center=TRUE))

rownames(log.top.drivers.all.cov.mat) <- top.drivers.all.cov2$ext_gene
#replace with zcores

log.top.drivers.all.cov.mat [,4:6] <- zscore.all.cov.mat
top.drivers.all.cov[duplicated(top.drivers.all.cov$ext_gene,),]
#why is Ddit3 still duplicated.


colnames(log.top.drivers.all.cov.mat)[4:6] <- c("E7.5_pagerank", "E8.5_pagerank", "E9.5_pagerank")

#i will reduce the top ranked tf's into 4 quantiles based on their coefficient of variance

qu_top_25 <- subset(top.drivers.all.cov2, top.drivers.all.cov2$cov > quantile(top.drivers.all.cov2
                                                                              $cov, 0.75))
# qu_a <-subset(top.drivers.all.cov2, 
#               top.drivers.all.cov2$cov < quantile(top.drivers.all.cov2$cov, 0.75) & 
#                 top.drivers.all.cov2$cov > quantile(top.drivers.all.cov2$cov, 0.50))
# qu_b <- subset(top.drivers.all.cov2, 
#                top.drivers.all.cov2$cov < quantile(top.drivers.all.cov2$cov, 0.50) & 
#                  top.drivers.all.cov2$cov > quantile(top.drivers.all.cov2$cov, 0.25))
# qu_bottom_25 <- subset(top.drivers.all.cov2, top.drivers.all.cov2$cov < quantile(top.drivers.all.cov2$cov, 0.25))

#pivot for first quadrant
top_25_log.top.drivers.all.cov.mat <- log.top.drivers.all.cov.mat[qu_top_25$ext_gene,]
top_25_log.top.drivers.all.cov.mat$tfs <- rownames(top_25_log.top.drivers.all.cov.mat)


#check the top drivers for interesting association between the means and pagerank scores
#filter for interesting associations.

E9.5_most_var_tfs <- top_25_log.top.drivers.all.cov.mat %>% 
  filter(E7.5_mean < E9.5_mean & E7.5_pagerank < E9.5_pagerank)

E7.5_most_var_tfs <- top_25_log.top.drivers.all.cov.mat %>% 
  filter(E7.5_mean > E9.5_mean & E7.5_pagerank > E9.5_pagerank)


tidy.all.pr.expr <- E7.5_most_var_tfs %>% 
  pivot_longer(-tfs,
               names_to = c("timepoint", ".value"), 
               names_sep="_" )
# after tidy dat
# plot dotplot

#rename col names
colnames(tidy.all.pr.expr) <- c("tfs", "timepoint","Log Scaled TPM", "Normalised Page Rank")
# plot: dot plot
pdf("./Figures/top_ranked_E7.5tfs.pdf", width=8, height=8, onefile=FALSE)
ggplot(data = tidy.all.pr.expr, aes(x = timepoint, y = tfs, 
                                    color = `Normalised Page Rank`, size = `Log Scaled TPM`)) + 
  geom_point() +
  scale_color_gradient(low = "white", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("Top ranked E9.5 TF's")
dev.off()
#ggsave("Top_ranked_tfs_dotplot.pdf")

plot_top_ranked_tfs(E7.5_most_var_tfs, "./Figures/top_ranked_E7.5tfs.pdf")
#ggsave("Top_ranked_tfs_dotplot.pdf")

#yes we identified new driver transcription factors for placenta development
#LHX6, Xbp1m Nrf1, Irf1 at E7.5
#we also confirmed already known placenta transcription factors from our analysis

#Gata, Foxd, ctcf.

#In conclusion, we would continue the analysis with other methods
