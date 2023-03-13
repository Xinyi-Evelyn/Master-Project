library("pheatmap")
library("dplyr")
library("dendextend")

Hgg_cluster <- readRDS("E:/Rlastest/project/the percentage of t cells in each sample cluster/all_HggSample_clusters.rds")
Lgg_cluster <- readRDS("E:/Rlastest/project/the percentage of t cells in each sample cluster/all_LggSample_clusters.rds")

Hgg_cluster <- transform(Hgg_cluster, glioma = "HGG")
Lgg_cluster <- transform(Lgg_cluster, glioma = "LGG")
all_clusters<- full_join(Hgg_cluster,Lgg_cluster)
ma <- as.matrix(all_clusters[,3:6])
colnames(ma) <- c("CTL", "DNT", "Thelper","Treg")
rownames(ma) <- as.vector(all_clusters$id)
ma <- apply(ma,2,as.numeric)
rownames(ma) <- as.vector(all_clusters$id)

annotation_row = data.frame(
  GliomaClass = factor(rep(c("HGG", "LGG"), c(106, 175)))
)
rownames(annotation_row) = rownames(ma)
ann_colors = list("HGG" = "#7570B3", "LGG" = "#E7298A")


#

pheatmap(ma,cellheight=10,cellwidth = 40,annotation_row = annotation_row, annotation_colors = ann_colors,main = "T cells in the each sample cluster",fontsize = 7)
hclust_sample <- hclust(dist(ma), method = "complete")
Tree <- as.dendrogram(hclust_sample) 
plot(Tree,horiz = TRUE)
clusters <- cutree(Tree, k = 4)
c1_data = ma[clusters == 2,]

dend <- ma %>% dist %>% hclust %>% as.dendrogram %>%  set("labels_to_character") %>% color_branches(k=4)
dend_list <- get_subdendrograms(dend, 4)

# Plotting the result
par(mfrow = c(2,3))
plot(dend, main = "Original dendrogram")
sapply(dend_list, plot)
sub_dend <- dend_list[[1]]
nleaves(sub_dend)
length(order.dendrogram(sub_dend))
subset_ma <- as.matrix(ma[order.dendrogram(sub_dend),])
order.dendrogram(sub_dend) <- rank(order.dendrogram(sub_dend))
heatmap(subset_ma, Rowv = sub_dend)

ctl_type <-data.frame(samples = rownames(subset_ma))
ctl_type$S <- grep("^_","",ctl_type$samples)


sub_dend <- dend_list[[2]]
nleaves(sub_dend)
length(order.dendrogram(sub_dend))
subset_ma <- as.matrix(ma[order.dendrogram(sub_dend),])
order.dendrogram(sub_dend) <- rank(order.dendrogram(sub_dend))
heatmap(subset_ma, Rowv = sub_dend)

sub_dend <- dend_list[[3]]
nleaves(sub_dend)
length(order.dendrogram(sub_dend))
subset_ma <- as.matrix(ma[order.dendrogram(sub_dend),])
order.dendrogram(sub_dend) <- rank(order.dendrogram(sub_dend))
heatmap(subset_ma, Rowv = sub_dend)

sub_dend <- dend_list[[4]]
nleaves(sub_dend)
length(order.dendrogram(sub_dend))
subset_ma <- as.matrix(ma[order.dendrogram(sub_dend),])
order.dendrogram(sub_dend) <- rank(order.dendrogram(sub_dend))
heatmap(subset_ma, Rowv = sub_dend)


