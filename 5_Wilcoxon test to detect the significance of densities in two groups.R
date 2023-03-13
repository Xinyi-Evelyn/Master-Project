setwd("E:/Rlastest/project/density_heatmap")
library("ggplot2")
library(pheatmap)


############################Step1 
df <- read.csv("pan_density.csv")
df2 <- df[1:24,c(-1,-2,-10)]
matrix <- as.matrix(df2)
matrix
rownames(matrix) <- c("_H_P001105","_H_P001205","_H_P001705","_H_P002802","_H_P001202","_H_P001401","_H_P002607-A1",
                      "_H_P002607-D1","_H_RA-024","_H_RA-030",'_H_318766-1','L_290866-1','L_290876-1','L_290878-1','L_290880-1','L_290884-1',
                      'L_290887-1','L_318768-1' ,'L_318774-1','L_318776-1','L_318778-1','L_318780-1','L_318782-1','L_318784-1')
idx <- which(matrix == 0)
matrix[idx] <- NA
matrix
data <- log(matrix)


##############################Step2 
df3 <- read.csv("cell_density.csv")
df4 <- df3[1:24,c(-1,-2,-11)]
matrix2 <- as.matrix(df4)
rownames(matrix2) <- df[1:24,1]
matrix2
rownames(matrix2) <- c("_H_P001105","_H_P001205","_H_P001705","_H_P002802","_H_P001202","_H_P001401","_H_P002607-A1",
                       "_H_P002607-D1","_H_RA-024","_H_RA-030",'_H_318766-1','L_290866-1','L_290876-1','L_290878-1','L_290880-1','L_290884-1',
                       'L_290887-1','L_318768-1' ,'L_318774-1','L_318776-1','L_318778-1','L_318780-1','L_318782-1','L_318784-1')
idx <- which(matrix2 == 0)
matrix2[idx] <- NA
matrix2
data2 <- log(matrix2)
data3 <- cbind(data,data2)
dim(data3)

#
rownames(data3)
rownames(data3) <- 1:24
annotation_row = data.frame(
  GliomaClass = factor(rep(c("HGG", "LGG"), c(11, 13)))
)
rownames(annotation_row) = rownames(data3)
ann_colors = list("HGG" = "#7570B3", "LGG" = "#E7298A")
pheatmap(data3,annotation_row = annotation_row, annotation_colors = ann_colors, main = "The log(cells density) from both T panel and Pan panel",cellwidth = 30, cellheight = 15)



######Wilcoxon test
data3 <- as.data.frame(data3)
data3$grade <- "HGG"
data3[12:24,16] <- "LGG"

wilcox.test(data3$T.cell~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = T.cell))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("T cells") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$Mac~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = Mac))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("Macrophages") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())

wilcox.test(data3$Mac.PDL1pos~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = Mac.PDL1pos))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("Macrophages PDL1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())

wilcox.test(data3$DC~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = DC))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("DC") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())

wilcox.test(data3$DC.PDL1pos~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = DC.PDL1pos))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("DC PDL1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())

wilcox.test(data3$Microglia~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = Microglia))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("Microglia") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())

wilcox.test(data3$Microglia.PDL1pos~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = Microglia.PDL1pos))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("Microglia PDL1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())

wilcox.test(data3$CTL~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = data3$CTL))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("CTL") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$CTL.PD1pos~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = CTL.PD1pos ))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("CTL PD1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$T.Helper~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = T.Helper))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("T helper") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$T.Helper.PD1~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = T.Helper.PD1))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("T helper PD1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$DNT~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = DNT))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("DNT") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$DNT.PD1pos~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = DNT.PD1pos))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("DNT PD1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$Treg~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = Treg))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("Treg") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


wilcox.test(data3$Treg.PD1pos~data3$grade,paired = FALSE)
ggplot(data3, aes(x = grade, y = Treg.PD1pos))+
  geom_violin(aes(fill = grade)) +
  xlab("Grade") + ylab("Treg PD1+") +
  theme_bw()+
  theme(panel.grid.minor = element_blank())


cell_type <- colnames(data3) 
cell_type <- cell_type[-16]
cell_type 
P_value <- c(0.03202,0.07243,0.6842,0.02045,0.1483,0.2709,0.6038,0.494, 0.1471, 0.4585,0.8286,0.6085,0.8749,0.6857,0.8)
data4 <- as.data.frame(matrix(data = NA, nrow = 15, ncol = 2))
data4[,1] <- cell_type
data4[,2] <- P_value
names(data4) <- c("cell_type","P_value")


data4$class <- as.data.frame(as.matrix(seq(1:15)))
colnames(data4[,3]) <- "class"
plot(x = data4$class, y = data4$P_value, xlab = "cell type",ylim = c(0,1), xlim=c(1,15),ylab = "P value", main ="Wilcoxon test")
abline(h=0.05, lwd=2, col="red")



