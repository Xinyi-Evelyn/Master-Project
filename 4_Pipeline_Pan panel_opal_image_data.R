
#Example: Pan panel LGG_290887 ####
library("SPIAT")
library("ggplot2")
library('data.table')
setwd("/home/dlabrstudio/working_dir/xinyi/Pan LGG")


image <- readRDS(file="PanLGGTI-290887-1.rds")
markers <- c("DAPI","PDL1","CD11c","CD68","TMEM119","CD3","GFAP")



##################################Clean the data and define the cell type


#Remove all cells without DAPI density
DAPI_non_zero_rows <- which(image[,10] != 0)
image <- image[DAPI_non_zero_rows,]
nrow(image)
Region <- which(image[,2] != "Cerebellum")
image <- image[Region,]
nrow(image)
class_image <- image[,c(10,17,22,27,32,37,42)]
names(class_image) <- markers

#Remove the cells just with DAPI marker
class_image$phenotype <- ''
class_image$phenotype[class_image$DAPI == 1] <- "OTHER"
for (i in 2:7){
  class_image$phenotype[class_image[,i] == 1]<- paste(class_image$phenotype[class_image[,i] == 1], markers[i], sep = ",")
}

#Remove the cells just with DAPI marker
class_image$phenotype[class_image$phenotype == "OTHER"] <- 0
class_image$phenotype <- gsub("OTHER,", "", class_image$phenotype)

#check the number of each phenotype after cleaning cells
table(class_image$phenotype)

#remove abnormal markers
class_image$phenotype <- gsub(",GFAP", "", class_image$phenotype)
table(class_image$phenotype)
phenotype_non_zero_rows <- which(class_image$phenotype != 0)
class_image <- class_image[phenotype_non_zero_rows,]
image <- image[phenotype_non_zero_rows,]
nrow(image)
nrow(class_image)



#Define the cell type
class_image$celltype <- "Undefined"

phenotypes =  c("GFAP","CD3","CD68","CD11c,CD68","PDL1,CD68","PDL1,CD11c,CD68","CD11c", "PDL1,CD11c",
                "TMEM119","CD68,TMEM119","CD11c,TMEM119","CD11c,CD68,TMEM119",
                "PDL1,TMEM119","PDL1,CD68,TMEM119","PDL1,CD11c,TMEM119","PDL1,CD11c,CD68,TMEM119")

names = c("Tumour", "T cell", "Mac","Mac", "Mac PDL1pos","Mac PDL1pos","DC", "DC PDL1pos",
          "Microglia","Microglia","Microglia","Microglia",
          "Microglia PDL1pos","Microglia PDL1pos","Microglia PDL1pos","Microglia PDL1pos")

all_phenotypes <- unique(class_image$phenotype)
inter_phenotypes <- intersect(phenotypes,all_phenotypes)
inter_phenotypes
names <- names[match(inter_phenotypes, phenotypes)]

for (i in 1:length(all_phenotypes)){
  class_image$celltype[class_image$phenotype == inter_phenotypes[i]] <- names[i]
}


#clean the undefined cells
phenotype_non_undefined_rows <- which(class_image$celltype != "Undefined")
class_image <- class_image[class_image$celltype != "Undefined",]
image <- image[phenotype_non_undefined_rows,]
nrow(class_image)
nrow(image)
type <- table(class_image$celltype)
class_image$phenotype[class_image$phenotype == "CD11c,CD68" ] <- "CD68"
class_image$phenotype[class_image$phenotype == "PDL1,CD11c,CD68" ] <- "PDL1,CD68"
class_image$phenotype[class_image$phenotype == "CD11c,TMEM119" | class_image$phenotype=="CD68,TMEM119"|class_image$phenotype =="CD11c,CD68,TMEM119"]  <- "TMEM119"
class_image$phenotype[class_image$phenotype == "PDL1,CD11c,TMEM119" | class_image$phenotype=="PDL1,CD68,TMEM119"|class_image$phenotype =="PDL1,CD11c,CD68,TMEM119"]  <- "PDL1,TMEM119"
phenotype <- table(class_image$phenotype)


##################################The abundance of each type of immune cell

#The percentage that each subtype of T cells occupied all of T cells 
df_type <- as.data.frame(type)
df_type <- df_type[-8,]
num <- df_type[,2]
names <- as.vector(df_type[,1])
cols <-  c("violetred1","black", "cyan","pink","white","yellowgreen", "cornsilk")
per.type<- paste(round(100 * num/ sum(df_type[,2]),2),"%")
pie(num,labels= per.type,col=cols,main = "The percentage of each type of immune cell")
legend("topright",names, cex=0.8, fill= cols)

#calculate the density of each cells
tumour <- which(class_image[,9] == "Tumour")
tumour_areas <- image[tumour,44]
tumour_areas <-  as.vector(as.matrix(tumour_areas))
length(tumour_areas)

#change the nanometer square to millimeter square

total_ares <- sum(tumour_areas)
total_ares <- total_ares/1000000
total_ares
cell_frq <- as.vector(as.matrix(df_type[,2]))
cell_frq 
cell_density <- cell_frq/total_ares
cell_density
df_type$cell_density <- as.matrix(cell_density)
colnames(df_type) <- c("celltype","number","cell_density")
df_type$cell_density_1 <- round(df_type$cell_density,3)
df_type$celltype <- reorder(df_type$celltype,df_type$cell_density)


#The number of each subtype of T cells occupied all of T cells 
df_type$celltype <- reorder(df_type$celltype,df_type$number)
ggplot(df_type, aes(x = celltype, y= number,fill = celltype)) + geom_bar(stat="identity") + 
  labs(title = paste("The number of each type of immune cells"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_text(aes(label=number), vjust=1.6, color="black", size=3.5)



##density of each cells  number/mm^2
#The density of each subtype of T cells occupied all of T cells
ggplot(df_type, aes(x = celltype, y= cell_density_1,fill = celltype)) + geom_bar(stat="identity") + 
  ylab("density")+
  xlab("cell type")+
  labs(title = paste("The density of each type of immune cells"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_text(aes(label=cell_density_1), vjust=1.6, color="black", size=3.5)



##################################Creating the SPIAT Objective
coord_x <- image$XMin
coord_y <- image$YMin
phenotypes <- class_image$phenotype


intensity_columns <- image[,c(11, 18, 23, 28, 33, 38, 43)]
intensity_matrix <- t(intensity_columns)

formatted_image <- format_image_to_spe(
  format="general",
  intensity_matrix=intensity_matrix,
  phenotypes=phenotypes,
  coord_x=coord_x,
  coord_y=coord_y)

numcol <- ncol(formatted_image)
id <- 1:numcol
colnames(formatted_image) <- paste("Cell_",id, sep="")
rownames(formatted_image) <- markers
assay(formatted_image)[,1:5]


#Specifying Cell types
formatted_image <- define_celltypes(formatted_image, 
                                    categories=c("TMEM119","GFAP","CD3","CD68","CD11c","PDL1,CD11c","PDL1,TMEM119","PDL1,CD68"), 
                                    category_colname = "Phenotype",
                                    names = c("Microglia","Tumour","T cell","Mac","DC","DC PDL1pos","Microglia PDL1pos","Mac PDL1pos"),new_colname = "Cell.Type")

p_cells <- calculate_cell_proportions(formatted_image, 
                                      reference_celltypes=NULL, 
                                      feature_colname ="Cell.Type",
                                      celltypes_to_exclude = "Tumour",
                                      plot.image = TRUE)
plot_cell_percentages(cell_proportions = p_cells, 
                      cells_to_exclude = "Tumour", cellprop_colname="Proportion_name")



##############Visualize the location of each cell
my_colors <- c("pink","red","blue","darkcyan", "yellow","purple","Green","black")
plot_cell_categories(formatted_image, c("Tumour","Microglia","T cell","Mac","DC","DC PDL1pos","Microglia PDL1pos","Mac PDL1pos"), my_colors, "Cell.Type")


# The percentage of target cells within the increasing distance from the reference cells

df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="Microglia",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="Mac",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="DC",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="T cell",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="Microglia PDL1pos",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="Mac PDL1pos",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="DC PDL1pos",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)


##########Identify the spatial T clusters
clusters <- identify_neighborhoods(formatted_image, 
                                   method = "hierarchical",
                                   min_neighborhood_size = 30,
                                   cell_types_of_interest = c("T cell"),radius = 200, feature_colname = "Cell.Type")


#The distance between APCs and clustering T cells
df1 <- number_within_radius_dataframe4(formatted_image,clusters = clusters ,target_cell_type = "DC",radius =  c(30,50,100,150,200))
df2 <- number_within_radius_dataframe4(formatted_image,clusters,target_cell_type = "Mac",radius =  c(30,50,100,150,200))
df3 <- number_within_radius_dataframe4(formatted_image,clusters= clusters ,target_cell_type = "Microglia",radius =  c(30,50,100,150,200))
plot1 <- plot_number(df1,df2,df3)
plot2 <- plot_percentge4(df1,df2,df3)
plot1
plot2














