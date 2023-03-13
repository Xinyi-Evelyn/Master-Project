#Example: T panel LGG_290887 ####

library("SPIAT")
library("ggplot2")
library('data.table')
setwd("/home/dlabrstudio/working_dir/xinyi/T_LGG_code")


image <- readRDS(file="TLGGTI-290887-1.rds")
markers <- c("DAPI","CD8","FOXP3","PD1","CD4","CD3","GFAP")


##################################Clean the data and define the cell type


#Remove all cells without DAPI density
DAPI_non_zero_rows <- which(image[,10] != 0)
image <- image[DAPI_non_zero_rows,]


#Label markers of each cell
class_image <- image[,c(10,17,20,27,32,37,42)]
names(class_image) <- markers
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
phenotypes = c("GFAP","CD8,CD3","CD8,PD1,CD3", "CD4,CD3", "PD1,CD4,CD3","CD3","PD1,CD3","FOXP3,CD4,CD3","FOXP3,PD1,CD4,CD3")
names = c("Tumour", "CTL", "CTL PD1pos", "T helper", "T helper PD1pos","DNT","DNT PD1pos","Treg","Treg PD1pos")

all_phenotypes <- unique(class_image$phenotype)
inter_phenotypes <- intersect(phenotypes,all_phenotypes)
inter_phenotypes
names <- names[match(inter_phenotypes, phenotypes)]

for (i in 1:length(all_phenotypes)){
  class_image$celltype[class_image$phenotype == inter_phenotypes[i]] <- names[i]
}


# clean the undefined cells
phenotype_non_undefined_rows <- which(class_image$celltype != "Undefined")
class_image <- class_image[class_image$celltype != "Undefined",]
image <- image[phenotype_non_undefined_rows,]
nrow(class_image)
nrow(image)
type <- table(class_image$celltype)


##################################The abundance of each type of immune cell

#The percentage that each subtype of T cells occupied all of T cells 
df_type <- as.data.frame(type)
df_type <- df_type[-9,]
num <- df_type[,2]
names <- as.vector(df_type[,1])
cols <-  c("violetred1","black", "cyan","pink","white","yellowgreen", "cornsilk","yellow")
per.type<- paste(round(100 * num/ sum(df_type[,2]),2),"%")
pie(num,labels= per.type,col=cols)
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
  labs(title = paste("The number of each subtype of T cells"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_text(aes(label=number), vjust=1.6, color="black", size=3.5)


#The density of each cells  number/mm^2
#The density of each subtype of T cells occupied all of T cells
ggplot(df_type, aes(x = celltype, y= cell_density_1,fill = celltype)) + geom_bar(stat="identity") + 
  ylab("density")+
  xlab("cell type")+
  labs(title = paste("The density of each subtype of T cells"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_text(aes(label=cell_density_1), vjust=1.6, color="black", size=3.5)



##################################Creating the SPIAT Objective

coord_x <- image$XMin
coord_y <- image$YMin
phenotypes <- class_image$phenotype
intensity_columns <- image[,c(11, 18, 21, 26, 31, 36, 41)]
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
                                    categories =  c("GFAP","CD8,CD3","CD8,PD1,CD3", "CD4,CD3", "PD1,CD4,CD3","CD3","PD1,CD3","FOXP3,CD4,CD3","FOXP3,PD1,CD4,CD3"), 
                                    category_colname = "Phenotype",
                                    names = c("Tumour", "CTL", "CTL PD1pos", "T helper", "T helper PD1pos","DNT","DNT PD1pos",'Treg','Treg PD1pos'),new_colname = "Cell.Type")

p_cells <- calculate_cell_proportions(formatted_image, 
                                      reference_celltypes=NULL, 
                                      feature_colname ="Cell.Type",
                                      celltypes_to_exclude = "Tumour",
                                      plot.image = TRUE)
plot_cell_percentages(cell_proportions = p_cells, 
                      cells_to_exclude = "Tumour", cellprop_colname="Proportion_name")




##############Visualize the location of each cell
my_colors <- c("pink","red","blue","darkcyan", "black","purple","Green","yellow","orange")
plot_cell_categories(formatted_image, c("Tumour", "CTL", "CTL PD1pos", "T helper", "T helper PD1pos","DNT","DNT PD1pos","Treg","Treg PD1pos"), my_colors, "Cell.Type")


# The percentage of target cells within the increasing distance from the reference cells

df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="DNT",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="T helper",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="CTL",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="DNT PD1pos",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="T helper PD1pos",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)
#

df <- number_within_radius_dataframe1(formatted_image, reference_cell= "Tumour",target_cell="CTL PD1pos",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell="DNT" ,target_cell="Tumour",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "T helper",target_cell="Tumour",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "CTL",target_cell="Tumour",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell="DNT PD1pos",target_cell="Tumour",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell="T helper PD1pos",target_cell="Tumour",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)

#
df <- number_within_radius_dataframe1(formatted_image, reference_cell= "CTL PD1pos",target_cell="Tumour",radii =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage1(df)



#########Identify the spatial T clusters
clusters <- identify_neighborhoods(formatted_image, 
                                   method = "hierarchical",
                                   min_neighborhood_size = 30,
                                   cell_types_of_interest = c("T helper","DNT","CTL","CTL PD1pos","T helper PD1pos","DNT PD1pos"),radius =200, feature_colname = "Cell.Type")

#identift the composition of each cluster
neighorhoods_vis <- composition_of_neighborhoods(clusters, feature_colname = "Cell.Type")
neighorhoods_vis <- neighorhoods_vis[neighorhoods_vis$Total_number_of_cells >=5,]
plot_composition_heatmap(neighorhoods_vis, feature_colname="Cell.Type")



#The distance between tumor cells and clustering and non-clustering T cells

#DNT

ids_cell_in_cluster <- id_cells_in_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "DNT")
ids_cell_not_in_cluster <- id_cells_notin_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "DNT")
df1 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_in_cluster,target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
df2 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_not_in_cluster, target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage3(df1,df2)


#"T helper"

ids_cell_in_cluster <- id_cells_in_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "T helper")
ids_cell_not_in_cluster <- id_cells_notin_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "T helper")
df3 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_in_cluster,target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
df4 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_not_in_cluster, target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage3(df3,df4)


#"CTL"

ids_cell_in_cluster <- id_cells_in_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "CTL")
ids_cell_not_in_cluster <- id_cells_notin_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "CTL")
df5 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_in_cluster,target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
df6 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_not_in_cluster, target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage3(df5,df6)


#"DNT PD1pos"


ids_cell_in_cluster <- id_cells_in_cluster (sce_object = formatted_image,clusters = clusters, T_cell ="DNT PD1pos")
ids_cell_not_in_cluster <- id_cells_notin_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "DNT PD1pos")
df7 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_in_cluster,target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
df8 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_not_in_cluster, target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage3(df7,df8)


###"T helper PD1pos"


ids_cell_in_cluster <- id_cells_in_cluster (sce_object = formatted_image,clusters = clusters, T_cell ="T helper PD1pos")
ids_cell_not_in_cluster <- id_cells_notin_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "T helper PD1pos")
df9 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_in_cluster,target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
df10 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_not_in_cluster, target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage3(df9,df10)


#"CTL PD1pos"

ids_cell_in_cluster <- id_cells_in_cluster (sce_object = formatted_image,clusters = clusters, T_cell ="CTL PD1pos")
ids_cell_not_in_cluster <- id_cells_notin_cluster (sce_object = formatted_image,clusters = clusters, T_cell = "CTL PD1pos")
df11 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_in_cluster,target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
df12 <- number_within_radius_dataframe3(formatted_image,Id_cells =ids_cell_not_in_cluster, target_cell_type="Tumour",radius =  c(30,50,100,200,300,400,500,600,700))
plot_number_percentage3(df11,df12)

