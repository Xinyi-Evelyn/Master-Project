################################### File1

library("tibble")
library("dbscan")



number_within_radius1 <- function(sce_object, reference_cell_type, target_cell_type, radius) {
  
  formatted_data <- data.frame(colData(sce_object))
  
  intensity_matrix <- assay(sce_object)
  
  markers <- rownames(intensity_matrix)
  cell_ids <- colnames(intensity_matrix)
  
  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_matrix_t <- t(intensity_matrix)
  intensity_df <- data.frame(intensity_matrix_t)
  colnames(intensity_df) <- markers
  
  formatted_image2 <- spatialCoords(formatted_image)
  numrow <- nrow(formatted_image2)
  id2 <- 1:numrow
  rownames(formatted_image2) <- paste("Cell_",id2, sep="")
  
  formatted_data <- cbind(formatted_data,formatted_image2)
  formatted_data <- cbind(formatted_data,intensity_df)
  formatted_data <- formatted_data[complete.cases(formatted_data), ]
  
  #Select  reference cells
  Ref_cells <- which(formatted_data[,3]== reference_cell_type)
  reference_cells <- formatted_data[Ref_cells,]
  if (nrow(reference_cells) == 0) {
    stop("There are no reference cells found for the marker")
  }
  
  
  #select Target cells 
  
  Ttype_cells <- which(formatted_data[,3]== target_cell_type)
  target_cells <- formatted_data[Ttype_cells,]
  if (nrow(target_cells) == 0) {
    return(0)
  }
  
  #Get the coordinates to find neighbours
  reference_cell_cords <- reference_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  target_cell_cords <- target_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  
  #frNN output ids, the rowid of reference_cell_cords matches the row number of target_cell_cords
  search_result <- frNN(target_cell_cords, eps = radius, query = reference_cell_cords, sort = FALSE)
  rownums <- unique(unlist(search_result$id))
  return(length(rownums))
}
```
#
number_within_radius_dataframe1 <- function(sce_object, reference_cell, target_cell, radii){
  radius_number <- data.frame()
  i <- 1
  for (radius1 in radii) {
    result <- number_within_radius1(sce_object, reference_cell, target_cell, radius = radius1)
    radius_number[i,1] <- radius1
    radius_number[i,2] <- result
    i <- i+1
  }
  
  names(radius_number) <- c("radius","number")
  n <- length(radii)
  total_num <- radius_number[n,2]
  for (i in n:1){
    if (i==1){
      diff <- radius_number[i,2] - 0
      radius_number[i,2] <- diff
    }
    else{
      diff <- radius_number[i,2] - radius_number[i-1,2]
      radius_number[i,2] <- diff
    }
  }
  for (i in 1:n){
    radius_number[i,2] <- round(radius_number[i,2]/ total_num,3)
  }
  names(radius_number) <- c("radius","number_percentage")
  return(radius_number)
}


#
plot_number_percentage1 <- function(df){
  plot <- ggplot(data = df,aes(x=radius,y=number_percentage),mian="the percentage of the number")+
    geom_point()+
    geom_line()+
    xlab("radius")+
    ylab("number_percentage")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())+
    scale_x_continuous(limits=c(0,700), breaks=seq(0,700,100))+
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1))
  return(plot)
}

################################### File2
#
id_cells_in_cluster <- function(sce_object,clusters, T_cell){
  df_cluster <- data.frame(colData(clusters))
  no_na <- which(df_cluster[,4] != 0)
  df_cluster2 <- df_cluster[no_na,]
  id_cell_in_cluster <- which(df_cluster2[,4] != "Free_cell")
  df_cluster2 <- df_cluster2[id_cell_in_cluster,]
  
  # find id  of specific cells - T_cell in clusters 
  ids_cell_in_cluster <- rownames(df_cluster2[df_cluster2$Cell.Type == T_cell,])
  ids_cell_in_cluster1 <- gsub("Cell_", "",      ids_cell_in_cluster)
  ids_cell_in_cluster1 <- sort(as.numeric(ids_cell_in_cluster1))
  return(ids_cell_in_cluster1)
}

###
id_cells_notin_cluster <- function(sce_object,clusters, T_cell){
  df_cluster <- data.frame(colData(clusters))
  no_na <- which(df_cluster[,4] != 0)
  df_cluster2 <- df_cluster[no_na,]
  id_cell_in_cluster <- which(df_cluster2[,4] != "Free_cell")
  df_cluster2 <- df_cluster2[id_cell_in_cluster,]
  
  # find id  of specific cells - T_cell in clusters 
  ids_cell_in_cluster <- rownames(df_cluster2[df_cluster2$Cell.Type == T_cell,])
  
  # find id  of specific cells - T_cell not in clusters 
  data_image <- data.frame(colData(sce_object))
  ids_cell <- rownames(data_image[data_image$Cell.Type == T_cell,])
  ids_cell_not_in_cluster <- setdiff(ids_cell,ids_cell_in_cluster)
  
  ids_cell_not_in_cluster1 <- gsub("Cell_", "", ids_cell_not_in_cluster)
  ids_cell_not_in_cluster1 <- sort(as.numeric(ids_cell_not_in_cluster1))
  return(ids_cell_not_in_cluster1)
}

################################### File3
library("reshape2")

##
number_within_radius3 <- function(sce_object,Id_cells,target_cell_type,radius) {
  
  formatted_data <- data.frame(colData(sce_object))
  intensity_matrix <- assay(sce_object)
  
  markers <- rownames(intensity_matrix)
  cell_ids <- colnames(intensity_matrix)
  
  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_matrix_t <- t(intensity_matrix)
  intensity_df <- data.frame(intensity_matrix_t)
  colnames(intensity_df) <- markers
  
  formatted_image2 <- spatialCoords(formatted_image)
  numrow <- nrow(formatted_image2)
  id2 <- 1:numrow
  rownames(formatted_image2) <- paste("Cell_",id2, sep="")
  
  formatted_data <- cbind(formatted_data,formatted_image2)
  formatted_data <- cbind(formatted_data,intensity_df)
  formatted_data <- formatted_data[complete.cases(formatted_data), ]
  
  #Select  reference cells
  reference_cells <- formatted_data[Id_cells,]
  if (nrow(reference_cells) == 0) {
    stop("There are no reference cells found for the marker")
  }
  
  
  #select Target cells 
  Ttype_cells <- which(formatted_data[,3]== target_cell_type)
  target_cells <- formatted_data[Ttype_cells,]
  if (nrow(target_cells) == 0) {
    return(0)
  }
  
  #Get the coordinates to find neighbours
  reference_cell_cords <- reference_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  target_cell_cords <- target_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  #frNN output ids, the rowid of reference_cell_cords matches the row number of target_cell_cords
  search_result <- frNN(target_cell_cords, eps = radius, query = reference_cell_cords, sort = FALSE)
  rownums <- unique(unlist(search_result$id))
  return(length(rownums))
}

##
number_within_radius_dataframe3 <- function(sce_object,Id_cells, target_cell_type, radius){
  radius_number <- data.frame()
  i <- 1
  for (radii in radius) {
    result <- number_within_radius3(sce_object,Id_cells,target_cell_type, radius = radii)
    radius_number[i,1] <- radii
    radius_number[i,2] <- result
    i <- i+1
  }
  
  names(radius_number) <- c("radius","number")
  
  total_num <- radius_number[9,2]
  for (i in 9:1){
    if (i==1){
      diff <- radius_number[i,2] - 0
      radius_number[i,2] <- diff
    }
    else{
      diff <- radius_number[i,2] - radius_number[i-1,2]
      radius_number[i,2] <- diff
    }
  }
  for (i in 1:9){
    radius_number[i,2] <- round(radius_number[i,2]/ total_num,3)
  }
  names(radius_number) <- c("radius","number_percentage")
  return(radius_number)
}

#
plot_number_percentage3 <- function(df1,df2){
  colnames(df1)[2] <- 'in clusters'
  colnames(df2)[2] <- 'not in clusters'
  df <- merge(df1,df2)
  data <- melt(df, id = "radius")
  
  plot <- ggplot(data,aes(x=radius,y=value, group = variable,
                          color = variable,
                          shape = variable),mian="the percentage of tumor cells in different radius")+
    geom_point()+
    geom_line()+
    xlab("radius")+
    ylab("the percentage of tumor cells")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())+
    scale_x_continuous(limits=c(0,700), breaks=seq(0,700,100))+
    scale_y_continuous(limits=c(0,0.5), breaks=seq(0,1,0.1))
  return(plot)
}

################################### File4
#
id_cells_in_cluster4 <- function(clusters){
  df_cluster <- data.frame(colData(clusters))
  no_na <- which(df_cluster[,4] != 0)
  df_cluster2 <- df_cluster[no_na,]
  id_cell_in_cluster <- which(df_cluster2[,4] != "Free_cell")
  df_cluster2 <- df_cluster2[id_cell_in_cluster,]
  
  # find id  of specific cells - T_cell in clusters 
  ids_cell_in_cluster <- rownames(df_cluster2)
  ids_cell_in_cluster1 <- gsub("Cell_", "", ids_cell_in_cluster)
  ids_cell_in_cluster1 <- sort(as.numeric(ids_cell_in_cluster1))
  return(ids_cell_in_cluster1)
}



#
number_within_radius4 <- function(sce_object,Id_cells,target_cell_type,radius) {
  
  formatted_data <- data.frame(colData(sce_object))
  intensity_matrix <- assay(sce_object)
  
  markers <- rownames(intensity_matrix)
  cell_ids <- colnames(intensity_matrix)
  
  rownames(intensity_matrix) <- NULL
  colnames(intensity_matrix) <- NULL
  intensity_matrix_t <- t(intensity_matrix)
  intensity_df <- data.frame(intensity_matrix_t)
  colnames(intensity_df) <- markers
  
  formatted_image2 <- spatialCoords(formatted_image)
  numrow <- nrow(formatted_image2)
  id2 <- 1:numrow
  rownames(formatted_image2) <- paste("Cell_",id2, sep="")
  
  formatted_data <- cbind(formatted_data,formatted_image2)
  formatted_data <- cbind(formatted_data,intensity_df)
  formatted_data <- formatted_data[complete.cases(formatted_data), ]
  
  #Select  reference cells
  reference_cells <- formatted_data[Id_cells,]
  if (nrow(reference_cells) == 0) {
    stop("There are no reference cells found for the marker")
  }
  
  
  #select Target cells 
  Ttype_cells <- which(formatted_data[,3]== target_cell_type)
  target_cells <- formatted_data[Ttype_cells,]
  if (nrow(target_cells) == 0) {
    return(0)
  }
  
  #Get the coordinates to find neighbours
  reference_cell_cords <- reference_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  target_cell_cords <- target_cells[,c("Cell.X.Position", "Cell.Y.Position")]
  #frNN output ids, the rowid of reference_cell_cords matches the row number of target_cell_cords
  search_result <- frNN(target_cell_cords, eps = radius, query = reference_cell_cords, sort = FALSE)
  rownums <- unique(unlist(search_result$id))
  return(length(rownums))
}

#
number_within_radius_dataframe4 <- function(sce_object,clusters,target_cell_type,radius){
  Id_cells <- id_cells_in_cluster(clusters)
  
  radius_number <- data.frame()
  i <- 1
  for (radii in radius) {
    result <- number_within_radius4(sce_object,Id_cells,target_cell_type, radius = radii)
    radius_number[i,1] <- radii
    radius_number[i,2] <- result
    i <- i+1
  }
  
  names(radius_number) <- c("radius","number")
  
  for (i in 5:1){
    if (i==1){
      diff <- radius_number[i,2] - 0
      radius_number[i,2] <- diff
    }
    else{
      diff <- radius_number[i,2] - radius_number[i-1,2]
      radius_number[i,2] <- diff
    }
  }
  
  names(radius_number) <- c("radius","cell_number")
  return(radius_number)
}

#
plot_number <- function(df1,df2,df3){
  colnames(df1)[2] <- 'DC'
  colnames(df2)[2] <- 'Mac'
  colnames(df3)[2] <- 'Microglia'
  df <- merge(df1,df2)
  df <- merge(df,df3)
  data <- melt(df, id = "radius")
  
  plot <- ggplot(data,aes(x=radius,y=value, group = variable,
                          color = variable,
                          shape = variable),mian="the number of cells in different radius")+
    geom_point()+
    geom_line()+
    xlab("radius")+
    ylab("the number of cells")+
    theme_bw()+
    theme(panel.grid.minor = element_blank())+
    scale_x_continuous(limits=c(0,200), breaks=seq(0,200,50))
  return(plot)
}

