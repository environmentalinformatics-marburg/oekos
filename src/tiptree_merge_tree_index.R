fp_data <- ("E:/analysis/tiptree/tiptree/data/")

library(rgdal)
library(rgeos)
library(raster)
library(sp)

# Function sdist ---------------------------------------------------------------
sdist <- function(east_i, north_i, east_j, north_j){
  sqrt((east_i - east_j)**2 + (north_i - north_j)**2)
}
  
# Function HgCI ----------------------------------------------------------------
hgci <- function(bhd_i, bhd_j, pos_ij){
  index <- lapply(seq(length(bhd_j)), function(x){
    sd_ij <- sdist(east_i = pos_ij$EAST[1], north_i = pos_ij$NORTH[1], 
                   east_j = pos_ij$EAST[x+1], north_j = pos_ij$NORTH[x+1])  
    bhd_j[x] / (bhd_i / sd_ij)
  })
  return(sum(unlist(index)))
}
  

# Function Schuetz index -------------------------------------------------------
schtzi <- function(height_i, carea_i,height_j, carea_j, pos_ij){
  index <- lapply(seq(length(height_j)), function(x){
    sd_ij <- sdist(east_i = pos_ij$EAST[1], north_i = pos_ij$NORTH[1], 
                   east_j = pos_ij$EAST[x+1], north_j = pos_ij$NORTH[x+1])  
    e_ij <- sd_ij - sqrt(carea_i / pi) + sqrt(carea_j[x] / pi)
    thv_ij <- 0.5 * (sqrt(carea_i / pi) + sqrt(carea_j[x] / pi)) +
      0.65 * (height_j[x] - height_i)
    if(e_ij <= thv_ij){
      si <- 0.65 * (height_j[x] - height_i) / 
        (sqrt(carea_i / pi) + sqrt(carea_j[x] / pi)) + 0.5 - 
        e_ij / (sqrt(carea_i / pi) + sqrt(carea_j[x] / pi))
    } else {
      si <- 0.0
    }
    return(si)
  })
  return(sum(unlist(index)))
}


# Read datasets ----------------------------------------------------------------
lidar_trees <- readOGR(paste0(fp_data, "tiptree_2015_lidar_tree_locations.shp"),
                       layer = "tiptree_2015_lidar_tree_locations")

field_info <- read.table(paste0(fp_data, "tiptree_zuordnung_final.csv"),
                         sep = ";", header = TRUE)

head(field_info)
head(lidar_trees)


# Merge datasets ---------------------------------------------------------------
# Convert lidar_trees to data frame, merge it with the field data
lidar_trees_df <- as.data.frame(lidar_trees)
merged_data <- merge(field_info, lidar_trees_df, 
                     by.x = "OID_DM", by.y = "OID_1")

nrow(merged_data)
nrow(field_info)

merged_data$EAST <- merged_data$coords.x1
merged_data$NORTH <- merged_data$coords.x2

# Create shape dataset ---------------------------------------------------------
head(merged_data)
coordinates(merged_data) <- ~coords.x1+coords.x2 
projection(merged_data) <- projection(lidar_trees)


# Compute HgCI using BHD and height --------------------------------------------
# Create buffer of 8 m and extract neighbouring trees
tree_buffer <- gBuffer(merged_data, width = 8.0, byid = TRUE)
tree_ngbh <- gContains(tree_buffer, merged_data, 
                       returnDense = FALSE, byid = TRUE)

# Compute HgCI index
hgci_index <- lapply(seq(length(tree_ngbh)), function(x){
  tree_i <- merged_data@data[x,]
  tree_j <- merged_data@data[tree_ngbh[[x]],]
  tree_j <- tree_j[tree_j$GID_1 != tree_i$GID_1,]
  pos_ij <- rbind(tree_i[, c("EAST", "NORTH")], tree_j[, c("EAST", "NORTH")])
  if(nrow(tree_j) > 0){
    index_bhd <- hgci(bhd_i = tree_i$BHD, bhd_j = tree_j$BHD, pos_ij = pos_ij)
    index_height <- hgci(bhd_i = tree_i$HEIGHT, 
                         bhd_j = tree_j$HEIGHT, pos_ij = pos_ij)
  } else {
    index_bhd <- 0.0
    index_height <- 0.0
  }
  return(data.frame(HGCI_BHD = index_bhd,
                    HGCI_HEIGHT = index_height))
})
merged_data@data <- cbind(merged_data@data, do.call("rbind", hgci_index))


# Compute Schuetz index using --------------------------------------------------
# Create buffer of 50 m and extract neighbouring trees
tree_buffer <- gBuffer(merged_data, width = 75.0, byid = TRUE)
tree_ngbh <- gContains(tree_buffer, merged_data, 
                       returnDense = FALSE, byid = TRUE)

# Compute Schuetz index
index_schtz <- lapply(seq(length(tree_ngbh)), function(x){
  tree_i <- merged_data@data[x,]
  tree_j <- merged_data@data[tree_ngbh[[x]],]
  tree_j <- tree_j[tree_j$GID_1 != tree_i$GID_1,]
  pos_ij <- rbind(tree_i[, c("EAST", "NORTH")], tree_j[, c("EAST", "NORTH")])
  if(nrow(tree_j) > 0){
    index <- schtzi(height_i = tree_i$HEIGHT, carea_i = tree_i$CrownArea, 
                    height_j = tree_j$HEIGHT, carea_j = tree_j$CrownArea, 
                    pos_ij = pos_ij)
  } else {
    index <- NA
  }
  return(index)
})
merged_data@data$SCHTZI <- unlist(index_schtz)


# Write data to shape file  ----------------------------------------------------
writeOGR(merged_data, 
         dsn = paste0(fp_data, "tiptree_2015_merged_tree_indices.shp"), 
         layer = "tiptree_2015_merged_tree_indices", driver = "ESRI Shapefile", 
         overwrite = TRUE)
