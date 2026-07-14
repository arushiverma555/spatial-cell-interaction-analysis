#load required packages
library(tidyverse)

#read in the cellular network dataset
tda <- read.csv("tda_ds.csv")

#Renamed the epithelial cells to tumor cells 
tda <- tda %>%
  mutate(CellType = ifelse(CellType == "Epithelial", "Tumor", CellType))

#chose 1 tissue core for testing 
core_id <- "1_001_01-S-17-02150"

#keep only cells from the first core and the important protein columns 
core <- tda %>%
  filter(TMA_ID == core_id)%>%
  select(Cell_ID:jm)
#need to know if there are any proteins we shouldn't examine so we can edit 

#function used in protein interactions to get close and far cells (based on radius)
#returns mean protein values for the close and far cells. Not as useful as the protein_t_test function
#cellsA is original cell type (ex: "Tumor"), cellsB is other cell type (ex: "B")
protein_means <- function(tma_ID, cellsA, cellsB, typeA, radius){
  #makes logical vector for if a cell has a neighbor or not (type B cell within radius or not) 
  has_neighborA <- logical(nrow(cellsA))
  r2 <- radius^2
  
  #iterates over cellsA and finds the distance from the cell to all cells in cellsB
  #if there is a cell within radius, adds true to has_neighbor, serving as an index of cells with neighbors
  for (i in seq_len(nrow(cellsA))) {
    x1 <- cellsA$im[i]; y1 <- cellsA$jm[i]
    d2 <- (cellsB$im - x1)^2 + (cellsB$jm - y1)^2
    has_neighborA[i] <- any(d2 <= r2)
  }
  
  #filters cellsA by has_neighbors to get close and far cells 
  close_cellsA <- cellsA[has_neighborA, , drop = FALSE]
  far_cellsA <- cellsA[!has_neighborA, , drop = FALSE]
  
  #filters close and far cells for the protein columns and calculates the mean values for the protein expression levels
  close_protein_mean <- close_cellsA %>%
    summarise(across(CD45:PM4, \(x) mean(x, na.rm = TRUE)))
  far_protein_mean <- far_cellsA %>%
    summarise(across(CD45:PM4, \(x) mean(x, na.rm = TRUE)))
  
  #sets the protein means to 0 if there are no cells in close_cells or far cells (so it isn't all NA if there are no cells)
  if(nrow(close_cellsA) == 0){
    close_protein_mean[1,] <- numeric(ncol(close_protein_mean))
  }
  if(nrow(far_cellsA) == 0){
    far_protein_mean[1,] <- numeric(ncol(far_protein_mean))
  }
  
  #creating tibble to represent the info for the close/far cells
  #contains TMA ID, cell types, counts and whether the cells types are close/far
  protein_tibble <- tibble(
    TMA_ID = tma_ID,
    CellType = typeA,
    OtherCellType = unique(cellsB$CellType),
    Count = c(nrow(close_cellsA),nrow(far_cellsA)),
    CloseOrFar = c("Close", "Far")
  )
  
  #combining the mean protein expression data for the close and far cells into one tibble by binding rows
  protein_mean_data <- bind_rows(close_protein_mean, far_protein_mean)
  
  #returns the cell info bound with the protein data into a formatted tibble 
  bind_cols(protein_tibble, protein_mean_data)
}

#more useful function for calculating mean protein expression levels based on cell proximity to other cell type
#also calculates t-test and p-test and this info is returned for each protein type 
#used in protein interactions function to get both cell type data because it only calculates data for cellsA (one way)
protein_t_test <- function(tma_ID, cellsA, cellsB, typeA, radius){
  #makes logical vector representing the index for each of the cells in cellsA 
  #true value in the vector means that the associated cell has a neighbor within its radius 
  has_neighborA <- logical(nrow(cellsA))
  r2 <- radius^2
  
  #iterates over the cells in cellsA and calculates the distance for each cell to all the cells in cellsB
  #if the cell in cellA is less than radius distance away from a cell in cellsB, it adds true to the has_neighbor index
  for (i in seq_len(nrow(cellsA))) {
    x1 <- cellsA$im[i]; y1 <- cellsA$jm[i]
    d2 <- (cellsB$im - x1)^2 + (cellsB$jm - y1)^2
    has_neighborA[i] <- any(d2 <= r2)
  }
  
  #filters the cellsA vector into close and far cells based on the has_neighbor index vector 
  #selects only the protein data for the cells
  close_cellsA <- cellsA[has_neighborA, , drop = FALSE]%>%
    select(CD45:PM4)
  far_cellsA <- cellsA[!has_neighborA, , drop = FALSE] %>%
    select(CD45:PM4)
  
  #making the counts of the cells by using nrow() on the close and far cells dataframes
  t_test_count <-c(nrow(close_cellsA),nrow(far_cellsA))
  
  #getting t-test data and formatting in a tibble 
  #this gets the protein names and the amount of proteins is used to preallocate memory for the protein means and t/p statistics
  pnames <- colnames(close_cellsA)
  P <- length(pnames)
  
  # preallocating for the means
  mean_close <- numeric(P)
  mean_far   <- numeric(P)
  
  #preallocating for the t/p statistics, initially set to NA so that if the t-test doesn't apply, it is already set to NA 
  tstat      <- rep(NA_real_, P)
  pval       <- rep(NA_real_, P)
  
  #creating a logical vector for whether or not the t-test applies
  #if the cell count is 0 or 1, the t-test does not apply and false is put into the vector 
  do_test <- !(t_test_count[1] %in% c(0, 1) || t_test_count[2] %in% c(0, 1))
  
  #iterates over all of the proteins and calculates the close and far data for each protein as well as t/p statistics
  for (j in seq_len(P)) {
    #determines if the t-test should be calculated or not based on do_test logical vector
    if (do_test) {
      #makes sure that the variance among the proteins is not 0 (will throw an error in the t-test)
      if(var(close_cellsA[[j]]) != 0 | var(far_cellsA[[j]]) != 0){
        #calculates the t-test for each protein by using the close and far cell columns
        tt <- t.test(close_cellsA[[j]], far_cellsA[[j]])
        
        #estimates are named; unname() keeps it numeric
        #adds the data to the corresponding index in the pre-allocated vector 
        mean_close[j] <- unname(tt$estimate[1])
        mean_far[j]   <- unname(tt$estimate[2])
        tstat[j] <- unname(tt$statistic)
        pval[j]  <- tt$p.value
      }
      #if the variance is 0, the tstat and pval are not updated and only the mean is calculated
      else{
        mean_close[j] <- mean(close_cellsA[[j]], na.rm = TRUE)
        mean_far[j]   <- mean(far_cellsA[[j]], na.rm = TRUE)
      }
    }
    #if there are only 0 or 1 cells in the close/far cells vector this if statement adds an appropriate
    #value for the means and avoids the t.test function
    else{
      #if there are no cells in close or far cells, the respective protein mean is set to 0 while the other mean is calculated
      if(nrow(close_cellsA) == 0){
        mean_close[j] <- 0
        mean_far[j] <- mean(far_cellsA[[j]], na.rm = TRUE)
      }else if(nrow(far_cellsA) == 0){
        mean_close[j] <- mean(close_cellsA[[j]], na.rm = TRUE)
        mean_far[j] <- 0
        #if there is only 1 cell in close or far cells, the respective protein mean is set 
        #to the value for the 1 cell and the other mean is calculuated
      }else if(nrow(close_cellsA) == 1){
        mean_close[j] <- close_cellsA[[j]]
        mean_far[j] <- mean(far_cellsA[[j]], na.rm = TRUE)
      }else{
        mean_far[j] <- far_cellsA[[j]]
        mean_close[j] <- mean(close_cellsA[[j]], na.rm = TRUE)
      }
    }
  }
  
  #builds the protein result tibble in one shot (2 rows: Close/Far)
  #iterates over the preallocated vectors and concatenates them so that a tibble can be formed
  t_test_tbl <- tibble::as_tibble(
    #this is used to set the names after the tibble is constructed
    setNames(
      #generating the list
      lapply(seq_len(P), function(j) {
        #stacks the close and far protein means in one column, the t-stats (both identical for the two rows) in another column
        #and the p-value in the third column (also identical for the rosw)
        list(
          c(mean_close[j], mean_far[j]),
          c(tstat[j], tstat[j]),
          c(pval[j],  pval[j])
        )
        #unlists the list so it becomes a total data frame
      }) |> unlist(recursive = FALSE),
      #makes the names for the three columns generated 
      as.vector(rbind(
        paste(pnames, "mean"),
        paste(pnames, "t-test"),
        paste(pnames, "p-value")
      ))
    )
  )
  
  #generating the cell info for the close and far cell data 
  #contains the TMA_ID, the cell type, the neighboring celltype, the count of close and far cells
  #and the cell distance (far or close). Then it gets bound to the t-test tibble and returned
  protein_tibble <- tibble(
    TMA_ID = tma_ID,
    CellType = typeA,
    OtherCellType = unique(cellsB$CellType),
    InteractionCount = t_test_count,
    CellDistance = c("Close","Far")
  )
  bind_cols(protein_tibble, t_test_tbl)
}

#testing 
typeA <- "CD4_Treg"
typeB <- "Macrophage"
tma_ID <- unique(core$TMA_ID)
cellsA <- core %>% filter(CellType == typeA)
#Get all cells belonging to cell type B
cellsB <- core %>% filter(CellType == typeB)
#print(protein_means(tma_ID,cellsA,cellsB,typeA,50))
test_t <- protein_t_test(tma_ID,cellsA, cellsB, typeA,50)

#calculate protein interactions between 2 cell types, 
#ex: Tumor vs CD8_T
#Two cells are considered interacting if they are within 50 microns of each other
#output:
# 4 row tibble containing TMA_ID, cell type of the original and neighboring cells
# counts of close and far cells, mean protein expression and t/p statistic for each protein (1 row close, 1 row far) 
# the columns are the mean protein expression, the t-test and the p-value for the respective protein
protein_interactions <- function(core, typeA, typeB, radius = 50) {
  #Get all cells belonging to cell type A
  cellsA <- core %>% filter(CellType == typeA)
  #Get all cells belonging to cell type B
  cellsB <- core %>% filter(CellType == typeB)
  
  #since protein_t_test function is directed, this function does A to B and B to A and binds outputs
  cellA_protein <- protein_t_test(unique(core$TMA_ID),cellsA, cellsB, typeA, radius)
  cellB_protein <- protein_t_test(unique(core$TMA_ID),cellsB, cellsA, typeB, radius)
  protein_data <-bind_rows(cellA_protein,cellB_protein)
}

#generate all cell-type combinations within a single core and return the protein interaction data for each combo
#ex: tumor - CD8-T, tumor - macrophage, CD8-T - macrophage etc
#returns eight rows per cell-type pair:
core_protein_interactions <- function(core) {
  
  #find all unique cell types present in this core
  cell_types <- sort(unique(core$CellType))
  
  #create every possible pair of cell types while avoiding duplicate pairs
  #ex: Tumor-CD8 included, CD8-Tumor excluded
  type_pairs <- expand.grid(
    CellTypeA = cell_types,
    CellTypeB = cell_types,
    stringsAsFactors = FALSE
  ) %>%
    filter(CellTypeA < CellTypeB)
  
  #calculation protein interaction stats for every cell type pair
  interaction_table <- purrr::map2_dfr(
    type_pairs$CellTypeA,
    type_pairs$CellTypeB,
    ~ protein_interactions(core, .x, .y, radius = 50)
  )
  interaction_table
}

#making the protein table for one core and sorting 
protein_test_table <- core_protein_interactions(core) %>%
  mutate(CellType = factor(CellType, levels = cell_order),
         OtherCellType = factor(OtherCellType, levels = cell_order)
  ) %>%
  arrange(TMA_ID, CellType, OtherCellType)

#getting the TMA IDs in the dataset
core_ids<- unique(tda$TMA_ID)

#loop through every core and build a cell-cell interaction table
#final output: 2 rows have the close and far cell data containing mean protein expression 
#and t/p statistics
all_protein_interactions <- purrr::map_dfr(
  
  #get IDs for every tissue core in dataset
  core_ids,
  function(id) {
    cat("Processing:", id, "\n")
    
    core_temp <- tda %>%
      filter(TMA_ID == id)
    #does it for the protein interactions
    core_protein_interactions(core_temp)
  }
)

#order cell types so output table is easier to read and consistent
cell_order <- c(
  "B",
  "CD4_T",
  "CD4_Treg",
  "CD8_T",
  "Tumor",
  "Fibroblast",
  "Granulocyte",
  "Macrophage",
  "Monocyte",
  "Myofibroblast"
)

#ordering alphabetically for each core
all_protein_interactions <- all_protein_interactions %>%
  mutate(
    CellType = factor(CellType, levels = cell_order),
    OtherCellType = factor(OtherCellType, levels = cell_order)
  ) %>%
  arrange(TMA_ID, CellType, OtherCellType)

#save final protein interaction table
write.csv(all_protein_interactions, "all_protein_interactions.csv", row.names = FALSE)