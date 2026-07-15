library("tidyverse")
#library("readxl")

#This is just reading the CSV and filtering into cores based on TMA_ID
tda <- read.csv("tda_ds.csv")
tma_id_factors <- levels(factor(tda[["TMA_ID"]]))
tda_core1 <- tda %>%
  filter(TMA_ID == "1_001_01-S-17-02150") %>%
  select(-(CD45:PM4),-(Patient_ID:overall_survival))
tda_core2 <- tda %>%
  filter(TMA_ID == tma_id_factors[2]) %>%
  select(-(CD45:PM4),-(Patient_ID:overall_survival))

# Compute Euclidean distances from ONE "source" cell to MANY "target" cells.
# Inputs:
#   cell1  : a 1-row data.frame representing the source cell
#   cells2 : a data.frame of target cells (one row per cell)
#   x_col  : name of the x-coordinate column (default "im")
#   y_col  : name of the y-coordinate column (default "jm")
#   id_col : name of the unique cell ID column used to name the output (default "Cell_ID")
#   max_distance : allows user to set maximum distance to be returned in the dataframe (default -1 means all cells included)
# Output:
#   A named numeric vector of distances (length = nrow(cells2)).
#   names(output) are the target cell IDs from cells2[[id_col]].
# Notes:
#   - This is vectorized: it computes all distances in one shot
#   - Assumes coordinates are numeric; if they are character/factor, convert first
distance_vec <- function(cell1, cells2, max_distance = -1, x_col = "im", y_col = "jm", id_col = "Cell_ID") {
  #Gets source cell x and y data 
  x1 <- cell1[[x_col]][[1]]
  y1 <- cell1[[y_col]][[1]]
  
  #Gets x and y data from the data frame of target cells
  x2 <- cells2[[x_col]]
  y2 <- cells2[[y_col]]
  
  #makes data frame of distances 
  d <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  
  #adds Cell ID names
  names(d) <- cells2[[id_col]]
  
  #filters d based on max_distance
  if(max_distance == -1){
    return (d)
  }else{
    d_filtered <- d[d <= max_distance]
    return (d_filtered)
  }
}

# Return the set of cell types present in the dataset, in a consistent order.
# Inputs:
#   core : data frame containing all of the core data (needs to have "CellType" column)
# Output:
#   A character vector of cell type names, e.g. c("B", "Myeloid", "T")
cell_types <- function(core) {
  #gets the distinct cell type labels and sorts them lexicographically
  sort(unique(core$CellType))
}

# Counts how many cells of each type appear in a single core by tabulating the 'CellType' column.
# Inputs:
#   core (data.frame): Must contain a column named `CellType`
#     (character or factor), with one row per cell.
# Outputs:
#   integer vector: Counts per cell type. The vector is named using the
#     cell-type labels produced by 'table(core$CellType)'.
cell_count <- function(core) {
  as.integer(table(core$CellType))
}



# Computes the mean Euclidean (or other) distance from a single reference cell
# ('cell_vector') to cells of each other cell type in the core, returning one
# mean distance per cell type.
# This is an older function used to calculate distance, cell_type_distance is the updated version
# Inputs:
#   core (data.frame): Dataset of cells that includes CellType
#     column and im/jm for  'distance_vec()'.
#   cell_vector (list / named vector / 1-row data.frame): Representation of a
#     single reference cell. Must include a CellType entry and im/jm for distance_vec()`.
#   max_distance : allows user to set maximum distance to be returned in the dataframe (default -1 means all cells included)
# Outputs:
#   numeric named vector: Names are the other cell types found in `core`
#     (excluding cell_vector's CellType). Values are the mean distance from
#     the reference cell to all cells of that type.
# Notes / More details:
#   - Extracts the reference cell type as `cell_original <- cell_vector[["CellType"]]`.
#   - Filters `core` to exclude cells with the same `CellType` as the reference.
#   - Splits the remaining rows by `CellType` using `split(...)`, producing a
#     list of data frames (one per cell type).
#   - For each cell-type data frame `df_cells`, calls `distance_vec(cell_vector, df_cells)`
#     to get distances to each cell, then averages them with `mean()`.
#   - Uses `vapply(..., numeric(1))` to enforce a numeric scalar return per group.
#   - Requires `%>%` and `dplyr::filter` (i.e., dplyr/magrittr). Also depends on
#     a user-defined `distance_vec()` function.
#   - Edge cases:
#       * If there are no “different” cells, `split_cells` is empty and an empty
#         numeric vector is returned.
#       * If any `distance_vec()` results contain `NA`, `mean()` will return `NA`
#         unless `na.rm = TRUE` is added.
euclidean_distance <- function(core,cell_vector, max_d = -1){
  cell_original <- cell_vector[["CellType"]]
  different_cells <- core %>%
    dplyr::filter(CellType != cell_original)
  split_cells <- split(different_cells, different_cells[["CellType"]])
  vapply(split_cells, function(df_cells){
    distance_cell <- distance_vec(cell_vector, df_cells, max_distance = max_d)
    mean(distance_cell, na.rm = TRUE)
  }, numeric(1))
}


# Computes the average distance between cells of cell_type1 and cells of
# cell_type2 within a core by (1) computing distances from each `cell_type1`
# cell to all cell_type2 cells and (2) averaging those distances across
# cell_type1 cells.
# Inputs:
#   core: data.frame/tibble with a `CellType` column and im/jm for `distance_vec()`.
#   cell_type1: character scalar naming the “source” cell type.
#   cell_type2: character scalar naming the “target” cell type.
#   max_distance: numeric scalar (default -1) passed through to `distance_vec()`
#     acts as a cutoff for cells when specified
# Outputs:
#   A numeric scalar: the mean distance from `cell_type1` to `cell_type2`,
#   computed as the mean of per-`cell_type1`-cell mean distances to `cells2`.
# Notes / More details:
#   - Subsets `core` into `cells1` and `cells2` using `filter(CellType == ...)`.
#     (`cells1` is created but not used later; `idx` is used instead.)
#   - Finds row indices of `cell_type1` cells via `which(core[["CellType"]] == cell_type1)`.
#   - Stops with an error if there are no `cell_type1` cells.
#   - For each index `i` in `idx`, calls:
#       `distance_vec(core[i, , drop = FALSE], cells2, max_distance)`
#     to compute distances from that single cell to all `cell_type2` cells.
#   - Computes the mean distance per `cell_type1` cell, then averages across
#     all `cell_type1` cells.
#   - Dependencies: `%>%` and `filter()` (dplyr/magrittr), plus a user-defined
#     `distance_vec()`.
#   - Edge cases:
#       * If `cells2` has 0 rows, `distance_vec()`/`mean()` behavior may error or
#         return `NA` depending on implementation.
cell_type_distance <- function(core, cell_type1,cell_type2, max_distance = -1){
  # cells1 <- core %>% 
  #   filter(CellType == cell_type1)
  cells2 <- core %>% 
    filter(CellType == cell_type2)
  idx <- which(core[["CellType"]] == cell_type1)
  if (length(idx) == 0) stop("No cells of type: ", cell_type1)
  cells2
  per_cell <- lapply(idx, function(i) {
    distance_vec(core[i, , drop = FALSE], cells2, max_distance)
  })
  mean(unlist(per_cell), na.rm = TRUE)
  # a <- mean(vapply(per_cell, mean, numeric(1)),na.rm=TRUE)
  # b <- mean(unlist(lapply(per_cell,unlist)),na.rm=TRUE)
  # print(a)
  # print(b)
}

# Computes a full (square) matrix of mean distances between all pairs of cell
# types in `core` by calling `cell_type_distance()` for every (row type, column
# type) combination.
# Inputs:
#   core: data.frame/tibble containing a `CellType` column and the coordinate/feature
#   columns required by `cell_type_distance()` / `distance_vec()`.
#   max_distance: numeric scalar (default -1) forwarded to `cell_type_distance()`
#   (and ultimately to `distance_vec()`); if not -1, distances may be filtered by
#   the cutoff, which can make results directional and produce `NaN` for some pairs.
# Outputs:
#   A numeric matrix with dimensions (#cell types) x (#cell types). Row and column
#   names are the cell type labels from `cell_types(core)`. Entry [i, j] is the
#   result of `cell_type_distance(core, type_i, type_j, max_distance)`.
# Notes / More details:
#   - Retrieves the set/order of cell types via `cell_types(core)` and uses it to
#     size and name the matrix.
#   - Fills every cell of the matrix using nested loops over all i and j.
#   - The diagonal [i, i] is computed as well; depending on `cell_type_distance()`,
#     this may be meaningful, `NA`, or `NaN` (especially if `max_distance` filters
#     out all within-type distances).
#   - If `cell_type_distance()` is not symmetric (e.g., due to `max_distance`
#     filtering), then `distance_mat[i, j]` may differ from `distance_mat[j, i]`.
#     The commented line would force symmetry by copying values, but would overwrite
#     potentially different directed values.
distance_matrix <- function(core, max_distance = -1){
  types_cells <- cell_types(core)
  num_cell_types = length(types_cells)
  distance_mat <- matrix(ncol=num_cell_types, nrow=num_cell_types, dimnames = list(cell_types(core), cell_types(core)))
  for(i in 1:length(types_cells)){
    for(j in i:length(types_cells)){
      distance_mat[i,j] <- cell_type_distance(core, types_cells[i], types_cells[j], max_distance)
      # distance_mat[j,i] <- distance_mat[i,j]
    }
  }
  distance_mat
}



euclidean_distance(tda_core1, tda_core1[1,])
print(cell_type_distance(tda_core1,"B","Epithelial"))

cell_type_distance(tda_core1,"B","Epithelial", 10)
cell_type_distance(tda_core1,"Epithelial","B", 10)
#is now symmetric when called with distance limit because redid averages
core1 <- distance_matrix(tda_core1,50)
core2 <- distance_matrix(tda_core2)