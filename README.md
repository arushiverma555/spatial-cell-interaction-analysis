# Spatial Cell Interaction Analysis
Computational framework for identifying spatial relationships between immune and tumor cells in multiplex tissue imaging datasets using graph-based network analysis.

# Technologies 
R • tidyverse • igraph • networkD3 • Statistical Analysis • Spatial Omics

# Why This Project Matters
The spatial organization of immune cells within tumors strongly influences disease progression and response to immunotherapy.

This project develops a computational workflow to quantify cell-cell neighborhoods, construct interaction networks, and identify protein-expression differences associated with spatial proximity.

# My Contributions
- Developed computational pipeline for pairwise spatial distance calculations
- Designed graph-based interaction network generation
- Built protein interaction analysis workflow
- Created interactive network visualizations
- Optimized visualization aesthetics and edge-weight filtering
- Performed statistical comparisons between neighboring and distant cell populations

## Workflow

1. Import cell-coordinate and cell-type data
2. Calculate pairwise Euclidean distances
3. Classify cell pairs as near or far
4. Summarize interactions by cell-type pair
5. Compare protein expression between near and far populations
6. Construct weighted interaction networks
7. Generate interactive visualizations

## Computational Methods

- Pairwise Euclidean distance calculations
- Cell-neighborhood classification
- Cell-type interaction aggregation
- Protein-expression comparison
- Statistical testing
- Weighted graph construction
- Interactive network visualization

## Example Outputs

### Overall Cell-Interaction Network


### Core-Level Interaction Network


### Interaction Summary Table

