# ================================================================
# Cell-Cell Spatial Interaction Network Visualization (v3)
# ================================================================
#
# ---------- THE BIG PICTURE (read this first) -------------------
#
# Imagine a slice of tumor tissue under a microscope. It is full of
# different kinds of cells: tumor cells, immune cells (B cells, T
# cells, macrophages...), and structural cells (fibroblasts). We took
# ~150 little circular samples of tissue ("cores"), and for every core
# we know where each cell is and which of 44 proteins it is making.
#
# We want to answer two questions and draw the answers as a network:
#
#   Q1. Do certain cell types like to sit NEXT TO each other more (or
#       less) than you'd expect if they were just scattered randomly?
#       -> this becomes the ARROW COLOR (red = they clump together,
#          blue = they avoid each other).
#
#   Q2. When two cell types ARE neighbors, does that change which
#       proteins they make? (A cell might "turn on" a protein only when
#       it's next to a certain other cell.)
#       -> this becomes the ARROW THICKNESS (thicker = more proteins
#          change when they're neighbors).
#
# The network = dots (cell types) connected by arrows (relationships).
#
# ---------- HOW TO READ THE FINISHED PICTURE --------------------
#   DOT (node)   = one cell type. BIGGER dot = more of that cell type.
#   ARROW        = points FROM the cell whose proteins we measured TO
#                  the neighbor we compared against.
#   ARROW COLOR  = a diverging scale that shows direction AND size:
#                  RED   -> found together MORE than chance (attract)
#                  BLUE  -> found together LESS than chance (avoid)
#                  GREY  -> about as often as chance (neutral)
#                  DARKER shade = bigger difference from expected
#                  (light red = a little more, dark red = way more)
#   ARROW WIDTH  = how many of the 44 proteins change when they're near
#   LOOP on a dot= that cell type clumps with ITSELF (an "aggregate")
#
# ---------- WHAT CHANGED FROM THE OLD VERSION (v2) --------------
# Our mentor pointed out the old picture was misleading, mainly because
# arrow thickness was just the raw number of contacts. Tumor cells are
# SO numerous that they touched everything a lot, so every tumor arrow
# looked huge -- even though that's just because there are many of them,
# not because anything special is happening. The fix (this version):
#   1. Compare OBSERVED contacts to EXPECTED contacts (what you'd get by
#      pure chance given how common each cell type is), so abundance no
#      longer fakes a strong signal.  -> arrow COLOR
#   2. Use the number of significant proteins for arrow WIDTH instead.
#   3. Only draw an arrow when the observed-vs-expected difference is
#      real (statistically significant AND a big enough change), so the
#      picture isn't a hairball of meaningless arrows.
#   4. Curve the arrows and allow self-loops so A->B and B->A are both
#      visible and "clumps with itself" can be shown.
#   5. Dot size = number of cells (shrunk with a log, see STEP 3b).
#
# OUTPUTS (two interactive web pages you can open in a browser):
#   overall_interaction_network.html -- averaged over all cores
#   core_interaction_network.html    -- a single core (change core_id)
# ================================================================


# ================================================================
# STEP 0: make sure the helper toolkits ("packages") are installed.
# A package is just a bundle of pre-written code someone else made so
# we don't have to reinvent it. This loop checks each one and installs
# it only if it's missing.
# ================================================================
packages <- c("tidyverse", "igraph", "visNetwork", "htmlwidgets")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {       # "is it already here?"
    install.packages(pkg, repos = "https://cloud.r-project.org")  # if not, get it
  }
}

# ================================================================
# STEP 1: load those toolkits so we can use them in this session.
# ================================================================
library(tidyverse)   # tools for slicing/summarizing tables (filter, group_by...)
library(igraph)      # does the math of arranging dots into a network layout
library(visNetwork)  # turns our dots+arrows table into an interactive picture
library(htmlwidgets) # saves that picture as a shareable .html file

# ================================================================
# STEP 2: read in the data table your partner produced.
# check.names = FALSE tells R "don't rename my columns" -- we want the
# column names to stay exactly as written (e.g. "CD45 t_test"), not be
# auto-changed to "CD45.t_test".
# ================================================================
df <- read.csv("all_protein_interactions_with_counts.csv", check.names = FALSE)

# ----------------------------------------------------------------
# TUNABLE SETTINGS ("knobs"). Change these numbers to change behavior;
# you never have to touch the code below.
# ----------------------------------------------------------------

# --- Knobs for deciding when a PROTEIN really changed near a neighbor.
# The data already ran a "t-test" for each protein. A t-test asks:
# "are two groups of numbers actually different, or could the gap just
# be luck?" It gives back two numbers:
#   * p-value : the chance the difference is just luck. SMALL = real.
#   * t-stat  : how big the difference is (further from 0 = bigger).
# We only trust a protein change if BOTH bars below are cleared.
P_THRESHOLD <- 0.05   # p must be under 5% (less than 1-in-20 chance it's luck)
T_THRESHOLD <- 2      # AND the effect size |t| must be over 2 (not tiny)
# Why both? With 44 proteins tested, ~2 will look "significant" by pure
# chance every time. Requiring a real-sized effect too weeds those out.

# --- Knobs for deciding when an ARROW (a proximity relationship) is
# worth drawing. Same idea, but now about "do they sit together?".
INT_PADJ       <- 0.05  # the difference from chance must be significant...
FOLD_THRESHOLD <- 1.5   # ...AND at least 1.5x more (or less) than expected
LOG2_FOLD      <- log2(FOLD_THRESHOLD)  # same threshold, converted to log2
# WHY log2? A ratio of 2 (twice as much) and 1/2 (half as much) are
# equally big changes, but as plain numbers 2 and 0.5 look lopsided.
# Taking log2 makes them symmetric: log2(2) = +1 and log2(1/2) = -1.
# So "+" = more than expected, "-" = fewer, and 0 = exactly expected.
#
# HEADS-UP: we pool counts over ~150 cores, so the counts get huge, and
# with huge counts the significance test says "yes, significant" for
# almost everything. So in practice FOLD_THRESHOLD is the knob that
# actually controls how many arrows you see. Raise it (e.g. 2.0) to
# show fewer, stronger arrows; lower it to show more.

# --- Knobs for the picture itself.
EDGE_CURVE <- 0.45            # how strongly the two arrows between a pair bow
# APART. Higher = the A->B and B->A arcs separate
# more (easier to tell apart); too high and they
# balloon outward. 0.45 keeps them clearly split.
EDGE_MAX_WIDTH <- 6          # thickest an arrow can get. Thinner arrows are
# easier to separate when two run between the same
# pair; raise for more emphasis, lower for clarity.
SHOW_PROTEIN_LABELS <- FALSE  # writing every protein name on every arrow
# turns the graph into unreadable spaghetti,
# so this is OFF. The proteins are still there
# -- just hover an arrow to see them. Turn ON
# only after you've thinned the graph.
N_LABEL_PROTEINS    <- 3      # if labels are ON, how many to print per arrow

LAYOUT_STYLE <- "circle"      # how to arrange the dots:
#  "circle" = evenly spaced ring (tidy, every
#             arrow is easy to see) -- recommended
#  "fr"     = a "spring" layout that pulls
#             strongly-related dots together
#             (meaningful, but can look clumped)

# --- Knobs for the DIVERGING arrow color (blue <- grey -> red).
# The color now shows BOTH direction AND magnitude:
#   dark blue  = way FEWER contacts than expected
#   light blue = a little fewer      (light) grey   = about expected
#   light red  = a little more
#   dark red   = way MORE contacts than expected
COLOR_CAP <- 1.2              # a log2 ratio of this size (~2.3x) or beyond gets
# the FULLEST color. Smaller cap = colors
# saturate sooner (more dramatic); larger =
# more gradual. 2.0 (4x) is a good default.
SHOW_NEUTRAL_EDGES <- FALSE   # FALSE = draw only the significant arrows (clean;
#         you'll mostly see light->dark reds + the
#         blue avoidance pair, little pure grey).
# TRUE  = ALSO draw the near-expected pairs as thin
#         faint grey lines, so the neutral
#         "backbone" is visible. More informative,
#         but busier (all ~100 arrows shown).

# ----------------------------------------------------------------
# STEP 2b: find the protein columns automatically.
# The table has 3 columns per protein: "<name> mean", "<name> t_test",
# "<name> p_value". We don't want to hand-type 44 protein names, so we
# let the code discover them: find every column ending in "mean", strip
# off the " mean", and that leaves the protein name (e.g. "CD45").
# ----------------------------------------------------------------
all_cols   <- colnames(df)                                  # all column names
mean_cols  <- all_cols[grepl("(\\.|[[:space:]])mean$", all_cols)]  # the "... mean" ones
protein_names <- sub("(\\.|[[:space:]])mean$", "", mean_cols)      # chop off " mean"

# Given a protein like "CD45", find its matching t_test/p_value column.
# The separator might be a space ("CD45 t_test") or a dot ("CD45.t_test")
# depending on how the file was saved, so we try both and take whichever exists.
find_matching_col <- function(stub, all_cols, suffix) {
  candidates <- c(paste0(stub, ".", suffix), paste0(stub, " ", suffix))
  hit <- candidates[candidates %in% all_cols]
  if (length(hit) == 0) return(NA_character_)   # found nothing -> flag it
  hit[1]
}
# vapply just runs find_matching_col once for every protein name.
tcols <- unname(vapply(protein_names, find_matching_col, character(1),
                       all_cols = all_cols, suffix = "t_test"))
pcols <- unname(vapply(protein_names, find_matching_col, character(1),
                       all_cols = all_cols, suffix = "p_value"))

# Safety check: if any protein's t_test/p_value column wasn't found, stop
# now with a clear message instead of producing a silently-wrong picture.
if (anyNA(tcols) || anyNA(pcols)) {
  stop("Could not match all t_test/p_value columns. Unmatched: ",
       paste(protein_names[is.na(tcols) | is.na(pcols)], collapse = ", "))
}
n_proteins <- length(protein_names)   # should be 44

# ----------------------------------------------------------------
# STEP 3: for each row, count how many proteins "really changed".
#
# Background on the table's shape: every directed pair in every core has
# TWO rows -- one labeled "Close" (neighbors within 50 microns) and one
# "Far". The t-test already compared Close vs Far, so the SAME t/p
# numbers appear on both rows. That means we only need the Close rows to
# read the test results (using both would double-count). The Close rows
# also conveniently carry the contact counts we need later.
# ----------------------------------------------------------------
close_idx <- df$CellDistance == "Close"   # TRUE/FALSE: which rows are "Close"

# Pull out the t-stats and p-values as plain number grids (matrices):
# one row per Close row, one column per protein.
t_mat <- as.matrix(df[close_idx, tcols])
p_mat <- as.matrix(df[close_idx, pcols])

# sig_mat: a grid of TRUE/FALSE -- TRUE where a protein cleared BOTH bars.
# (!is.na guards against missing values, which we treat as "not significant".)
sig_mat <- !is.na(p_mat) & (p_mat < P_THRESHOLD) & (abs(t_mat) > T_THRESHOLD)

# For each Close row, tally up how many proteins were significant.
close_rows <- df[close_idx, ] %>%
  mutate(
    n_significant = rowSums(sig_mat, na.rm = TRUE),  # count of TRUEs in that row
    has_signal    = n_significant > 0                # did ANY protein change? Y/N
  ) %>%
  # keep only the columns we actually need going forward
  select(TMA_ID, CellType, OtherCellType,
         CellCount, ActualInteractions, ExpectedInteractions,
         n_significant, has_signal)

# ----------------------------------------------------------------
# STEP 3b: how many cells of each type are there? (drives DOT SIZE)
#
# In the table, CellCount is split by distance: the "Close" row counts
# cells that HAVE a neighbor of the other type, the "Far" row counts
# cells that DON'T. Add the two together and you get every cell of that
# type in that core. We compute that per core, then average over cores.
#
# WHY log later? Tumor might have ~4000 cells while a rare immune type
# has ~40. If dot size were the raw count, Tumor would be a beach ball
# and everything else a speck. A logarithm squashes big numbers gently
# (10 -> 1, 100 -> 2, 1000 -> 3, 10000 -> 4), so all dots stay visible
# while still showing "more cells = bigger". (The log is applied in STEP 5.)
# ----------------------------------------------------------------
cells_per_type <- df %>%
  group_by(TMA_ID, CellType, OtherCellType) %>%
  summarise(tot = sum(CellCount, na.rm = TRUE), .groups = "drop") %>%  # Close+Far
  group_by(TMA_ID, CellType) %>%
  summarise(tot = first(tot), .groups = "drop") %>%   # same for every neighbor
  group_by(CellType) %>%
  summarise(mean_cells = mean(tot, na.rm = TRUE), .groups = "drop")    # avg / core

# ----------------------------------------------------------------
# STEP 4: turn ~150 per-core rows for each pair into ONE arrow, and
# decide whether that arrow is worth drawing.
#
# For each directed pair (e.g. "B cells measured next to T cells"):
#   obs = total OBSERVED close contacts, added up over all cores
#   exp = total EXPECTED close contacts (what pure chance predicts,
#         given how many of each cell type there are). Your partner
#         computed this with a Poisson/CSR formula -- "CSR" = Complete
#         Spatial Randomness, i.e. "if cells were sprinkled at random".
#
#   fold = obs / exp. Above 1 = more contacts than chance (attract);
#          below 1 = fewer (avoid). We add 0.5 to top and bottom so we
#          never divide by zero.
#   log2_ratio = log2(fold), the symmetric "+/-" version explained above.
#
#   Poisson test = a stats test that asks "is this observed count
#          surprisingly far from the expected count, or normal wiggle?"
#          It returns a p-value (small = surprising = real signal).
#   BH adjustment = because we run this test ~100 times (once per pair),
#          some will look significant by luck. Benjamini-Hochberg (BH)
#          tightens the p-values to account for all those repeated tests.
#
#   keep = TRUE only if the difference is significant AND big enough
#          (>= 1.5x). That's the filter that removes boring arrows like
#          "we saw 95 where we expected 100".
# ----------------------------------------------------------------
build_directed_edges <- function(rows) {
  rows %>%
    group_by(CellType, OtherCellType) %>%   # one group per directed pair
    summarise(
      obs        = sum(ActualInteractions, na.rm = TRUE),    # observed, pooled
      exp        = sum(ExpectedInteractions, na.rm = TRUE),  # expected, pooled
      mean_obs   = mean(ActualInteractions, na.rm = TRUE),   # per-core avg (for tooltip)
      mean_exp   = mean(ExpectedInteractions, na.rm = TRUE),
      mean_n_sig = mean(n_significant, na.rm = TRUE),        # avg # proteins changed
      frac_signal= mean(has_signal, na.rm = TRUE),           # fraction of cores w/ any change
      n_cores    = n(),                                      # how many cores contributed
      .groups = "drop"
    ) %>%
    mutate(
      log2_ratio = ifelse(exp > 0, log2((obs + 0.5) / (exp + 0.5)), NA_real_),
      fold       = ifelse(exp > 0, (obs + 0.5) / (exp + 0.5), NA_real_),
      # Run the Poisson test once per pair. mapply feeds each pair's obs
      # and exp into the little function. round(o) because the test wants
      # a whole number of events.
      p_int = mapply(function(o, e) {
        if (is.na(e) || e <= 0) return(NA_real_)    # can't test against 0 expected
        poisson.test(round(o), T = 1, r = e, alternative = "two.sided")$p.value
      }, obs, exp)
    ) %>%
    mutate(
      padj = p.adjust(p_int, method = "BH"),   # correct for running many tests
      # THE FILTER: keep an arrow only if it's both real and sizeable.
      keep = !is.na(padj) & padj < INT_PADJ & abs(log2_ratio) >= LOG2_FOLD
    ) %>%
    rename(from = CellType, to = OtherCellType)  # network language: from -> to
}

overall_edges <- build_directed_edges(close_rows)   # the pooled-over-all-cores arrows

# ----------------------------------------------------------------
# STEP 4b: for each arrow, figure out WHICH proteins are doing the work.
# This is the mentor's "later, label the exact proteins" step.
#
# For a given pair, for each protein we compute:
#   frac = the FRACTION of cores where that protein was significant.
#          (Significant in 140 of 150 cores = very reliable; 3 of 150 =
#           probably a fluke.) We rank proteins by this.
#   dir  = the average t-stat in the cores where it WAS significant.
#          Positive means the protein is HIGHER when the cell is near the
#          neighbor (we show an up-arrow); negative means LOWER (down-arrow).
#
# We loop over each unique pair (~100 of them) and build two text labels:
# a short one (top few proteins) for drawing on the arrow, and a full
# ranked one for the hover tooltip.
# ----------------------------------------------------------------
build_protein_labels <- function(row_idx) {
  ct  <- close_rows$CellType[row_idx]       # "from" cell type for each row
  oct <- close_rows$OtherCellType[row_idx]  # "to" cell type for each row
  up_arrow <- "↑"; down_arrow <- "↓"
  
  # the list of unique directed pairs present in these rows
  pairs <- unique(data.frame(from = ct, to = oct, stringsAsFactors = FALSE))
  
  # build the label text for ONE pair (from = f, to = t)
  make_one <- function(f, t) {
    sel <- row_idx[ct == f & oct == t]      # rows (cores) for just this pair
    sm <- sig_mat[sel, , drop = FALSE]      # their significance grid
    tm <- t_mat[sel, , drop = FALSE]        # their t-stats
    sig_num <- ifelse(is.na(sm), 0, sm * 1) # TRUE/FALSE -> 1/0 for counting
    n_sig   <- colSums(sig_num)             # per protein: # of cores significant
    frac    <- n_sig / nrow(sm)             # ...as a fraction of cores
    # average t-stat among ONLY the cores where it was significant -> direction
    tsig    <- ifelse(is.na(tm), 0, tm) * sig_num
    dir_t   <- ifelse(n_sig > 0, colSums(tsig) / n_sig, NA_real_)
    
    keep <- which(frac > 0)                 # proteins that mattered in >=1 core
    if (!length(keep)) return(c(short = "", full = ""))  # none (e.g. self-pairs)
    # sort best-first: most-consistent proteins, ties broken by bigger effect
    keep <- keep[order(-frac[keep], -abs(replace(dir_t[keep], is.na(dir_t[keep]), 0)))]
    arr  <- ifelse(is.na(dir_t[keep]) | dir_t[keep] >= 0, up_arrow, down_arrow)
    
    # short label: e.g. "CD8↑, CD3↑, FOXP3↓"
    short <- paste0(protein_names[keep][seq_len(min(N_LABEL_PROTEINS, length(keep)))],
                    arr[seq_len(min(N_LABEL_PROTEINS, length(keep)))], collapse = ", ")
    # full label with the % of cores: e.g. "CD8↑ (92%), CD3↑ (88%), ..."
    full  <- paste0(protein_names[keep], arr,
                    " (", round(frac[keep] * 100), "%)", collapse = ", ")
    c(short = short, full = full)
  }
  
  # run make_one for every pair, collect into a tidy table
  lab <- t(mapply(make_one, pairs$from, pairs$to))
  tibble(from = pairs$from, to = pairs$to,
         label_short = lab[, "short"], label_full = lab[, "full"])
}

overall_labels <- build_protein_labels(seq_len(nrow(close_rows)))

# ----------------------------------------------------------------
# STEP 5: draw the network.
# This one function takes the arrow table + cell counts + protein labels
# and produces the interactive picture. It's reused for the all-cores
# view and the single-core view.
# ----------------------------------------------------------------

# CONTINUOUS DIVERGING palette. We build two smooth color ramps that
# both START at the same light-grey neutral and get more saturated as the
# ratio grows: one toward red (more than expected), one toward blue (fewer).
# So the SHADE encodes magnitude and the HUE encodes direction.
COL_NEUTRAL <- "#e6e6e6"  # light grey = about expected (the center)
ramp_more <- colorRampPalette(c(COL_NEUTRAL, "#fb6a4a", "#de2d26", "#a50f15"))(100) # -> vivid red
ramp_less <- colorRampPalette(c(COL_NEUTRAL, "#6baed6", "#2171b5", "#08519c"))(100) # -> steel/deep blue

# turn a log2 ratio into a color: pick red or blue ramp by sign, then pick
# how deep into the ramp by how big the ratio is (capped at COLOR_CAP).
ratio_to_color <- function(log2r) {
  vapply(log2r, function(x) {
    if (is.na(x)) return(COL_NEUTRAL)
    mag <- min(abs(x) / COLOR_CAP, 1)          # 0 (neutral) .. 1 (fully saturated)
    idx <- max(1, round(mag * 99) + 1)         # position along the 100-step ramp
    if (x >= 0) ramp_more[idx] else ramp_less[idx]
  }, character(1))
}

make_network_plot <- function(edge_summary, cells_tbl, protein_labels,
                              title_text, filter_edges = TRUE) {
  # every cell type that appears anywhere becomes a dot
  nodes_vec <- sort(unique(c(edge_summary$from, edge_summary$to)))
  
  # which arrows to draw:
  #   - normally, only the significant ones (keep == TRUE)
  #   - if SHOW_NEUTRAL_EDGES is on, draw ALL of them (the neutral ones will
  #     be thin + light grey so they read as a faint background)
  disp <- if (filter_edges && !SHOW_NEUTRAL_EDGES) filter(edge_summary, keep) else edge_summary
  disp <- disp %>% left_join(protein_labels, by = c("from", "to")) %>%
    mutate(label_short = ifelse(is.na(label_short), "", label_short),
           label_full  = ifelse(is.na(label_full),  "", label_full))
  
  # ---- WHERE to place each dot ----------------------------------
  # We build a lightweight network object (igraph) just to compute
  # positions. Self-loops don't help positioning, so we leave them out
  # of the layout math (they still get DRAWN later).
  layout_edges <- disp %>% filter(from != to)
  g_ig <- graph_from_data_frame(
    layout_edges %>% transmute(from, to, w = pmax(mean_n_sig, 1e-6)),
    directed = TRUE,
    vertices = data.frame(name = nodes_vec)
  )
  set.seed(42)  # fixes randomness so the picture looks the same every run
  if (LAYOUT_STYLE == "circle") {
    # place the dots evenly around a circle (tidy and readable)
    layout_raw <- layout_in_circle(g_ig, order = order(nodes_vec))
    layout_scaled <- layout_raw * 500   # spread them out to fill the canvas
  } else {
    # "spring" layout: strongly-related dots get pulled together
    layout_raw <- layout_with_fr(g_ig, weights = E(g_ig)$w)
    layout_scaled <- layout_raw * 600
  }
  
  # ---- the DOTS (nodes) -----------------------------------------
  cell_lookup <- setNames(cells_tbl$mean_cells, cells_tbl$CellType)  # name -> count
  nodes <- tibble(id = nodes_vec) %>%
    mutate(
      label = id,
      raw_cells = as.numeric(cell_lookup[id]),
      value = log10(ifelse(is.na(raw_cells), 0, raw_cells) + 1),  # <-- the log squash
      x = layout_scaled[match(id, nodes_vec), 1],   # position from the layout
      y = layout_scaled[match(id, nodes_vec), 2],
      title = paste0("<b>", id, "</b><br>Mean cells/core: ",
                     round(raw_cells), " (log-scaled node)")   # hover text
    )
  
  # ---- the ARROWS (edges) ---------------------------------------
  edges <- disp %>%
    mutate(
      color   = ratio_to_color(log2_ratio),   # diverging blue<-grey->red (shade = magnitude)
      # thickness = # proteins changed for real arrows; neutral backbone stays thin
      value   = ifelse(!is.na(keep) & keep, mean_n_sig, 0.4),
      arrows  = "to",                          # draw an arrowhead at the "to" end
      label   = if (SHOW_PROTEIN_LABELS) label_short else NA_character_,
      # curve the arrows; A->B and B->A bow opposite ways so both are visible
      smooth.enabled   = TRUE,
      smooth.type      = ifelse(from < to, "curvedCW", "curvedCCW"),
      smooth.roundness = EDGE_CURVE,
      direction = ifelse(log2_ratio >= 0, "MORE than expected (attract)",
                         "FEWER than expected (avoid)"),
      # the hover tooltip: the full story for this one arrow
      title = paste0(
        "<b>", from, " → ", to, "</b><br>",
        "Interactions: ", direction, "<br>",
        "Observed/Expected fold: ", round(fold, 2),
        "  (log2 = ", round(log2_ratio, 2), ")<br>",
        "Mean observed/core: ", round(mean_obs, 1),
        " vs expected ", round(mean_exp, 1), "<br>",
        "Poisson p(adj): ", signif(padj, 2), "<br>",
        "Avg # significant proteins/core: ", round(mean_n_sig, 1),
        " / ", n_proteins, "<br>",
        "Cores with >=1 shifted protein: ", round(frac_signal * 100), "%<br>",
        "Top proteins (", up_or_down_legend(), "): ",
        ifelse(label_full == "", "none (self-pair has no near/far contrast)",
               label_full)
      )
    ) %>%
    select(from, to, value, color, arrows, label,
           smooth.enabled, smooth.type, smooth.roundness, title)
  
  # ---- hand the dots + arrows to visNetwork to render ------------
  visNetwork(nodes, edges, main = title_text,
             width = "100%", height = "850px") %>%
    visNodes(shape = "dot", scaling = list(min = 15, max = 55),   # dot size range
             color = list(background = "#a6cee3", border = "#3182bd"),
             font = list(size = 18, color = "#222222",
                         strokeWidth = 3, strokeColor = "#ffffff")) %>%
    visEdges(scaling = list(min = 1, max = EDGE_MAX_WIDTH),                    # arrow width range
             font = list(size = 10, color = "#333333",
                         strokeWidth = 4, strokeColor = "#ffffff", align = "top"),
             arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5))) %>%
    # clicking/hovering a dot highlights just its arrows
    visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
               nodesIdSelection = TRUE) %>%
    visPhysics(enabled = FALSE) %>%   # freeze the layout so it doesn't jiggle
    visInteraction(tooltipDelay = 100, dragNodes = TRUE,
                   dragView = TRUE, zoomView = TRUE) %>%
    visLayout(randomSeed = 42)
}

# tiny helper used in the tooltip text above
up_or_down_legend <- function() paste0("↑ higher / ↓ lower near neighbor")

# prints a plain-text legend + how many arrows survived the filter
print_legend <- function(edge_summary) {
  n_disp <- sum(edge_summary$keep, na.rm = TRUE)
  n_all  <- nrow(edge_summary)
  cat("\nLegend:\n",
      "  Arrow points FROM the cell type whose proteins were tested\n",
      "  TO the neighbor it was tested against.\n",
      "  Color : diverging scale. RED = more than expected (attract),\n",
      "          BLUE = fewer than expected (avoid), GREY = about\n",
      "          expected. DARKER shade = bigger difference from expected.\n",
      "  Width : average # of significant proteins (near vs far)\n",
      "  Label : top shifted proteins (up = higher, down = lower near)\n",
      "  Node  : log-scaled number of cells of that type\n",
      "  Curve : reciprocal arrows (A->B, B->A) bow opposite ways\n",
      "  Loop  : cell type aggregating with itself\n",
      "  Shown : ", n_disp, " of ", n_all, " directed pairs pass the\n",
      "          significance + fold-change filter.\n\n", sep = "")
}

# ================================================================
# STEP 6: build the OVERALL picture (averaged over every core) and save it.
# ================================================================
overall_plot <- make_network_plot(
  overall_edges, cells_per_type, overall_labels,
  paste0("Cell-Type Interaction Network – All ",
         n_distinct(df$TMA_ID), " Cores (obs vs expected)")
)
saveWidget(overall_plot, "overall_interaction_network.html", selfcontained = TRUE)
print("Saved: overall_interaction_network.html")
print_legend(overall_edges)

# ================================================================
# STEP 7: build a picture for ONE core only. Change core_id to inspect a
# different tissue sample (each patient/sample can look different).
# ================================================================
core_id    <- unique(df$TMA_ID)[1]              # <-- change [1] to [2], [3]... etc
core_mask  <- close_rows$TMA_ID == core_id      # which rows belong to this core
core_rows  <- close_rows[core_mask, ]
core_edges <- build_directed_edges(core_rows)   # arrows for just this core
core_labels<- build_protein_labels(which(core_mask))
# this core's own cell counts (so dot sizes match this specific sample)
core_cells <- df %>%
  filter(TMA_ID == core_id) %>%
  group_by(CellType, OtherCellType) %>%
  summarise(tot = sum(CellCount, na.rm = TRUE), .groups = "drop") %>%
  group_by(CellType) %>%
  summarise(mean_cells = first(tot), .groups = "drop")
core_plot  <- make_network_plot(
  core_edges, core_cells, core_labels,
  paste0("Cell-Type Interaction Network – Core: ", core_id)
)
saveWidget(core_plot, "core_interaction_network.html", selfcontained = TRUE)
print(paste("Saved: core_interaction_network.html  (core:", core_id, ")"))

# ================================================================
# STEP 8: print a plain-text ranked table of the arrows we actually drew,
# strongest first, with the proteins driving each. Handy for a quick
# sanity check or to paste into a slide/message without opening the html.
# ================================================================
overall_edges %>%
  filter(keep) %>%                                       # only the drawn arrows
  left_join(overall_labels, by = c("from", "to")) %>%    # attach protein names
  arrange(desc(abs(log2_ratio))) %>%                     # biggest change first
  mutate(direction = ifelse(log2_ratio >= 0, "MORE", "FEWER")) %>%
  select(from, to, direction, fold, log2_ratio, padj,
         mean_n_sig, top_proteins = label_short) %>%
  print(n = 200)